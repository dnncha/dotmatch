from __future__ import annotations

import json
from pathlib import Path

import pytest

from assaycode import cli
from dotmatch.assayscript import AssayScriptError, load_and_compile, write_compiled_plan


def _write_library(path: Path, rows: list[tuple[str, str]]) -> Path:
    path.write_text(
        "target_id\ttarget_seq\n"
        + "".join(f"{target_id}\t{sequence}\n" for target_id, sequence in rows),
        encoding="utf-8",
    )
    return path


def _write_spec(tmp_path: Path, *, duplicate: bool = False) -> Path:
    guides = _write_library(
        tmp_path / "guides.tsv",
        [("g1", "ACGT"), ("g2", "ACGT" if duplicate else "TTTT")],
    )
    samples = _write_library(
        tmp_path / "samples.tsv",
        [("s1", "AACCGGTT"), ("s2", "TTGGCCAA")],
    )
    combinations = tmp_path / "pairs.tsv"
    combinations.write_text("sample\tguide\ns1\tg1\ns2\tg2\n", encoding="utf-8")
    spec = tmp_path / "assay.toml"
    spec.write_text(
        f"""
schema_version = 2
name = "dual-guide-assay"
assay_type = "crispr"

[[segments]]
name = "sample"
read = "I1"
library = "{samples.name}"
start = 0
length = 8
metric = "hamming"
k = 1

[[segments]]
name = "guide"
read = "R2"
library = "{guides.name}"
anchor = "GTTT"
after_anchor = 1
jitter = 2
length = 4
metric = "levenshtein"
k = 1
orientation = "auto"

[constraints]
allowed_combinations = "{combinations.name}"
""".lstrip(),
        encoding="utf-8",
    )
    return spec


def test_compiler_selects_segment_specific_strategies(tmp_path: Path) -> None:
    plan = load_and_compile(_write_spec(tmp_path))

    by_name = {segment.name: segment for segment in plan.segments}
    assert by_name["sample"].strategy == "packed_hamming_neighborhood"
    assert by_name["guide"].strategy.startswith(
        "orientation_dispatch+jitter_scan+anchor_scan+"
    )
    assert plan.execution_order == ["sample", "guide"]
    assert plan.allowed_combinations_sha256
    assert all(len(segment.library_sha256) == 64 for segment in plan.segments)


def test_compiler_records_reviewable_runtime_findings(tmp_path: Path) -> None:
    plan = load_and_compile(_write_spec(tmp_path))

    assert any("orientation auto" in finding for finding in plan.findings)
    assert any("jitter search" in finding for finding in plan.findings)
    assert any("allowed combinations contain 2 tuples" in finding for finding in plan.findings)


def test_compiler_marks_duplicate_library_unsafe(tmp_path: Path) -> None:
    plan = load_and_compile(_write_spec(tmp_path, duplicate=True))
    guide = next(segment for segment in plan.segments if segment.name == "guide")

    assert guide.safety_status == "unsafe"
    assert any("duplicate target" in finding for finding in guide.findings)


def test_compiler_rejects_ambiguous_extraction_definition(tmp_path: Path) -> None:
    spec = _write_spec(tmp_path)
    spec.write_text(
        spec.read_text(encoding="utf-8").replace(
            'anchor = "GTTT"',
            'start = 4\nanchor = "GTTT"',
        ),
        encoding="utf-8",
    )

    with pytest.raises(AssayScriptError, match="exactly one of start or anchor"):
        load_and_compile(spec)


def test_compiler_rejects_unknown_combination_segment(tmp_path: Path) -> None:
    spec = _write_spec(tmp_path)
    combinations = tmp_path / "pairs.tsv"
    combinations.write_text("sample\tunknown\ns1\tx\n", encoding="utf-8")

    with pytest.raises(AssayScriptError, match="unknown segments"):
        load_and_compile(spec)


def test_compiled_plan_is_portable_json_with_fingerprints(tmp_path: Path) -> None:
    plan = load_and_compile(_write_spec(tmp_path))
    output = write_compiled_plan(plan, tmp_path / "build" / "assay.plan.json")
    data = json.loads(output.read_text(encoding="utf-8"))

    assert data["compiler_schema_version"] == 1
    assert data["source_schema_version"] == 2
    assert len(data["source_sha256"]) == 64
    assert data["name"] == "dual-guide-assay"


def test_assaycode_compile_and_inspect_commands(tmp_path: Path, capsys) -> None:
    spec = _write_spec(tmp_path)
    output = tmp_path / "assay.plan.json"

    assert cli.main(["compile", str(spec), "--out", str(output)]) == 0
    compile_summary = json.loads(capsys.readouterr().out)
    assert compile_summary["segments"] == 2
    assert output.exists()

    assert cli.main(["inspect", str(output)]) == 0
    inspect_summary = json.loads(capsys.readouterr().out)
    assert inspect_summary["name"] == "dual-guide-assay"
    assert {segment["name"] for segment in inspect_summary["segments"]} == {
        "sample",
        "guide",
    }


def test_assaycode_compile_rejects_v1_without_breaking_dotmatch(tmp_path: Path, capsys) -> None:
    spec = tmp_path / "old.toml"
    spec.write_text('schema_version = 1\n', encoding="utf-8")

    assert cli.main(["compile", str(spec), "--out", str(tmp_path / "plan.json")]) == 2
    assert "schema_version must be 2" in capsys.readouterr().err
