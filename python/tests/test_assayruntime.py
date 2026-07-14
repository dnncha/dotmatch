from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from assaycode import cli
from dotmatch.assayruntime import run_compiled_plan
from dotmatch.assayscript import AssayScriptError, load_and_compile, write_compiled_plan


def _library(path: Path, rows: list[tuple[str, str]]) -> Path:
    path.write_text(
        "target_id\ttarget_seq\n" + "".join(f"{name}\t{sequence}\n" for name, sequence in rows),
        encoding="utf-8",
    )
    return path


def _fastq(path: Path, records: list[tuple[str, str]]) -> Path:
    path.write_text(
        "".join(f"@{read_id}\n{sequence}\n+\n{'I' * len(sequence)}\n" for read_id, sequence in records),
        encoding="utf-8",
    )
    return path


def _compiled_fixture(tmp_path: Path, *, combinations: str | None = None) -> Path:
    _library(tmp_path / "samples.tsv", [("s1", "AAAA"), ("s2", "CCCC")])
    _library(tmp_path / "guides.tsv", [("g1", "ACGT"), ("g2", "ACGA")])
    constraint = ""
    if combinations is not None:
        (tmp_path / "pairs.tsv").write_text(combinations, encoding="utf-8")
        constraint = '\n[constraints]\nallowed_combinations = "pairs.tsv"\n'
    spec = tmp_path / "assay.toml"
    spec.write_text(
        """
schema_version = 2
name = "joint-screen"
assay_type = "crispr"

[[segments]]
name = "sample"
read = "I1"
library = "samples.tsv"
start = 0
length = 4
metric = "hamming"
k = 0

[[segments]]
name = "guide"
read = "R2"
library = "guides.tsv"
start = 0
length = 4
metric = "hamming"
k = 1
""".lstrip()
        + constraint,
        encoding="utf-8",
    )
    return write_compiled_plan(load_and_compile(spec), tmp_path / "assay.plan.json")


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_allowed_tuple_resolves_ambiguous_segment_calls(tmp_path: Path) -> None:
    plan = _compiled_fixture(tmp_path, combinations="sample\tguide\ns1\tg1\ns2\tg2\n")
    i1 = _fastq(tmp_path / "I1.fastq", [("r1/1", "AAAA"), ("r2/1", "CCCC")])
    r2 = _fastq(tmp_path / "R2.fastq", [("r1/2", "ACGT"), ("r2/2", "ACGA")])

    result = run_compiled_plan(
        plan, {"I1": i1, "R2": r2}, tmp_path / "out", accept_findings=True
    )

    rows = _rows(result.assignments)
    assert [row["status"] for row in rows] == ["unique", "unique"]
    assert [row["reason"] for row in rows] == ["allowed_combination_rescue"] * 2
    assert rows[0]["guide_status"] == "ambiguous"
    assert rows[0]["guide_candidates"] == "g1;g2"
    assert rows[0]["guide_target"] == "g1"
    assert rows[1]["guide_target"] == "g2"
    assert _rows(result.counts) == [
        {"combination": "guide=g1|sample=s1", "count": "1"},
        {"combination": "guide=g2|sample=s2", "count": "1"},
    ]


def test_disallowed_tuple_is_not_forced(tmp_path: Path) -> None:
    plan = _compiled_fixture(tmp_path, combinations="sample\tguide\ns1\tg1\n")
    i1 = _fastq(tmp_path / "I1.fastq", [("r1", "CCCC")])
    r2 = _fastq(tmp_path / "R2.fastq", [("r1", "ACGT")])

    result = run_compiled_plan(
        plan, {"I1": i1, "R2": r2}, tmp_path / "out", accept_findings=True
    )

    row = _rows(result.assignments)[0]
    assert row["status"] == "none"
    assert row["reason"] == "disallowed_combination"
    assert _rows(result.counts) == []


def test_runtime_preserves_ambiguity_without_constraints(tmp_path: Path) -> None:
    plan = _compiled_fixture(tmp_path)
    i1 = _fastq(tmp_path / "I1.fastq", [("r1", "AAAA")])
    r2 = _fastq(tmp_path / "R2.fastq", [("r1", "ACGT")])

    result = run_compiled_plan(
        plan, {"I1": i1, "R2": r2}, tmp_path / "out", accept_findings=True
    )

    assert _rows(result.assignments)[0]["status"] == "ambiguous"
    summary = json.loads(result.summary.read_text(encoding="utf-8"))
    assert summary["status_counts"]["ambiguous"] == 1
    assert summary["rates"]["ambiguous"] == 1.0


def test_runtime_fails_before_outputs_when_fastqs_desynchronize(tmp_path: Path) -> None:
    plan = _compiled_fixture(tmp_path)
    i1 = _fastq(tmp_path / "I1.fastq", [("r1", "AAAA"), ("r2", "AAAA")])
    r2 = _fastq(tmp_path / "R2.fastq", [("r1", "ACGT")])
    output = tmp_path / "out"

    with pytest.raises(AssayScriptError, match="record counts differ"):
        run_compiled_plan(plan, {"I1": i1, "R2": r2}, output, accept_findings=True)

    assert not output.exists()


def test_runtime_rejects_changed_library_fingerprint(tmp_path: Path) -> None:
    plan = _compiled_fixture(tmp_path)
    (tmp_path / "guides.tsv").write_text("target_id\ttarget_seq\ng1\tTTTT\n", encoding="utf-8")

    with pytest.raises(AssayScriptError, match="fingerprint changed"):
        run_compiled_plan(plan, {}, tmp_path / "out")


def test_runtime_rejects_unknown_allowed_target(tmp_path: Path) -> None:
    plan = _compiled_fixture(tmp_path, combinations="sample\tguide\ns1\tmissing\n")
    i1 = _fastq(tmp_path / "I1.fastq", [("r1", "AAAA")])
    r2 = _fastq(tmp_path / "R2.fastq", [("r1", "ACGT")])

    with pytest.raises(AssayScriptError, match="unknown guide target missing"):
        run_compiled_plan(
            plan, {"I1": i1, "R2": r2}, tmp_path / "out", accept_findings=True
        )


def test_anchor_jitter_and_reverse_complement_extraction(tmp_path: Path) -> None:
    _library(tmp_path / "targets.tsv", [("t1", "ACGT")])
    spec = tmp_path / "anchor.toml"
    spec.write_text(
        """
schema_version = 2
name = "anchor"
assay_type = "panel"

[[segments]]
name = "target"
read = "R1"
library = "targets.tsv"
anchor = "GGGG"
after_anchor = 0
jitter = 1
length = 4
orientation = "reverse_complement"
metric = "hamming"
k = 0
""".lstrip(),
        encoding="utf-8",
    )
    plan = write_compiled_plan(load_and_compile(spec), tmp_path / "anchor.plan.json")
    # Reverse complement is GGGGTACGT; jitter finds ACGT one base after the anchor.
    r1 = _fastq(tmp_path / "R1.fastq", [("r1", "ACGTACCCC")])

    result = run_compiled_plan(plan, {"R1": r1}, tmp_path / "out", accept_findings=True)

    row = _rows(result.assignments)[0]
    assert row["status"] == "unique"
    assert row["target_target"] == "t1"
    assert row["target_start"] == "5"
    assert row["target_orientation"] == "reverse_complement"


def test_max_reads_is_recorded_and_events_are_watch_compatible(tmp_path: Path) -> None:
    plan = _compiled_fixture(tmp_path, combinations="sample\tguide\ns1\tg1\n")
    i1 = _fastq(tmp_path / "I1.fastq", [("r1", "AAAA"), ("r2", "AAAA")])
    r2 = _fastq(tmp_path / "R2.fastq", [("r1", "ACGT"), ("r2", "ACGT")])

    result = run_compiled_plan(
        plan,
        {"I1": i1, "R2": r2},
        tmp_path / "out",
        max_reads=1,
        accept_findings=True,
    )

    assert result.total_reads == 1
    event = json.loads(result.events.read_text(encoding="utf-8"))
    assert event == {
        "read_id": "r1",
        "status": "unique",
        "target": "guide=g1|sample=s1",
    }


def test_assaycode_execute_cli_runs_compiled_plan(tmp_path: Path, capsys) -> None:
    plan = _compiled_fixture(tmp_path, combinations="sample\tguide\ns1\tg1\n")
    i1 = _fastq(tmp_path / "I1.fastq", [("r1", "AAAA")])
    r2 = _fastq(tmp_path / "R2.fastq", [("r1", "ACGT")])

    assert cli.main(
        [
            "execute",
            str(plan),
            "--i1",
            str(i1),
            "--r2",
            str(r2),
            "--out",
            str(tmp_path / "run"),
            "--accept-findings",
        ]
    ) == 0

    summary = json.loads(capsys.readouterr().out)
    assert summary["status"] == "experimental"
    assert summary["total_reads"] == 1
    assert summary["status_counts"] == {"unique": 1}
    assert Path(summary["assignments"]).is_file()


def test_runtime_requires_explicit_acceptance_of_compiler_findings(tmp_path: Path) -> None:
    plan = _compiled_fixture(tmp_path)

    with pytest.raises(AssayScriptError, match="explicitly accept_findings"):
        run_compiled_plan(plan, {}, tmp_path / "out")


def test_runtime_rejects_non_integer_max_reads(tmp_path: Path) -> None:
    plan = _compiled_fixture(tmp_path)

    with pytest.raises(ValueError, match="positive integer"):
        run_compiled_plan(plan, {}, tmp_path / "out", max_reads=1.5, accept_findings=True)
