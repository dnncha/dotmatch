from __future__ import annotations

import json
import socket
from pathlib import Path

import pytest

try:
    from jsonschema import Draft202012Validator
except ModuleNotFoundError:
    Draft202012Validator = None  # type: ignore[assignment,misc]

from dotmatch.agent_tools import (
    _allowed_fixes,
    _artifact,
    _export_skill,
    _write_candidate_spec,
    invoke_tool,
    list_tools,
)


def test_contract_defines_exact_six_tools() -> None:
    contract = list_tools()
    schema_path = Path(__file__).resolve().parents[2] / "agent-tools.schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    if Draft202012Validator is not None:
        Draft202012Validator.check_schema(schema)
        Draft202012Validator(schema).validate(contract)
    else:
        assert schema["$schema"] == "https://json-schema.org/draft/2020-12/schema"
        assert "$defs" in schema
    assert contract["schema_version"] == "1.0"
    from dotmatch import __version__
    assert contract["generated_for_version"] == __version__
    assert {item["name"] for item in contract["tools"]} == {
        "discover",
        "prepare_assay",
        "inspect_assay",
        "run_assay",
        "review_assay",
        "handoff_assay",
    }


def test_discover_returns_stable_local_envelope() -> None:
    result = invoke_tool("discover", {"intent": "crispr-guide-counting"})
    assert result["schema_version"] == "1.0"
    assert result["tool_contract_version"] == "1.0"
    assert result["status"] == "passed"
    assert result["exit_code"] == 0
    assert result["result"]["local_only"] is True
    assert result["result"]["network_requests"] == "none"
    assert result["spec"] == {"revision": 0, "path": "", "sha256": ""}
    schema_path = Path(__file__).resolve().parents[2] / "agent-tools.schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    if Draft202012Validator is not None:
        envelope_schema = {"$ref": "#/$defs/envelope", "$defs": schema["$defs"]}
        Draft202012Validator(envelope_schema).validate(result)
    else:
        assert set(schema["$defs"]["envelope"]["required"]) <= set(result)


def test_unknown_and_shell_shaped_fields_are_rejected(tmp_path: Path) -> None:
    marker = tmp_path / "must-not-exist"
    result = invoke_tool("discover", {"command": f"touch {marker}"})
    assert result["status"] == "invalid_input"
    assert result["exit_code"] == 2
    assert not marker.exists()


def test_python_api_keeps_a_schema_valid_envelope_for_an_invalid_tool_type() -> None:
    result = invoke_tool(7, {})  # type: ignore[arg-type]
    assert result["tool"] == ""
    assert result["status"] == "invalid_input"
    assert result["exit_code"] == 2
    assert "tool must be a string" in result["findings"][0]["message"]


def test_contract_types_are_not_coerced(tmp_path: Path) -> None:
    targets = tmp_path / "targets.tsv"
    reads_dir = tmp_path / "reads"
    targets.write_text("guide_a\tACGT\n", encoding="utf-8")
    reads_dir.mkdir()
    (reads_dir / "sample.fastq").write_text("@r0\nACGT\n+\nIIII\n", encoding="utf-8")

    string_boolean = invoke_tool(
        "prepare_assay",
        {
            "intent": "crispr-guide-counting",
            "targets": str(targets),
            "reads_dir": str(reads_dir),
            "output_dir": str(tmp_path / "output-a"),
            "link_reads": "false",
            "minimum_free_bytes": 0,
        },
    )
    assert string_boolean["status"] == "invalid_input"
    assert "link_reads must be a boolean" in string_boolean["findings"][0]["message"]

    string_integer = invoke_tool(
        "prepare_assay",
        {
            "intent": "crispr-guide-counting",
            "targets": str(targets),
            "reads_dir": str(reads_dir),
            "output_dir": str(tmp_path / "output-b"),
            "threads": "4",
            "minimum_free_bytes": 0,
        },
    )
    assert string_integer["status"] == "invalid_input"
    assert "threads must be an integer" in string_integer["findings"][0]["message"]


def test_prepare_multithreaded_count_uses_aggregate_diagnostics(tmp_path: Path) -> None:
    targets = tmp_path / "targets.tsv"
    reads_dir = tmp_path / "reads"
    targets.write_text("guide_a\tACGT\n", encoding="utf-8")
    reads_dir.mkdir()
    (reads_dir / "sample.fastq").write_text("@r0\nACGT\n+\nIIII\n", encoding="utf-8")

    result = invoke_tool(
        "prepare_assay",
        {
            "intent": "crispr-guide-counting",
            "targets": str(targets),
            "reads_dir": str(reads_dir),
            "output_dir": str(tmp_path / "output"),
            "threads": 4,
            "minimum_free_bytes": 0,
        },
    )

    assert result["status"] in {"passed", "needs_review"}
    assert any(item["finding_id"] == "aggregate_diagnostics_only" for item in result["findings"])
    spec = (tmp_path / "output" / "assay.toml").read_text(encoding="utf-8")
    assert "threads = 4" in spec
    assert "assignments = false" in spec
    assert "ambiguous = false" in spec
    assert "unmatched = false" in spec
    assert not any("--assignments" in item for item in result["result"]["plan"])


def test_structured_paths_are_not_shell_interpreted(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    targets = tmp_path / "targets;touch pwned"
    reads_dir = tmp_path / "reads"
    targets.write_text("guide_a\tACGT\n", encoding="utf-8")
    reads_dir.mkdir()
    (reads_dir / "sample.fastq").write_text("@r0\nACGT\n+\nIIII\n", encoding="utf-8")

    result = invoke_tool(
        "prepare_assay",
        {
            "intent": "crispr-guide-counting",
            "targets": str(targets),
            "reads_dir": str(reads_dir),
            "output_dir": str(tmp_path / "output"),
            "minimum_free_bytes": 0,
        },
    )
    assert result["status"] in {"passed", "needs_review"}
    assert not (tmp_path / "pwned").exists()


def test_prepare_refuses_non_empty_output(tmp_path: Path) -> None:
    targets = tmp_path / "targets.tsv"
    reads_dir = tmp_path / "reads"
    output = tmp_path / "output"
    targets.write_text("guide_a\tACGT\n", encoding="utf-8")
    reads_dir.mkdir()
    (reads_dir / "sample.fastq").write_text("@r0\nACGT\n+\nIIII\n", encoding="utf-8")
    output.mkdir()
    (output / "owned.txt").write_text("keep\n", encoding="utf-8")

    result = invoke_tool(
        "prepare_assay",
        {
            "intent": "crispr-guide-counting",
            "targets": str(targets),
            "reads_dir": str(reads_dir),
            "output_dir": str(output),
        },
    )
    assert result["status"] == "invalid_input"
    assert (output / "owned.txt").read_text(encoding="utf-8") == "keep\n"


def test_prepare_refuses_symlink_output(tmp_path: Path) -> None:
    targets = tmp_path / "targets.tsv"
    reads_dir = tmp_path / "reads"
    actual_output = tmp_path / "actual-output"
    linked_output = tmp_path / "linked-output"
    targets.write_text("guide_a\tACGT\n", encoding="utf-8")
    reads_dir.mkdir()
    (reads_dir / "sample.fastq").write_text("@r0\nACGT\n+\nIIII\n", encoding="utf-8")
    actual_output.mkdir()
    linked_output.symlink_to(actual_output, target_is_directory=True)

    result = invoke_tool(
        "prepare_assay",
        {
            "intent": "crispr-guide-counting",
            "targets": str(targets),
            "reads_dir": str(reads_dir),
            "output_dir": str(linked_output),
        },
    )
    assert result["status"] == "invalid_input"
    assert "symlink" in result["findings"][0]["message"]


def test_insufficient_disk_is_a_block_not_invalid_input(tmp_path: Path) -> None:
    spec = tmp_path / "assay.toml"
    spec.write_text("schema_version = 1\n", encoding="utf-8")
    result = invoke_tool(
        "inspect_assay",
        {"spec": str(spec), "minimum_free_bytes": 10**30},
    )
    assert result["status"] == "blocked"
    assert "insufficient disk" in result["findings"][0]["message"]


def test_prepare_rejects_a_truncated_fastq(tmp_path: Path) -> None:
    targets = tmp_path / "targets.tsv"
    reads_dir = tmp_path / "reads"
    targets.write_text("guide_a\tACGT\n", encoding="utf-8")
    reads_dir.mkdir()
    (reads_dir / "truncated.fastq").write_text("@r0\nACGT\n+\n", encoding="utf-8")
    result = invoke_tool(
        "prepare_assay",
        {
            "intent": "crispr-guide-counting",
            "targets": str(targets),
            "reads_dir": str(reads_dir),
            "output_dir": str(tmp_path / "output"),
            "minimum_free_bytes": 0,
        },
    )
    assert result["status"] == "invalid_input"
    assert "truncated FASTQ" in result["findings"][0]["message"]


def test_inspect_blocks_an_unsafe_target_library(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = Path(__file__).resolve().parents[2]
    native = root / "dotmatch"
    if not native.is_file():
        pytest.skip("candidate native CLI is not built")
    monkeypatch.setenv("DOTMATCH_NATIVE_CLI", str(native))

    targets = tmp_path / "unsafe-targets.tsv"
    reads = tmp_path / "reads.fastq"
    targets.write_text("guide_a\tAAAAAA\tGENEA\nguide_b\tAAAAAT\tGENEB\n", encoding="utf-8")
    reads.write_text("@r0\nAAAAAA\n+\nIIIIII\n", encoding="utf-8")
    spec = tmp_path / "assay.toml"
    spec.write_text(
        f'''schema_version = 1
status = "ready"
mode = "count"
assay_type = "crispr"
targets = "{targets}"

[[samples]]
id = "sample"
fastq = "{reads}"

[run]
out_dir = "{tmp_path / 'assay_out'}"

[extract]
start = 0
length = 6

[assignment]
k = 1
metric = "hamming"
ambiguity_policy = "radius"
ambiguous = "discard"
''',
        encoding="utf-8",
    )

    result = invoke_tool("inspect_assay", {"spec": str(spec), "minimum_free_bytes": 0})
    assert result["status"] == "blocked"
    assert any(item["finding_id"] == "unsafe_targets" for item in result["findings"])


def test_forbidden_scientific_remediations_are_never_returned_as_allowed() -> None:
    allowed, forbidden = _allowed_fixes(
        {
            "assay_fixes": [
                {"section": "extract", "key": "start", "suggested_value": "2"},
                {"section": "targets", "key": "library", "suggested_value": "changed.tsv"},
                {"section": "assignment", "key": "ambiguous", "suggested_value": "count"},
                {"section": "reliability", "key": "profile", "suggested_value": "exploratory"},
                {"section": "backend", "key": "mode", "suggested_value": "gpu-metal-experimental"},
            ]
        }
    )
    assert [(item["section"], item["key"]) for item in allowed] == [("extract", "start")]
    assert {(item["section"], item["key"]) for item in forbidden} == {
        ("targets", "library"),
        ("assignment", "ambiguous"),
        ("reliability", "profile"),
        ("backend", "mode"),
    }


def test_export_skill_requires_empty_directory_and_has_metadata(tmp_path: Path) -> None:
    target = tmp_path / "skill"
    result = _export_skill(str(target))
    assert result["status"] == "passed"
    assert (target / "SKILL.md").is_file()
    assert (target / "agents" / "openai.yaml").is_file()
    assert sorted(path.name for path in (target / "references").iterdir()) == [
        "crispr.md",
        "evidence-policy.md",
        "perturb-seq.md",
    ]
    with pytest.raises(ValueError, match="must be empty"):
        _export_skill(str(target))


def test_candidate_revision_is_numbered_and_does_not_edit_original(tmp_path: Path) -> None:
    targets = tmp_path / "targets.tsv"
    reads = tmp_path / "reads.fastq"
    targets.write_text("guide_a\tACGT\n", encoding="utf-8")
    reads.write_text("@r0\nNNACGT\n+\nIIIIII\n", encoding="utf-8")
    original = tmp_path / "assay.toml"
    original_text = f'''schema_version = 1
status = "ready"
mode = "count"
assay_type = "crispr"
targets = "{targets}"

[[samples]]
id = "sample"
fastq = "{reads}"

[run]
out_dir = "{tmp_path / 'assay_out'}"

[extract]
start = 0
length = 4

[assignment]
k = 1
metric = "hamming"
ambiguous = "discard"
'''
    original.write_text(original_text, encoding="utf-8")
    candidate = _write_candidate_spec(
        original,
        [{"section": "extract", "key": "start", "suggested_value": "2"}],
        1,
    )
    assert candidate.name == "assay.agent-r1.toml"
    assert "start = 2" in candidate.read_text(encoding="utf-8")
    assert "assay_out.agent-r1" in candidate.read_text(encoding="utf-8")
    assert original.read_text(encoding="utf-8") == original_text
    with pytest.raises(ValueError, match="refusing to overwrite"):
        _write_candidate_spec(
            original,
            [{"section": "extract", "key": "start", "suggested_value": "2"}],
            1,
        )


def test_run_stops_after_three_numbered_revisions(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import dotmatch.agent_tools as tools_module

    targets = tmp_path / "targets.tsv"
    reads = tmp_path / "reads.fastq"
    targets.write_text("guide_a\tACGT\n", encoding="utf-8")
    reads.write_text("@r0\nACGT\n+\nIIII\n", encoding="utf-8")
    spec = tmp_path / "assay.toml"
    spec.write_text(
        f'''schema_version = 1
status = "ready"
mode = "count"
assay_type = "crispr"
targets = "{targets}"

[[samples]]
id = "sample"
fastq = "{reads}"

[run]
out_dir = "{tmp_path / 'assay_out'}"

[extract]
start = 0
length = 4

[assignment]
k = 1
metric = "hamming"
ambiguous = "discard"
''',
        encoding="utf-8",
    )

    monkeypatch.setattr(tools_module, "_captured_assay", lambda _command: (2, "", ""))

    def failed_summary(path: Path) -> tuple[dict[str, object], dict[str, Path]]:
        revision = tools_module._revision_from_path(path)
        return (
            {
                "overall_status": "failed",
                "findings": [
                    {
                        "finding_id": "assignment_rate_below_min",
                        "severity": "error",
                        "message": "simulated bounded revision",
                    }
                ],
                "assay_fixes": [
                    {
                        "section": "extract",
                        "key": "start",
                        "suggested_value": str(revision + 1),
                    }
                ],
            },
            {},
        )

    monkeypatch.setattr(tools_module, "_load_reliability", failed_summary)
    result = invoke_tool(
        "run_assay",
        {"spec": str(spec), "max_revisions": 3, "minimum_free_bytes": 0},
    )

    assert result["status"] == "blocked"
    assert result["spec"]["revision"] == 3
    assert [item["revision"] for item in result["result"]["revision_history"]] == [0, 1, 2, 3]
    assert any(item["finding_id"] == "revision_limit_exhausted" for item in result["findings"])


def test_contract_json_is_plain_json() -> None:
    # Guard against accidentally placing comments or host-specific paths in the public contract.
    encoded = json.dumps(list_tools(), sort_keys=True)
    assert "/" + "Users/" not in encoded
    assert "shell_commands_accepted" in encoded


def test_artifact_hashing_is_stable(tmp_path: Path) -> None:
    path = tmp_path / "artifact.tsv"
    path.write_text("a\tb\n1\t2\n", encoding="utf-8")
    first = _artifact(path, role="fixture")
    second = _artifact(path, role="fixture")
    assert first["sha256"] == second["sha256"]
    assert len(first["sha256"]) == 64
    assert first["bytes"] == path.stat().st_size


def test_discover_makes_no_network_connection(monkeypatch: pytest.MonkeyPatch) -> None:
    def refuse_connect(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("network connection attempted")

    monkeypatch.setattr(socket.socket, "connect", refuse_connect)
    result = invoke_tool("discover", {})
    assert result["status"] == "passed"


def test_interrupted_execution_has_stable_status(monkeypatch: pytest.MonkeyPatch) -> None:
    import dotmatch.agent_tools as tools_module

    def interrupt(_data: dict[str, object], _envelope: dict[str, object]) -> None:
        raise KeyboardInterrupt

    monkeypatch.setitem(tools_module._TOOL_HANDLERS, "discover", interrupt)
    result = invoke_tool("discover", {})
    assert result["status"] == "interrupted"
    assert result["exit_code"] == 130
    assert result["findings"][0]["finding_id"] == "interrupted_execution"


def test_crispr_agent_prepare_run_review_and_handoff(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = Path(__file__).resolve().parents[2]
    native = root / "dotmatch"
    if not native.is_file():
        pytest.skip("candidate native CLI is not built")
    monkeypatch.setenv("DOTMATCH_NATIVE_CLI", str(native))

    def refuse_connect(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("ordinary agent workflow attempted an outbound connection")

    monkeypatch.setattr(socket.socket, "connect", refuse_connect)

    sequences = [
        "AAAAAA",
        "CCCCCC",
        "GGGGGG",
        "TTTTTT",
        "ACGTAC",
        "CATGCA",
        "GACTGA",
        "TGACGT",
        "AGCTTC",
        "CTAGAG",
    ]
    targets = tmp_path / "guides.tsv"
    targets.write_text(
        "guide_id\tsequence\tgene\n"
        + "".join(f"guide_{index}\t{sequence}\tGENE{index}\n" for index, sequence in enumerate(sequences)),
        encoding="utf-8",
    )
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    reads = reads_dir / "sample.fastq"
    reads.write_text(
        "".join(
            f"@read_{index}_{repeat}\nNN{sequence}AA\n+\nIIIIIIIIII\n"
            for index, sequence in enumerate(sequences)
            for repeat in range(10)
        ),
        encoding="utf-8",
    )
    project = tmp_path / "agent-project"

    prepared = invoke_tool(
        "prepare_assay",
        {
            "intent": "crispr-guide-counting",
            "targets": str(targets),
            "reads_dir": str(reads_dir),
            "output_dir": str(project),
            "minimum_free_bytes": 0,
        },
    )
    assert prepared["status"] == "passed"
    spec = prepared["spec"]["path"]

    inspected = invoke_tool("inspect_assay", {"spec": spec, "minimum_free_bytes": 0})
    assert inspected["status"] == "passed"
    completed = invoke_tool("run_assay", {"spec": spec, "minimum_free_bytes": 0})
    assert completed["status"] == "passed"
    assert completed["result"]["revision_history"][-1]["revision"] == 0
    assert completed["spec"]["sha256"]

    reviewed = invoke_tool("review_assay", {"spec": spec})
    assert reviewed["status"] == "passed"
    assert reviewed["result"]["ambiguity_policy"] == "discard"
    assert any(item["role"] == "counts" and item["sha256"] for item in reviewed["artifacts"])

    handoff_dir = tmp_path / "handoff"
    handed_off = invoke_tool("handoff_assay", {"spec": spec, "output_dir": str(handoff_dir)})
    assert handed_off["status"] == "passed"
    assert handed_off["result"]["raw_data_included"] is False
    assert not list(handoff_dir.rglob("*.fastq"))
    assert (handoff_dir / "handoff_manifest.json").is_file()
    assert (handoff_dir / "SHA256SUMS").is_file()

    repeated = invoke_tool("run_assay", {"spec": spec, "minimum_free_bytes": 0})
    assert repeated["status"] == "invalid_input"
    assert "refusing to overwrite completed assay outputs" in repeated["findings"][0]["message"]


def test_agent_revises_a_wrong_reverse_orientation_without_editing_original(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = Path(__file__).resolve().parents[2]
    native = root / "dotmatch"
    if not native.is_file():
        pytest.skip("candidate native CLI is not built")
    monkeypatch.setenv("DOTMATCH_NATIVE_CLI", str(native))

    sequences = [
        "AACCGGTA",
        "CCGTAAGC",
        "GGTACCAT",
        "TTGCAACG",
        "ACAGTCCA",
        "CAGGATTC",
        "GTCACTAG",
        "TACGAGGT",
        "AGTCGCAA",
        "GACATGTA",
    ]
    complement = str.maketrans("ACGTN", "TGCAN")

    def reverse_complement(sequence: str) -> str:
        return sequence.translate(complement)[::-1]

    targets = tmp_path / "guides.tsv"
    targets.write_text(
        "guide_id\tsequence\tgene\n"
        + "".join(f"guide_{index}\t{sequence}\tGENE{index}\n" for index, sequence in enumerate(sequences)),
        encoding="utf-8",
    )
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    (reads_dir / "sample.fastq").write_text(
        "".join(
            f"@read_{index}_{repeat}\n{reverse_complement('NN' + sequence + 'AA')}\n+\nIIIIIIIIIIII\n"
            for index, sequence in enumerate(sequences)
            for repeat in range(10)
        ),
        encoding="utf-8",
    )
    project = tmp_path / "reverse-project"
    prepared = invoke_tool(
        "prepare_assay",
        {
            "intent": "crispr-guide-counting",
            "targets": str(targets),
            "reads_dir": str(reads_dir),
            "output_dir": str(project),
            "minimum_free_bytes": 0,
        },
    )
    assert prepared["status"] == "passed"
    original = Path(prepared["spec"]["path"])
    inferred = original.read_text(encoding="utf-8")
    assert 'orientation = "reverse_complement"' in inferred

    # Simulate a reviewed-but-wrong starting hypothesis. The agent must keep
    # this file immutable once execution begins and write assay.agent-r1.toml.
    original.write_text(
        inferred.replace('orientation = "reverse_complement"', 'orientation = "forward"'),
        encoding="utf-8",
    )
    wrong_hash = original.read_bytes()
    completed = invoke_tool(
        "run_assay",
        {"spec": str(original), "minimum_free_bytes": 0, "max_revisions": 3},
    )

    assert completed["status"] == "passed"
    assert completed["spec"]["revision"] == 1
    candidate = Path(completed["spec"]["path"])
    assert candidate.name == "assay.agent-r1.toml"
    assert 'orientation = "reverse_complement"' in candidate.read_text(encoding="utf-8")
    assert original.read_bytes() == wrong_hash
    assert [item["revision"] for item in completed["result"]["revision_history"]] == [0, 1]
    oriented = candidate.parent / "assay_out.agent-r1" / "oriented_inputs" / "sample.reverse-complement.fastq"
    assert oriented.is_file()
    handed_off = invoke_tool(
        "handoff_assay",
        {
            "spec": str(candidate),
            "output_dir": str(tmp_path / "reverse-handoff"),
            "minimum_free_bytes": 0,
        },
    )
    assert handed_off["status"] == "passed"
    assert not list((tmp_path / "reverse-handoff").rglob("*.fastq"))
