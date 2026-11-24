from __future__ import annotations

import json
from pathlib import Path

import pytest

from dotmatch.triage import triage_run


def _write_targets(path: Path, rows: list[tuple[str, str]]) -> Path:
    path.write_text(
        "id\tsequence\n" + "".join(f"{tid}\t{seq}\n" for tid, seq in rows),
        encoding="utf-8",
    )
    return path


def _write_fastq(path: Path, seqs: list[str]) -> Path:
    lines = [f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n" for i, seq in enumerate(seqs)]
    path.write_text("".join(lines), encoding="utf-8")
    return path


def test_healthy_run_produces_healthy_verdict_and_json_artifact(tmp_path: Path) -> None:
    targets = _write_targets(
        tmp_path / "targets.tsv",
        [
            ("g1", "AAAAAAAA"),
            ("g2", "CCCCCCCC"),
            ("g3", "GGGGGGGG"),
            ("g4", "TTTTTTTT"),
        ],
    )
    reads = _write_fastq(
        tmp_path / "reads.fastq",
        ["NNAAAAAAAA", "NNCCCCCCCC", "NNGGGGGGGG", "NNTTTTTTTT"] * 5,
    )
    out_dir = tmp_path / "triage_out"

    report = triage_run(
        targets=targets,
        reads=reads,
        target_start=2,
        target_length=8,
        out_dir=out_dir,
    )

    assert report.verdict == "healthy"
    assert report.read_count == 20
    assert report.target_count == 4
    assert report.exact_assignment_rate == pytest.approx(1.0)
    assert report.radius_k1_assignment_rate == pytest.approx(1.0)
    assert report.best_k1_assignment_rate == pytest.approx(1.0)
    assert report.changed_reads == 0
    assert report.findings == ()

    report_json = out_dir / "triage_report.json"
    assert report_json.is_file()
    payload = json.loads(report_json.read_text(encoding="utf-8"))
    assert payload["schema_version"] == "dotmatch.triage.v1"
    assert payload["verdict"] == "healthy"
    assert payload["read_count"] == 20


def test_offset_shift_is_diagnosed_as_likely_misconfigured(tmp_path: Path) -> None:
    targets = _write_targets(
        tmp_path / "targets.tsv",
        [
            ("g1", "ACGTTAGC"),
            ("g2", "TGCATCGA"),
            ("g3", "CTAGCTAG"),
            ("g4", "GGAATTCC"),
        ],
    )
    reads = _write_fastq(
        tmp_path / "reads.fastq",
        ["NNNNNACGTTAGCAA", "NNNNNTGCATCGAAA", "NNNNNCTAGCTAGAA", "NNNNNGGAATTCCAA"] * 5,
    )

    report = triage_run(
        targets=targets,
        reads=reads,
        target_start=2,
        target_length=8,
        out_dir=tmp_path / "triage_offset",
    )

    assert report.verdict == "likely_misconfigured"
    codes = {finding.code: finding for finding in report.findings}
    assert "offset_shift_detected" in codes
    assert codes["offset_shift_detected"].severity == "critical"
    assert codes["offset_shift_detected"].evidence["suggested_target_start"] == 5


def test_reverse_complement_is_diagnosed_as_likely_misconfigured(tmp_path: Path) -> None:
    targets = _write_targets(
        tmp_path / "targets.tsv",
        [
            ("g1", "AACGTACC"),
            ("g2", "TTGCATGG"),
            ("g3", "CCATGCTA"),
            ("g4", "GGCTTAAG"),
        ],
    )
    # Reverse complements of g1..g4 placed at offset 2
    revcomps = ["GGTACGTT", "CCATGCAA", "TAGCATGG", "CTTAAGCC"]
    reads = _write_fastq(
        tmp_path / "reads.fastq",
        [f"NN{rc}NN" for rc in revcomps] * 5,
    )

    report = triage_run(
        targets=targets,
        reads=reads,
        target_start=2,
        target_length=8,
        out_dir=tmp_path / "triage_rc",
    )

    assert report.verdict == "likely_misconfigured"
    codes = {finding.code: finding for finding in report.findings}
    assert "reverse_complement_detected" in codes
    assert codes["reverse_complement_detected"].severity == "critical"


def test_duplicate_targets_flagged_critical(tmp_path: Path) -> None:
    targets = _write_targets(
        tmp_path / "targets.tsv",
        [
            ("g1", "AAAAAAAA"),
            ("g2", "AAAAAAAA"),
            ("g3", "CCCCCCCC"),
        ],
    )
    reads = _write_fastq(tmp_path / "reads.fastq", ["NNCCCCCCCC"] * 10)

    report = triage_run(
        targets=targets,
        reads=reads,
        target_start=2,
        target_length=8,
        out_dir=tmp_path / "triage_dup",
    )

    assert report.verdict == "likely_misconfigured"
    codes = {finding.code: finding for finding in report.findings}
    assert "duplicate_targets" in codes
    assert codes["duplicate_targets"].severity == "critical"
    assert codes["duplicate_targets"].evidence["duplicate_target_count"] == 2


def test_k1_collision_and_policy_sensitivity_flag_review_required(tmp_path: Path) -> None:
    targets = _write_targets(
        tmp_path / "targets.tsv",
        [
            ("g1", "AAAAAAAA"),
            ("g2", "AAAAAAAC"),
            ("g3", "GGGGGGGG"),
        ],
    )
    reads = _write_fastq(
        tmp_path / "reads.fastq",
        ["NNAAAAAAAA", "NNAAAAAAAC", "NNAAAAAAAT", "NNGGGGGGGG"] * 5,
    )

    report = triage_run(
        targets=targets,
        reads=reads,
        target_start=2,
        target_length=8,
        out_dir=tmp_path / "triage_policy",
    )

    assert report.verdict == "review_required"
    codes = {finding.code: finding for finding in report.findings}
    assert "library_collision_k1" in codes
    assert codes["library_collision_k1"].severity == "warning"
    assert "policy_sensitive_counts" in codes
    assert codes["policy_sensitive_counts"].severity == "warning"


def test_high_invalid_window_rate_flagged_critical(tmp_path: Path) -> None:
    targets = _write_targets(
        tmp_path / "targets.tsv",
        [("g1", "AAAAAAAA"), ("g2", "CCCCCCCC")],
    )
    reads = _write_fastq(tmp_path / "reads.fastq", ["ACGT"] * 10)

    report = triage_run(
        targets=targets,
        reads=reads,
        target_start=2,
        target_length=8,
        out_dir=tmp_path / "triage_invalid",
    )

    assert report.verdict == "likely_misconfigured"
    codes = {finding.code: finding for finding in report.findings}
    assert "high_invalid_window_rate" in codes
    assert codes["high_invalid_window_rate"].severity == "critical"


def test_refuses_non_empty_output_directory(tmp_path: Path) -> None:
    targets = _write_targets(tmp_path / "targets.tsv", [("g1", "AAAAAAAA")])
    reads = _write_fastq(tmp_path / "reads.fastq", ["NNAAAAAAAA"])
    out_dir = tmp_path / "existing"
    out_dir.mkdir()
    (out_dir / "keep.txt").write_text("existing", encoding="utf-8")

    with pytest.raises(ValueError, match="empty directory"):
        triage_run(
            targets=targets,
            reads=reads,
            target_start=2,
            target_length=8,
            out_dir=out_dir,
        )
