from __future__ import annotations

import json
from pathlib import Path

from dotmatch.segments import ReadGrammar, SegmentSpec, decode_segments


def _write_table(path: Path, rows: list[tuple[str, str]]) -> Path:
    path.write_text(
        "id\tsequence\n" + "".join(f"{tid}\t{seq}\n" for tid, seq in rows),
        encoding="utf-8",
    )
    return path


def _write_fastq(path: Path, seqs: list[str]) -> Path:
    lines = [f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n" for i, seq in enumerate(seqs)]
    path.write_text("".join(lines), encoding="utf-8")
    return path


def test_perturb_seq_cell_umi_guide_deduplicates_molecules(tmp_path: Path) -> None:
    cells = _write_table(
        tmp_path / "cells.tsv",
        [("C1", "AAAA"), ("C2", "CCCC")],
    )
    guides = _write_table(
        tmp_path / "guides.tsv",
        [("G1", "GGGGGG"), ("G2", "TTTTTT")],
    )
    # Layout: [0:4]=cell, [4:8]=umi, [8:14]=guide
    reads = _write_fastq(
        tmp_path / "reads.fastq",
        [
            "AAAAACGTGGGGGG",  # C1 + UMI ACGT + G1 (PCR dup 1)
            "AAAAACGTGGGGGG",  # C1 + UMI ACGT + G1 (PCR dup 2)
            "AAAATGCAGGGGGg",  # C1 + UMI TGCA + G1 (new molecule, lowercase base)
            "CCCCACGTTTTTTT",  # C2 + UMI ACGT + G2 (same UMI in different cell -> distinct molecule)
        ],
    )
    grammar = ReadGrammar(
        segments=(
            SegmentSpec(name="cell", kind="target", start=0, length=4, targets=cells, k=0),
            SegmentSpec(name="umi", kind="umi", start=4, length=4, scope=("cell",)),
            SegmentSpec(name="guide", kind="target", start=8, length=6, targets=guides, k=0),
        )
    )
    out_dir = tmp_path / "decoded"

    report = decode_segments(grammar=grammar, reads=reads, out_dir=out_dir)

    assert report.total_reads == 4
    assert report.joint_unique_reads == 4
    assert report.joint_ambiguous_reads == 0
    assert report.joint_unmatched_reads == 0
    assert report.joint_invalid_reads == 0
    assert report.umi_dedup_enabled is True
    assert report.unique_molecules == 3
    assert report.umi_conflict_molecules == 0
    assert report.umi_conflict_reads == 0

    counts_lines = (out_dir / "segment_counts.tsv").read_text(encoding="utf-8").splitlines()
    assert counts_lines[0] == "cell\tguide\tread_count\tumi_count"
    rows = {tuple(line.split("\t")[:2]): tuple(int(x) for x in line.split("\t")[2:]) for line in counts_lines[1:]}
    assert rows == {
        ("C1", "G1"): (3, 2),
        ("C2", "G2"): (1, 1),
    }

    summary = json.loads((out_dir / "segments_summary.json").read_text(encoding="utf-8"))
    assert summary["schema_version"] == "dotmatch.segments.v1"
    assert summary["unique_molecules"] == 3


def test_umi_conflict_within_scope_is_explicitly_excluded_and_logged(tmp_path: Path) -> None:
    cells = _write_table(tmp_path / "cells.tsv", [("C1", "AAAA")])
    guides = _write_table(tmp_path / "guides.tsv", [("G1", "GGGGGG"), ("G2", "TTTTTT")])
    reads = _write_fastq(
        tmp_path / "reads.fastq",
        [
            "AAAAACGTGGGGGG",  # C1 + ACGT -> G1
            "AAAAACGTTTTTTT",  # C1 + ACGT -> G2 (conflict on same cell+UMI!)
            "AAAATGCAGGGGGG",  # C1 + TGCA -> G1 (clean molecule)
        ],
    )
    grammar = ReadGrammar(
        segments=(
            SegmentSpec(name="cell", kind="target", start=0, length=4, targets=cells, k=0),
            SegmentSpec(name="umi", kind="umi", start=4, length=4, scope=("cell",)),
            SegmentSpec(name="guide", kind="target", start=8, length=6, targets=guides, k=0),
        )
    )
    out_dir = tmp_path / "conflicts"

    report = decode_segments(grammar=grammar, reads=reads, out_dir=out_dir)

    assert report.joint_unique_reads == 3
    assert report.unique_molecules == 1
    assert report.umi_conflict_molecules == 1
    assert report.umi_conflict_reads == 2

    conflicts = (out_dir / "umi_conflicts.tsv").read_text(encoding="utf-8").splitlines()
    assert conflicts[0] == "cell\tumi\tconflicting_targets\tread_count"
    assert len(conflicts) == 2
    assert conflicts[1].startswith("C1\tACGT\t")
    assert conflicts[1].endswith("\t2")


def test_paired_end_del_decoding_and_explicit_joint_outcomes(tmp_path: Path) -> None:
    cycle1 = _write_table(tmp_path / "cycle1.tsv", [("A1", "AAAA"), ("A2", "AACC")])
    cycle2 = _write_table(tmp_path / "cycle2.tsv", [("B1", "GGGG"), ("B2", "TTTT")])
    # R1: [0:4]=cycle1, [4:8]=umi; R2: [0:4]=cycle2
    r1 = _write_fastq(
        tmp_path / "r1.fastq",
        [
            "AAAACCCC",  # unique A1 (d=0 to A1, d=2 to A2) + valid UMI
            "AAACCCCC",  # ambiguous at k=1 (d=1 to A1 and d=1 to A2)
            "TTTTCCCC",  # unmatched cycle1
            "AAAACCNN",  # invalid UMI (contains N)
        ],
    )
    r2 = _write_fastq(
        tmp_path / "r2.fastq",
        [
            "GGGG",  # unique B1
            "GGGG",  # unique B1
            "GGGG",  # unique B1
            "GGGG",  # unique B1
        ],
    )
    grammar = ReadGrammar(
        segments=(
            SegmentSpec(name="cycle1", kind="target", start=0, length=4, read=1, targets=cycle1, k=1),
            SegmentSpec(name="umi", kind="umi", start=4, length=4, read=1),
            SegmentSpec(name="cycle2", kind="target", start=0, length=4, read=2, targets=cycle2, k=0),
        )
    )
    out_dir = tmp_path / "del_out"

    report = decode_segments(grammar=grammar, reads=r1, reads_r2=r2, out_dir=out_dir)

    assert report.total_reads == 4
    assert report.joint_unique_reads == 1
    assert report.joint_ambiguous_reads == 1
    assert report.joint_unmatched_reads == 1
    assert report.joint_invalid_reads == 1
    assert report.unique_molecules == 1


def test_paired_end_rejects_mismatched_record_counts_and_invalid_grammar(tmp_path: Path) -> None:
    import pytest

    cycle1 = _write_table(tmp_path / "c1.tsv", [("A1", "AAAA")])
    cycle2 = _write_table(tmp_path / "c2.tsv", [("B1", "GGGG")])
    r1 = _write_fastq(tmp_path / "r1.fastq", ["AAAA", "AAAA"])
    r2 = _write_fastq(tmp_path / "r2.fastq", ["GGGG"])
    grammar = ReadGrammar(
        segments=(
            SegmentSpec(name="c1", kind="target", start=0, length=4, read=1, targets=cycle1),
            SegmentSpec(name="c2", kind="target", start=0, length=4, read=2, targets=cycle2),
        )
    )

    with pytest.raises(ValueError, match="reads_r2 is required"):
        decode_segments(grammar=grammar, reads=r1, out_dir=tmp_path / "missing_r2")

    with pytest.raises(ValueError, match="different record counts"):
        decode_segments(grammar=grammar, reads=r1, reads_r2=r2, out_dir=tmp_path / "mismatch")

    with pytest.raises(ValueError, match="unknown target segment"):
        ReadGrammar(
            segments=(
                SegmentSpec(name="c1", kind="target", start=0, length=4, targets=cycle1),
                SegmentSpec(name="umi", kind="umi", start=4, length=4, scope=("nonexistent",)),
            )
        )
