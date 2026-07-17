"""Light tests for tl and pure parsers (no heavy optional deps required for basic paths)."""

import gzip

import pytest

from dotmatch.multiqc import (
    DOTMATCH_SEARCH_PATTERNS,
    _multiqc_file_path,
    parse_sample_qc_tsv,
    parse_crispr_qc_summary_tsv,
    parse_assay_manifest_summary_tsv,
    parse_summary_json,
    parse_panel_summary_json,
    parse_dotmatch_artifacts,
)

# Use the real fixtures that the workflow checks also validate
ROOT = __import__("pathlib").Path(__file__).resolve().parents[2]


def test_pure_parsers_on_fixtures():
    sample = parse_sample_qc_tsv(ROOT / "examples/workflows/multiqc/data/sample_qc.tsv")
    assert "plasmid" in sample or len(sample) > 0
    sample_row = next(iter(sample.values()))
    assert "assignment_rate" in sample_row
    assert isinstance(sample_row["assignment_rate"], float)
    assert isinstance(sample_row["total_reads"], int)

    crispr = parse_crispr_qc_summary_tsv(ROOT / "examples/workflows/multiqc/data/crispr_qc.summary.tsv")
    row = next(iter(crispr.values()))
    assert "qc_status" in row
    assert "zero_count_fraction" in row

    manifest = parse_assay_manifest_summary_tsv(ROOT / "examples/workflows/multiqc/data/assay_manifest.summary.tsv")
    row = next(iter(manifest.values()))
    assert "status" in row
    assert "primary_report" in row

    summary = parse_summary_json(ROOT / "examples/workflows/fixtures/assay_out/summary.json")
    row = next(iter(summary.values()))
    assert "assigned_unique" in row
    assert "ambiguous" in row

    panel = parse_panel_summary_json(ROOT / "examples/workflows/multiqc/data/panel_summary.json")
    row = next(iter(panel.values()))
    assert "status" in row
    assert "minimum_hamming_distance" in row


def test_native_multiqc_parser_contract() -> None:
    assert DOTMATCH_SEARCH_PATTERNS["dotmatch/sample_qc"]["fn"] == "*sample_qc.tsv"
    assert DOTMATCH_SEARCH_PATTERNS["dotmatch/crispr_qc"]["fn"] == "*crispr_qc.summary.tsv"
    assert DOTMATCH_SEARCH_PATTERNS["dotmatch/top_unmatched"]["fn"] == "*top_unmatched.tsv"

    parsed = parse_dotmatch_artifacts(
        [
            ROOT / "examples/workflows/multiqc/data/sample_qc.tsv",
            ROOT / "examples/workflows/multiqc/data/crispr_qc.summary.tsv",
            ROOT / "examples/workflows/multiqc/data/assay_manifest.summary.tsv",
            ROOT / "examples/workflows/multiqc/data/panel_summary.json",
        ]
    )
    assert {"sample_qc", "crispr_qc", "assay_manifest", "panel_summary"} <= set(parsed)
    assert isinstance(next(iter(parsed["crispr_qc"].values()))["coverage_fraction"], float)


def test_multiqc_file_path_accepts_common_record_shapes() -> None:
    assert _multiqc_file_path({"root": "/tmp/run", "fn": "sample_qc.tsv"}).as_posix() == "/tmp/run/sample_qc.tsv"
    assert _multiqc_file_path({"path": "/tmp/run/crispr_qc.summary.tsv"}).as_posix() == "/tmp/run/crispr_qc.summary.tsv"


def test_tl_importable():
    # tl requires anndata at runtime for most calls, but the module must import
    import dotmatch.tl as tl
    assert hasattr(tl, "assign_features")
    assert hasattr(tl, "feature_counts")
    assert hasattr(tl, "crispr_guide_assignment")


def test_assign_features_requires_anndata(monkeypatch):
    import dotmatch.tl as tl
    # Force the guard
    monkeypatch.setattr(tl, "_HAS_ANNDATA", False)
    with pytest.raises(ImportError, match="anndata"):
        tl.assign_features(None, library=[], seq_col="x")


def test_stream_assign_loads_targets_and_summarizes_without_loading_reads(tmp_path):
    import dotmatch

    targets = tmp_path / "targets.csv"
    targets.write_text("target_id,target_seq\na,ACGT\nb,TTTT\n", encoding="utf-8")
    reads = tmp_path / "reads.fastq.gz"
    with gzip.open(reads, "wt", encoding="utf-8") as fh:
        fh.write(
            "@exact\nACGT\n+\nIIII\n"
            "@corrected\nACGA\n+\nIIII\n"
            "@none\nGGGG\n+\nIIII\n"
            "@short\nAC\n+\nII\n"
        )

    rows = list(dotmatch.stream_assign(reads, targets, target_length=4, k=1))

    assert [row.status_name for row in rows] == ["unique", "unique", "none", "invalid"]
    assert rows[0].target_name == "a"
    assert rows[1].best_distance == 1
    summary = dotmatch.assignment_summary(rows)
    assert summary["total_reads"] == 4
    assert summary["assigned_unique"] == 2
    assert summary["assigned_exact"] == 1
    assert summary["assigned_corrected"] == 1
    assert summary["unmatched"] == 1
    assert summary["invalid"] == 1
    assert summary["assignment_rate"] == 0.5


def test_stream_assign_metric_hamming_does_not_rescue_indels(tmp_path):
    import dotmatch

    reads = tmp_path / "reads.fastq"
    reads.write_text("@indel\nACGTT\n+\nIIIII\n", encoding="utf-8")

    hamming = list(dotmatch.stream_assign(reads, [("guide_1", "ACGT")], target_length=5, k=1, metric="hamming"))
    levenshtein = list(dotmatch.stream_assign(reads, [("guide_1", "ACGT")], target_length=5, k=1))

    assert hamming[0].status_name == "none"
    assert levenshtein[0].status_name == "unique"


def test_write_assignments_tsv_returns_summary(tmp_path):
    import dotmatch

    reads = tmp_path / "reads.fastq"
    reads.write_text("@r0\nACGT\n+\nIIII\n@r1\nCCCC\n+\nIIII\n", encoding="utf-8")
    out = tmp_path / "assignments.tsv"

    summary = dotmatch.write_assignments_tsv(
        dotmatch.stream_assign(reads, [("guide_1", "ACGT")], target_length=4, k=0),
        out,
    )

    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "read_id\tobserved_seq\ttarget_id\ttarget_seq\tdistance\tstatus\tmatch_count\tsecond_best_distance"
    assert "r0\tACGT\tguide_1\tACGT\t0\tunique\t1\t-1" in lines
    assert "r1\tCCCC\t\t\t-1\tnone\t0\t-1" in lines
    assert summary["assigned_unique"] == 1
    assert summary["unmatched"] == 1
