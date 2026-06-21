"""Light tests for tl and pure parsers (no heavy optional deps required for basic paths)."""

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
