"""Light tests for tl and pure parsers (no heavy optional deps required for basic paths)."""

import pytest

from dotmatch.multiqc import (
    parse_sample_qc_tsv,
    parse_crispr_qc_summary_tsv,
    parse_assay_manifest_summary_tsv,
    parse_summary_json,
    parse_panel_summary_json,
)

# Use the real fixtures that the workflow checks also validate
ROOT = __import__("pathlib").Path(__file__).resolve().parents[2]


def test_pure_parsers_on_fixtures():
    sample = parse_sample_qc_tsv(ROOT / "examples/workflows/multiqc/data/sample_qc.tsv")
    assert "plasmid" in sample or len(sample) > 0
    assert "assignment_rate" in next(iter(sample.values()))

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
