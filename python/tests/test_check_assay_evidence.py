import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_assay_evidence.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_assay_evidence", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _valid_manifest() -> dict:
    return {
        "schema_version": 1,
        "assays": [
            {
                "id": "crispr_guide_counting",
                "label": "CRISPR guide counting",
                "status": "supported",
                "raw_artifacts": [
                    "benchmarks/raw/public_crispr_repeated.csv",
                    "benchmarks/raw/public_crispr_edlib_validation.csv",
                ],
                "reports": ["docs/benchmarks/public_crispr/README.md"],
                "gates": ["make public-crispr-evidence-gate"],
                "claim_boundary": "Public MAGeCK/Yusa rows only.",
                "commands": [
                    "make bench-public-crispr-repeated",
                    "make public-crispr-evidence-gate",
                ],
                "comparator_semantics": "MAGeCK exact-count agreement plus Edlib assignment validation.",
                "validation": "Gate requires Edlib mismatches equal zero.",
            },
            {
                "id": "inline_barcode",
                "label": "Inline barcode demultiplexing",
                "status": "gated",
                "raw_artifacts": ["benchmarks/raw/barcode_demux.csv"],
                "reports": ["docs/benchmarks/barcode_demux/README.md"],
                "gates": ["make barcode-comparison-gate"],
                "claim_boundary": "Fixture rows cover workflow smoke evidence only.",
                "next_public_evidence": "Add public barcode-sheet rows with comparator agreement.",
                "commands": ["make bench-barcode-demux"],
                "comparator_semantics": "Fixture rows record deterministic demux execution.",
                "validation": "Smoke gate checks deterministic demux execution.",
            },
            {
                "id": "perturb_seq",
                "label": "Perturb-seq guide or feature assignment",
                "status": "planned",
                "raw_artifacts": [],
                "reports": [],
                "gates": [],
                "claim_boundary": "No public perturb-seq evidence claim yet.",
                "next_public_evidence": "Add a public Perturb-seq guide or feature-barcode FASTQ fixture and oracle.",
            },
            {
                "id": "feature_barcode",
                "label": "Feature barcode assignment",
                "status": "planned",
                "raw_artifacts": [],
                "reports": [],
                "gates": [],
                "claim_boundary": "No public feature-barcode evidence claim yet.",
                "next_public_evidence": "Add a public cell-hashing or CITE-seq barcode fixture and comparator semantics.",
            },
            {
                "id": "amplicon_panel",
                "label": "Amplicon or panel target assignment",
                "status": "planned",
                "raw_artifacts": [],
                "reports": [],
                "gates": [],
                "claim_boundary": "No public amplicon/panel evidence claim yet.",
                "next_public_evidence": "Add a public panel-style FASTQ fixture and validation oracle.",
            },
            {
                "id": "oligo_adapter",
                "label": "Oligo or adapter target assignment",
                "status": "planned",
                "raw_artifacts": [],
                "reports": [],
                "gates": [],
                "claim_boundary": "No public oligo/adapter evidence claim yet.",
                "next_public_evidence": "Add a public oligo or adapter FASTQ fixture and validation oracle.",
            },
        ],
    }


def _write_assay_repo(root: Path, manifest=None) -> None:
    files = {
        "Makefile": (
            "public-crispr-evidence-gate:\n\ttrue\n"
            "barcode-comparison-gate:\n\ttrue\n"
            "bench-public-crispr-repeated:\n\ttrue\n"
            "bench-barcode-demux:\n\ttrue\n"
        ),
        "docs/assay-evidence.json": json.dumps(manifest or _valid_manifest(), indent=2) + "\n",
        "benchmarks/raw/public_crispr_repeated.csv": (
            "tool,workflow,command,exit_code\n"
            "dotmatch_hamming_k1,public_crispr_yusa,dotmatch count --targets guides.tsv,0\n"
        ),
        "benchmarks/raw/public_crispr_edlib_validation.csv": "sample,mismatches,checked_reads\nplasmid,0,100\n",
        "benchmarks/raw/barcode_demux.csv": (
            "tool,workflow,command,exit_code\n"
            "dotmatch_demux,real_public_inline_barcode,dotmatch demux --barcodes barcodes.tsv,0\n"
        ),
        "docs/benchmarks/public_crispr/README.md": "# Public CRISPR\n",
        "docs/benchmarks/barcode_demux/README.md": "# Barcode Demux\n",
    }
    for path, text in files.items():
        full = root / path
        full.parent.mkdir(parents=True, exist_ok=True)
        full.write_text(text, encoding="utf-8")


def test_assay_evidence_accepts_manifest_with_required_lanes(tmp_path):
    checker = _load_checker()
    _write_assay_repo(tmp_path)

    result = checker.audit(tmp_path)

    assert result.failures == []
    assert any("required assay lanes" in item for item in result.passed)


def test_assay_evidence_rejects_missing_required_lane(tmp_path):
    checker = _load_checker()
    manifest = _valid_manifest()
    manifest["assays"] = [assay for assay in manifest["assays"] if assay["id"] != "perturb_seq"]
    _write_assay_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("missing required assay lane: perturb_seq" in failure for failure in result.failures)


def test_assay_evidence_rejects_supported_lane_without_raw_artifact(tmp_path):
    checker = _load_checker()
    manifest = _valid_manifest()
    manifest["assays"][0]["raw_artifacts"].append("benchmarks/raw/missing.csv")
    _write_assay_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("missing raw artifact" in failure and "missing.csv" in failure for failure in result.failures)


def test_assay_evidence_rejects_gate_without_make_target(tmp_path):
    checker = _load_checker()
    manifest = _valid_manifest()
    manifest["assays"][0]["gates"].append("make missing-gate")
    _write_assay_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("missing make target" in failure and "missing-gate" in failure for failure in result.failures)


def test_assay_evidence_requires_next_public_evidence_for_planned_lanes(tmp_path):
    checker = _load_checker()
    manifest = _valid_manifest()
    del manifest["assays"][2]["next_public_evidence"]
    _write_assay_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("next_public_evidence" in failure and "perturb_seq" in failure for failure in result.failures)


def test_assay_evidence_requires_evidence_discipline_for_non_planned_lanes(tmp_path):
    checker = _load_checker()
    manifest = _valid_manifest()
    del manifest["assays"][0]["commands"]
    del manifest["assays"][1]["comparator_semantics"]
    del manifest["assays"][1]["validation"]
    _write_assay_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("commands" in failure and "crispr_guide_counting" in failure for failure in result.failures)
    assert any("comparator_semantics" in failure and "inline_barcode" in failure for failure in result.failures)
    assert any("validation" in failure and "inline_barcode" in failure for failure in result.failures)


def test_assay_evidence_rejects_raw_artifact_rows_without_command(tmp_path):
    checker = _load_checker()
    _write_assay_repo(tmp_path)
    (tmp_path / "benchmarks" / "raw" / "barcode_demux.csv").write_text(
        "tool,workflow,command,exit_code\n"
        "dotmatch_demux,real_public_inline_barcode,,0\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("barcode_demux.csv" in failure and "command" in failure for failure in result.failures)


def test_assay_evidence_rejects_nonzero_validation_mismatches(tmp_path):
    checker = _load_checker()
    _write_assay_repo(tmp_path)
    (tmp_path / "benchmarks" / "raw" / "public_crispr_edlib_validation.csv").write_text(
        "sample,mismatches,checked_reads\nplasmid,1,100\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("public_crispr_edlib_validation.csv" in failure and "mismatches" in failure for failure in result.failures)


def test_assay_evidence_rejects_empty_raw_artifact(tmp_path):
    checker = _load_checker()
    _write_assay_repo(tmp_path)
    (tmp_path / "benchmarks" / "raw" / "public_crispr_repeated.csv").write_text(
        "tool,workflow,command,exit_code\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("public_crispr_repeated.csv" in failure and "must contain at least one data row" in failure for failure in result.failures)
