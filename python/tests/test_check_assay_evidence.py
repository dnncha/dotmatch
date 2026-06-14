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
                    "benchmarks/raw/crispr_comparison_count_agreement_summary.csv",
                    "benchmarks/raw/crispr_comparison_hamming_k23_comparators.csv",
                ],
                "reports": ["docs/benchmarks/public_crispr/README.md"],
                "gates": ["make public-crispr-evidence-gate"],
                "claim_boundary": "Public MAGeCK/Yusa rows only; not universal CRISPR superiority.",
                "biological_unit": "per-read fixed-window guide assignment and per-sample guide counts",
                "unsupported_claims": ["screen interpretation", "gene essentiality inference"],
                "minimum_public_evidence": ["public FASTQ rows", "zero-mismatch assignment validation"],
                "commands": [
                    "make bench-public-crispr-repeated",
                    "make public-crispr-evidence-gate",
                ],
                "comparator_semantics": "MAGeCK exact-count agreement plus DotMatch-vs-Bowtie 1 Hamming k=2/k=3 rows and bounded Edlib assignment validation.",
                "validation": "Gate requires bounded Edlib mismatches equal zero; Hamming k=2 must clear >=8x vs Bowtie 1 and Hamming k=3 must clear >=2x vs Bowtie 1. This is not universal CRISPR superiority evidence.",
            },
            {
                "id": "inline_barcode",
                "label": "Inline barcode demultiplexing",
                "status": "gated",
                "raw_artifacts": ["benchmarks/raw/barcode_demux.csv"],
                "reports": ["docs/benchmarks/barcode_demux/README.md"],
                "gates": ["make barcode-comparison-gate"],
                "claim_boundary": "Fixed 8 bp Hamming k=1 public lane plus exact-prefix public lane; Levenshtein one-edit barcode lane remains synthetic fixture evidence.",
                "biological_unit": "per-read fixed-position barcode assignment",
                "unsupported_claims": ["arbitrary adapter trimming", "BCL demultiplexing"],
                "minimum_public_evidence": ["public barcode sheet", "exact-prefix comparator agreement", "fixed 8 bp Hamming k=1 Hamming-radius splitter agreement"],
                "next_public_evidence": "Add public barcode-sheet rows with comparator agreement.",
                "commands": ["make bench-barcode-demux"],
                "comparator_semantics": "Public exact-prefix lane plus fixed 8 bp Hamming k=1 lane with Hamming-radius splitter comparator.",
                "validation": "Gate checks deterministic demux execution, Hamming-radius splitter agreement, and real-data speed floors: k=0 must beat Cutadapt by >=5x, exact hash splitting by >=3x, fixed 8 bp Hamming k=1 must beat Cutadapt by >=5x, and Hamming-radius splitter by >=12x.",
            },
            {
                "id": "raw_bcl_demux",
                "label": "Classic per-cycle BCL demultiplexing",
                "status": "gated",
                "raw_artifacts": ["benchmarks/raw/bcl_demux.csv"],
                "reports": ["docs/benchmarks/bcl_demux/README.md"],
                "gates": ["make bcl-tiny-public-gate"],
                "claim_boundary": "Tiny public BCL rows cover a parser milestone only.",
                "biological_unit": "per-cluster classic per-cycle BCL parser output counts",
                "unsupported_claims": ["production demultiplexing replacement", "CBCL support"],
                "minimum_public_evidence": ["real run-folder rows", "production comparator validation"],
                "next_public_evidence": "Add real run-folder rows with production comparator validation.",
                "commands": ["make bench-bcl-10x", "make bcl-tiny-public-gate"],
                "comparator_semantics": "Tiny BCL rows validate count totals only.",
                "validation": "Tiny public gate checks parser milestone rows.",
            },
            {
                "id": "paired_combinatorial",
                "label": "Paired or combinatorial target assignment",
                "status": "smoke",
                "raw_artifacts": [],
                "reports": [],
                "gates": ["make cli-test"],
                "claim_boundary": "CLI regression only; no public pair comparator.",
                "biological_unit": "per-read paired fixed-window target assignment",
                "unsupported_claims": ["guide-pair effect inference", "cell-level calls"],
                "minimum_public_evidence": ["public dual-target FASTQ", "pair-level oracle validation"],
                "next_public_evidence": "Add public dual-guide rows with pair-level validation.",
                "commands": ["make cli-test"],
                "comparator_semantics": "Native CLI pair-count semantics regression.",
                "validation": "CLI test checks paired assignment diagnostics.",
            },
            {
                "id": "perturb_seq",
                "label": "Perturb-seq guide or feature assignment",
                "status": "planned",
                "raw_artifacts": [],
                "reports": [],
                "gates": [],
                "claim_boundary": "No public perturb-seq evidence claim yet.",
                "biological_unit": "per-read fixed-window guide assignment",
                "unsupported_claims": ["cell-level quantification", "perturbation-effect inference"],
                "minimum_public_evidence": ["public guide-capture FASTQ", "cell-level comparator before quantification claims"],
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
                "biological_unit": "per-read fixed-window feature barcode assignment",
                "unsupported_claims": ["cell hashing calls", "UMI/cell quantification"],
                "minimum_public_evidence": ["public feature FASTQ", "Cell Ranger-compatible comparator before cell-level claims"],
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
                "biological_unit": "per-read fixed-window primer or panel target assignment",
                "unsupported_claims": ["variant calling", "clinical interpretation"],
                "minimum_public_evidence": ["public amplicon FASTQ", "full-assay comparator before diagnostic claims"],
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
                "biological_unit": "per-read fixed-window oligo or adapter assignment",
                "unsupported_claims": ["adapter trimming", "read merging"],
                "minimum_public_evidence": ["public oligo/adapter FASTQ", "trimming comparator before trimming claims"],
                "next_public_evidence": "Add a public oligo or adapter FASTQ fixture and validation oracle.",
            },
        ],
    }


def _write_assay_repo(root: Path, manifest=None) -> None:
    files = {
        "Makefile": (
            "public-crispr-evidence-gate:\n\ttrue\n"
            "barcode-comparison-gate:\n\ttrue\n"
            "bcl-tiny-public-gate:\n\ttrue\n"
            "bench-bcl-10x:\n\ttrue\n"
            "cli-test:\n\ttrue\n"
            "bench-public-crispr-repeated:\n\ttrue\n"
            "bench-barcode-demux:\n\ttrue\n"
        ),
        "docs/assay-evidence.json": json.dumps(manifest or _valid_manifest(), indent=2) + "\n",
        "benchmarks/raw/public_crispr_repeated.csv": (
            "tool,workflow,command,exit_code\n"
            "dotmatch_hamming_k1,public_crispr_yusa,dotmatch count --targets guides.tsv,0\n"
        ),
        "benchmarks/raw/public_crispr_edlib_validation.csv": "sample,mismatches,checked_reads\nplasmid,0,100\n",
        "benchmarks/raw/crispr_comparison_count_agreement_summary.csv": "comparison,status,total_delta,differing_guides\nexact,ok,0,0\n",
        "benchmarks/raw/crispr_comparison_hamming_k23_comparators.csv": (
            "dataset,k,comparison,status,n_targets,dotmatch_reads_per_sec,bowtie1_reads_per_sec,dotmatch_assigned_reads,bowtie1_assigned_reads\n"
            "d,2,dotmatch_hamming_k2_vs_bowtie1,ok,50000,100,10,9,9\n"
        ),
        "benchmarks/raw/barcode_demux.csv": (
            "tool,workflow,command,exit_code\n"
            "dotmatch_demux,real_public_inline_barcode,dotmatch demux --barcodes barcodes.tsv,0\n"
        ),
        "benchmarks/raw/bcl_demux.csv": (
            "tool,workflow,command,exit_code\n"
            "dotmatch_bcl,public_10x_tiny_bcl,dotmatch bcl-demux --run-folder run,0\n"
        ),
        "docs/benchmarks/public_crispr/README.md": "# Public CRISPR\n",
        "docs/benchmarks/barcode_demux/README.md": "# Barcode Demux\n",
        "docs/benchmarks/bcl_demux/README.md": "# BCL Demux\n",
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


def test_assay_evidence_accepts_checked_in_manifest():
    checker = _load_checker()

    result = checker.audit(ROOT)

    assert result.failures == []


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
    perturb = next(assay for assay in manifest["assays"] if assay["id"] == "perturb_seq")
    del perturb["next_public_evidence"]
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


def test_assay_evidence_rejects_missing_strong_contract_fragments(tmp_path):
    checker = _load_checker()
    manifest = _valid_manifest()
    crispr = manifest["assays"][0]
    crispr["raw_artifacts"] = [
        artifact for artifact in crispr["raw_artifacts"]
        if artifact != "benchmarks/raw/crispr_comparison_hamming_k23_comparators.csv"
    ]
    barcode = manifest["assays"][1]
    barcode["validation"] = "Gate checks deterministic demux execution."
    _write_assay_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("crispr_guide_counting" in failure and "hamming_k23_comparators" in failure for failure in result.failures)
    assert any("inline_barcode" in failure and "real-data speed floors" in failure for failure in result.failures)
