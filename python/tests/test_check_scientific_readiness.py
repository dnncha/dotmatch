import importlib.util
import json
from pathlib import Path
from typing import Optional


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_scientific_readiness.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_scientific_readiness", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _manifest() -> dict:
    return {
        "schema_version": 1,
        "status": "evidence_bounded",
        "scope": "Known-target assignment only.",
        "not_validated_for": ["clinical diagnostics"],
        "controls": [
            {
                "id": "input_integrity",
                "status": "required",
                "evidence": ["tests/test_cli_fastq.sh"],
                "gates": ["make cli-test"],
                "acceptance": "Malformed inputs fail.",
            },
            {
                "id": "memory_safety",
                "status": "required",
                "evidence": ["src/qda.c"],
                "gates": ["make asan"],
                "acceptance": "Sanitizers pass.",
            },
            {
                "id": "oracle_validation",
                "status": "required",
                "evidence": ["benchmarks/raw/public_crispr_edlib_validation.csv"],
                "gates": ["make public-crispr-evidence-gate"],
                "acceptance": "Oracle mismatches are zero.",
            },
            {
                "id": "public_assay_evidence",
                "status": "required",
                "evidence": ["docs/assay-evidence.json"],
                "gates": ["make assay-evidence-ready"],
                "acceptance": "Claims are evidence-bounded.",
            },
            {
                "id": "public_claim_wording",
                "status": "required",
                "evidence": ["README.md"],
                "gates": ["make scientific-readiness-ready"],
                "acceptance": "Public wording avoids unsupported broad claims.",
            },
            {
                "id": "distribution_reproducibility",
                "status": "required",
                "evidence": ["docs/packaging.md"],
                "gates": ["make python-package-test"],
                "acceptance": "Packages are smoke-tested.",
            },
            {
                "id": "release_governance",
                "status": "required",
                "evidence": ["docs/release-process.md"],
                "gates": ["make release-ready"],
                "acceptance": "Release gates pass.",
            },
        ],
    }


def _write_repo(root: Path, manifest: Optional[dict] = None) -> None:
    files = {
        "Makefile": (
            "cli-test:\n\ttrue\n"
            "asan:\n\ttrue\n"
            "public-crispr-evidence-gate:\n\ttrue\n"
            "assay-evidence-ready:\n\ttrue\n"
            "scientific-readiness-ready:\n\ttrue\n"
            "python-package-test:\n\ttrue\n"
            "release-ready:\n\ttrue\n"
        ),
        "README.md": (
            "# DotMatch\n\n"
            "Deterministic known-target assignment.\n\n"
            "Evidence boundary: see DotMatch Evidence Notes. The strongest current evidence is native fixed-window indexed assignment, "
            "public CRISPR guide-counting comparisons, and checked public inline-barcode lanes; broader alignment, demultiplexing, "
            "screen-analysis, or BCL replacement claims need their own gates.\n"
        ),
        "docs/scientific-readiness.json": json.dumps(manifest or _manifest(), indent=2) + "\n",
        "tests/test_cli_fastq.sh": "#!/bin/sh\n",
        "src/qda.c": "int main(void) { return 0; }\n",
        "benchmarks/raw/public_crispr_edlib_validation.csv": "sample,mismatches\ns,0\n",
        "docs/assay-evidence.json": "{}\n",
        "docs/packaging.md": "# Packaging\n",
        "docs/release-process.md": "# Release\n",
        "docs/scientific-claims.md": (
            "# DotMatch Evidence Notes\n\n"
            "## Strongest Scoped Performance Evidence\n\n"
            "These are scoped to benchmark rows, not to global short-read alignment or all demultiplexing tasks.\n"
            "`make native-exact-gate` records large-library fixed-length `k=2` substitution rows and "
            "Levenshtein `k=2` insertion/deletion rows with a minimum `10.83x` speedup over Edlib, "
            "a minimum `9.58x` speedup over Edlib, and at most `1.00` verified candidates/read.\n"
            "`make crispr-comparison-gate` requires two real public datasets and "
            "DotMatch-vs-Bowtie 1 Hamming `k=2`/`k=3` rows. "
            "The current Bowtie 1 artifact records `9.71x` for Hamming `k=2` and "
            "the current Bowtie 1 artifact records `2.36x` for Hamming `k=3`.\n"
            "The fair guide-counter compatibility lane is Hamming `k=1`, no indels.\n"
            "`make public-crispr-evidence-gate` supports public rows, not universal CRISPR superiority.\n"
            "`make barcode-comparison-gate` requires public barcode evidence; "
            "the Levenshtein one-edit barcode lane remains synthetic fixture evidence.\n"
        ),
        "docs/benchmarks/native/README.md": (
            "# Native Edlib Benchmark Report\n\n"
            "## Gated Native Scaling Claims\n\n"
            "| claim | large_library_rows | min_speedup_vs_edlib | median_speedup_vs_edlib | max_verified_per_read | min_speedup_required | max_verified_required |\n"
            "| --- | --- | --- | --- | --- | --- | --- |\n"
            "| k=2 substitution indexed rows | 36 | 10.83 | 20.54 | 1.00 | 10.00 | 1.05 |\n"
            "| Levenshtein k=2 insertion/deletion rows | 18 | 9.58 | 15.58 | 1.00 | 8.00 | 25.00 |\n"
        ),
        "docs/benchmarks/crispr_comparison/README.md": (
            "# CRISPR Comparison Evidence\n\n"
            "## Hamming k2/k3 External Comparator Rows\n\n"
            "|dataset|k|records_per_sample|dotmatch_tool|bowtie1_tool|dotmatch_reads_per_sec|bowtie1_reads_per_sec|speedup|status|semantics|artifact|\n"
            "|---|---|---|---|---|---|---|---|---|---|---|\n"
            "|sanson_brunello|2|1000000|dotmatch_hamming_k2|bowtie1_crispr_hamming_k2|132272.5|13622.2|9.71|ok|Hamming k=2, no indels|k2.csv|\n"
            "|sanson_brunello|3|1000000|dotmatch_hamming_k3|bowtie1_crispr_hamming_k3|11900.7|5046.7|2.36|ok|Hamming k=3, no indels|k3.csv|\n"
        ),
        "docs/index.md": (
            "# DotMatch Documentation\n\n"
            "Evidence boundary: see DotMatch Evidence Notes. The strongest current evidence is native fixed-window indexed assignment, "
            "public CRISPR guide-counting comparisons, and checked public inline-barcode lanes; broader alignment, demultiplexing, "
            "screen-analysis, or BCL replacement claims need their own gates.\n"
        ),
    }
    for path, text in files.items():
        full = root / path
        full.parent.mkdir(parents=True, exist_ok=True)
        full.write_text(text, encoding="utf-8")


def test_scientific_readiness_accepts_minimal_manifest(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)

    result = checker.audit(tmp_path)

    assert result.failures == []


def test_scientific_readiness_accepts_checked_in_manifest():
    checker = _load_checker()

    result = checker.audit(ROOT)

    assert result.failures == []


def test_scientific_readiness_rejects_missing_required_control(tmp_path):
    checker = _load_checker()
    manifest = _manifest()
    manifest["controls"] = [control for control in manifest["controls"] if control["id"] != "public_claim_wording"]
    _write_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("missing required scientific readiness control: public_claim_wording" in failure for failure in result.failures)


def test_scientific_readiness_rejects_missing_gate_target(tmp_path):
    checker = _load_checker()
    manifest = _manifest()
    manifest["controls"][0]["gates"].append("make missing-gate")
    _write_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("missing make target" in failure and "missing-gate" in failure for failure in result.failures)


def test_scientific_readiness_rejects_missing_evidence(tmp_path):
    checker = _load_checker()
    manifest = _manifest()
    manifest["controls"][0]["evidence"].append("missing.tsv")
    _write_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("missing input_integrity evidence: missing.tsv" in failure for failure in result.failures)


def test_scientific_readiness_rejects_unsupported_public_claim_wording(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    (tmp_path / "README.md").write_text(
        "# DotMatch\n\nDotMatch is SOTA for fast exact short-DNA known-target assignment.\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("fast-exact positioning" in failure for failure in result.failures)
    assert any("SOTA language" in failure for failure in result.failures)


def test_scientific_readiness_rejects_missing_scoped_claim_contract(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    (tmp_path / "docs" / "scientific-claims.md").write_text(
        "# DotMatch Evidence Notes\n\n## Strongest Scoped Performance Evidence\n\n`make native-exact-gate`\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("global short-read" in failure for failure in result.failures)
    assert any("DotMatch-vs-Bowtie 1 Hamming" in failure for failure in result.failures)


def test_scientific_readiness_rejects_stale_native_claim_values(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    claims = tmp_path / "docs" / "scientific-claims.md"
    claims.write_text(
        claims.read_text(encoding="utf-8").replace("minimum `10.83x`", "minimum `10.99x`"),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("current native benchmark minimum" in failure for failure in result.failures)


def test_scientific_readiness_rejects_stale_crispr_claim_values(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    claims = tmp_path / "docs" / "scientific-claims.md"
    claims.write_text(
        claims.read_text(encoding="utf-8").replace("records `9.71x`", "records `9.99x`"),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("current CRISPR Bowtie 1 observed speedup" in failure for failure in result.failures)


def test_scientific_readiness_rejects_missing_readme_evidence_boundary(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    (tmp_path / "README.md").write_text("# DotMatch\n\nDeterministic known-target assignment.\n", encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("README.md" in failure and "Evidence boundary:" in failure for failure in result.failures)
    assert any("README.md" in failure and "public CRISPR" in failure for failure in result.failures)
    assert any("README.md" in failure and "guide-counting comparisons" in failure for failure in result.failures)


def test_scientific_readiness_rejects_missing_docs_index_evidence_boundary(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    (tmp_path / "docs" / "index.md").write_text("# DotMatch Documentation\n", encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("docs/index.md" in failure and "Evidence boundary:" in failure for failure in result.failures)
    assert any("docs/index.md" in failure and "public CRISPR" in failure for failure in result.failures)
    assert any("docs/index.md" in failure and "BCL replacement claims" in failure for failure in result.failures)
