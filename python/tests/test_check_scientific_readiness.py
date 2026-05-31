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
            "python-package-test:\n\ttrue\n"
            "release-ready:\n\ttrue\n"
        ),
        "docs/scientific-readiness.json": json.dumps(manifest or _manifest(), indent=2) + "\n",
        "tests/test_cli_fastq.sh": "#!/bin/sh\n",
        "src/qda.c": "int main(void) { return 0; }\n",
        "benchmarks/raw/public_crispr_edlib_validation.csv": "sample,mismatches\ns,0\n",
        "docs/assay-evidence.json": "{}\n",
        "docs/packaging.md": "# Packaging\n",
        "docs/release-process.md": "# Release\n",
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
    manifest["controls"] = [control for control in manifest["controls"] if control["id"] != "memory_safety"]
    _write_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("missing required scientific readiness control: memory_safety" in failure for failure in result.failures)


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
