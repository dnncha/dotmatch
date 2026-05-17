import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATOR = ROOT / "scripts" / "generate_evidence_gallery.py"
CHECKER = ROOT / "scripts" / "check_evidence_gallery.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _scenario(scenario_id: str, role: str) -> dict:
    report = f"docs/benchmarks/{scenario_id}/README.md"
    raw = f"benchmarks/raw/{scenario_id}.csv"
    return {
        "id": scenario_id,
        "title": f"{scenario_id} evidence",
        "category": "Public data" if role == "known_good" else "Report zoo",
        "assay_type": "inline_barcode",
        "status": "supported",
        "dataset": f"{scenario_id} dataset",
        "condition": role.replace("_", " "),
        "gallery_roles": [role],
        "primary_report": report,
        "raw_artifacts": [raw],
        "commands": [f"make {scenario_id}"],
        "comparator_semantics": "Fixture comparator semantics for gallery tests.",
        "validation": "Fixture validation row for gallery tests.",
        "proves": ["The report is linked from a generated scenario page."],
        "limits": ["The fixture does not represent a biological claim."],
        "report_examples": [{"label": "Primary report", "path": report, "kind": "markdown"}],
    }


def _write_gallery_repo(root: Path) -> None:
    scenarios = [
        _scenario("public_crispr_yusa", "known_good"),
        _scenario("barcode_autopsy_review", "low_confidence"),
        _scenario("barcode_wrong_offset_fixture", "wrong_offset"),
        _scenario("barcode_unsafe_correction", "unsafe_correction"),
        _scenario("feature_barcode_10x", "known_good"),
        _scenario("perturb_seq_10x_guide_capture", "known_good"),
        _scenario("amplicon_artic_primer_start", "known_good"),
        _scenario("oligo_adapter_truseq_prefix", "known_good"),
        _scenario("bcl_tiny_classic", "gated_parser_milestone"),
    ]
    manifest = {
        "schema_version": 1,
        "title": "DotMatch Evidence Gallery",
        "summary": "Generated gallery for public evidence and report examples.",
        "scenarios": scenarios,
    }
    manifest_path = root / "docs" / "evidence-gallery" / "manifest.json"
    manifest_path.parent.mkdir(parents=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (root / "Makefile").write_text(
        "\n".join(f"{scenario['commands'][0].split()[1]}:\n\ttrue" for scenario in scenarios) + "\n",
        encoding="utf-8",
    )
    for scenario in scenarios:
        for path in [scenario["primary_report"], *scenario["raw_artifacts"]]:
            full = root / path
            full.parent.mkdir(parents=True, exist_ok=True)
            if full.suffix == ".csv":
                full.write_text("tool,command,exit_code\nfixture,dotmatch count,0\n", encoding="utf-8")
            else:
                full.write_text(f"# {scenario['title']}\n", encoding="utf-8")


def test_evidence_gallery_generator_writes_index_scenarios_and_report_zoo(tmp_path):
    generator = _load(GENERATOR, "generate_evidence_gallery")
    _write_gallery_repo(tmp_path)

    generator.write_gallery(tmp_path)

    assert (tmp_path / "docs" / "evidence-gallery" / "README.md").is_file()
    assert (tmp_path / "docs" / "evidence-gallery" / "report-zoo" / "README.md").is_file()
    assert (tmp_path / "docs" / "evidence-gallery" / "scenarios" / "barcode_autopsy_review.md").is_file()
    index = (tmp_path / "docs" / "evidence-gallery" / "README.md").read_text(encoding="utf-8")
    assert "barcode_autopsy_review evidence" in index


def test_evidence_gallery_checker_accepts_generated_gallery(tmp_path):
    generator = _load(GENERATOR, "generate_evidence_gallery")
    checker = _load(CHECKER, "check_evidence_gallery")
    _write_gallery_repo(tmp_path)
    generator.write_gallery(tmp_path)

    result = checker.audit(tmp_path)

    assert result.failures == []
    assert any("evidence gallery" in item for item in result.passed)


def test_evidence_gallery_checker_rejects_missing_artifact(tmp_path):
    generator = _load(GENERATOR, "generate_evidence_gallery")
    checker = _load(CHECKER, "check_evidence_gallery")
    _write_gallery_repo(tmp_path)
    generator.write_gallery(tmp_path)
    (tmp_path / "docs" / "benchmarks" / "barcode_autopsy_review" / "README.md").unlink()

    result = checker.audit(tmp_path)

    assert any("missing primary_report" in failure for failure in result.failures)


def test_evidence_gallery_checker_rejects_stale_generated_page(tmp_path):
    generator = _load(GENERATOR, "generate_evidence_gallery")
    checker = _load(CHECKER, "check_evidence_gallery")
    _write_gallery_repo(tmp_path)
    generator.write_gallery(tmp_path)
    readme = tmp_path / "docs" / "evidence-gallery" / "README.md"
    readme.write_text(readme.read_text(encoding="utf-8") + "\nmanual edit\n", encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("generated evidence gallery file is stale" in failure for failure in result.failures)
