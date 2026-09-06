import json
from pathlib import Path

from editwitness import analyze
from editwitness.io import verify_result
from editwitness.models import Manifest
from editwitness.report import render_report
from test_cli_io import cli
from test_planner import assay
from editwitness.planner import plan_panel


def test_old_result_integrity_survives_new_optional_schema_fields():
    source = Path(__file__).parent / "fixtures/legacy-analysis-0.1.0a1.json"
    result = verify_result(source)
    assert result.package_version == "0.1.0a1"
    assert result.schema_version == "1.0"
    assert cli("verify", source).returncode == 0


def test_old_result_replay_requests_exact_original_package_version(tmp_path):
    source = Path(__file__).parent / "fixtures/legacy-analysis-0.1.0a1.json"
    manifest = tmp_path / "unused.json"
    manifest.write_text("{}", encoding="utf-8")
    result = cli("verify", source, "--manifest", manifest)
    assert result.returncode == 5
    assert "0.1.0a1" in json.loads(result.stderr)["message"]


def test_cli_init_never_silently_assumes_full_insert_observation():
    result = cli("init", "--fasta", "unused.fasta", "--left-primer", "AAC",
                 "--right-primer", "GCT", "--edit-position", 10, "--alternate", "T")
    assert result.returncode == 2
    assert "required" in json.loads(result.stderr)["message"]


def test_hypothesis_definition_and_candidate_evidence_appear_in_report(demo):
    report = render_report(analyze(demo))
    assert "Allele definition: window_deleted" in report
    assert "delete reference[0:900)" in report
    assert "outer / intended" in report
    assert "separating candidate" in report
    assert "exact-local-sequence-presence-v2" in report


def test_equivalent_candidate_dominance_enables_exact_search():
    candidates = tuple(assay(f"a{i:02d}", i + 1) for i in range(24))
    result = plan_panel(candidates, {a.id: {"h"} for a in candidates}, {"h"})
    assert result.algorithm == "exhaustive_minimum_cost"
    assert result.selected_assays == ("a00",)
    assert len(result.dominated_candidates) == 23


def test_cli_model_comparison_exposes_legacy_representation_dependence(demo, tmp_path):
    data = demo.model_dump(mode="json")
    sequence = demo.reference.sequence
    alternative = next(b for b in "ACGT" if b != sequence[500])
    expanded_sequence = sequence[195:450] + data["alleles"][1]["edits"][0]["sequence"] + sequence[451:500] + alternative
    data["alleles"].append({"id": "rescued", "edits": [{"start": 195, "end": 501, "sequence": expanded_sequence}]})
    data["hypotheses"].append({"id": "rescue_sensitive", "alleles": ["intended", "rescued"]})
    source = tmp_path / "input.json"
    source.write_text(Manifest.model_validate(data).model_dump_json(), encoding="utf-8")
    result = cli("compare-models", source)
    assert result.returncode == 0, result.stderr
    assert "rescue_sensitive" in json.loads(result.stdout)["equivalent_only_in_legacy"]
    assert "hidden_window_deletion" in json.loads(result.stdout)["equivalent_in_both"]


def test_cli_generated_challenges_are_actual_analysis_inputs(demo, tmp_path):
    data = demo.model_dump(mode="json")
    data["deletion_scan"] = dict(start_min=195, start_max=201, end_min=215, end_max=221, step=3)
    source = tmp_path / "original.json"
    source.write_text(json.dumps(data), encoding="utf-8")
    expanded = cli("expand-deletions", source)
    assert expanded.returncode == 0, expanded.stderr
    result = cli("analyze", "-", "--compact", stdin=expanded.stdout)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["generation"]["added_hypotheses"] == 9


def test_complete_witness_includes_reconstructable_final_alleles():
    source = Path(__file__).resolve().parents[1] / "examples/demo.json"
    result = cli("witness", source, "--hypothesis", "hidden_window_deletion", "--include-sequences")
    assert result.returncode == 0, result.stderr
    allele = next(a for a in json.loads(result.stdout)["alleles"] if a["id"] == "window_deleted")
    assert allele["final_sequence"] == ""
    assert allele["edits"] == [{"start": 0, "end": 900, "sequence": ""}]
