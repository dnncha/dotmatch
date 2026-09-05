import copy
import socket

from editwitness import analyze
from editwitness.models import Manifest


def test_demo_demonstrates_dropout_and_unresolvable_window_loss(demo):
    result = analyze(demo)
    assert [w.hypothesis_id for w in result.witnesses] == ["hidden_primer_deletion", "hidden_window_deletion"]
    assert result.plan.selected_assays == ("outer",)
    assert result.plan.resolved_hypotheses == ("hidden_primer_deletion",)
    assert result.plan.unresolved_hypotheses == ("hidden_window_deletion",)
    assert result.plan.cost_units == 2


def test_paired_end_read_gap_does_not_observe_product_length(paired):
    result = analyze(paired)
    assert len(result.witnesses) == 4
    obs = {(o.allele_id, o.assay_id): o for o in result.allele_observations}
    a, b = obs["intended", "inner"], obs["interior_deletion", "inner"]
    assert a.product_length != b.product_length
    assert a.signal_id == b.signal_id  # Latent product length is NOT an observed feature.


def test_allele_multiplicity_does_not_create_dosage_evidence(demo):
    result = analyze(demo)
    obs = {(o.hypothesis_id, o.assay_id): o.signal_ids for o in result.hypothesis_observations}
    assert obs["intended_biallelic", "inner"] == obs["hidden_primer_deletion", "inner"]
    assert len(obs["intended_biallelic", "inner"]) == 1


def test_more_assays_cannot_increase_equivalence_class(demo):
    before = analyze(demo)
    data = demo.model_dump(mode="json")
    data["assays"].append(data["candidates"].pop(0))
    after = analyze(Manifest.model_validate(data))
    assert {w.hypothesis_id for w in after.witnesses} < {w.hypothesis_id for w in before.witnesses}


def test_input_unchanged_and_deterministic(demo):
    before = copy.deepcopy(demo.model_dump(mode="json"))
    a = analyze(demo)
    b = analyze(demo)
    assert a == b
    assert len(a.result_sha256) == 64
    assert before == demo.model_dump(mode="json")


def test_no_alternatives_never_means_validated(demo):
    data = demo.model_dump(mode="json")
    data["hypotheses"] = data["hypotheses"][:1]
    result = analyze(Manifest.model_validate(data))
    assert result.conclusion == "distinguishable_only_within_declared_model"
    assert result.validation_status == "software-tested; not empirically validated"
    assert "completeness" in result.plan.note


def test_no_expected_signal_is_explicit(demo):
    data = demo.model_dump(mode="json")
    data["hypotheses"][0]["alleles"] = ["window_deleted", "window_deleted"]
    result = analyze(Manifest.model_validate(data))
    assert "NO_BASELINE_POSITIVE_SIGNAL" in {n.code for n in result.notices}
    assert "NO_EXPECTED_SIGNAL" in {n.code for n in result.notices}


def test_no_network_required(demo, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("network access is forbidden")
    monkeypatch.setattr(socket, "create_connection", forbidden)
    analyze(demo)


def test_html_sequence_preview_is_bounded_and_explicit(demo):
    from editwitness.report import render_report
    result = analyze(demo)
    observations = tuple(
        o.model_copy(update={"reads": ("A" * 4000,)}) if o.reads else o
        for o in result.allele_observations
    )
    report = render_report(result.model_copy(update={"allele_observations": observations}))
    assert "preview: first 2400 of 4000 bases; full sequence in JSON" in report
    assert "A" * 2401 not in report


def test_html_witness_preview_cap_does_not_change_complete_evidence(demo):
    from editwitness.models import Witness
    from editwitness.report import render_report
    result = analyze(demo)
    witnesses = tuple(result.witnesses[0].model_copy(update={"hypothesis_id": f"alternative_{i}"}) for i in range(60))
    result = result.model_copy(update={"witnesses": witnesses})
    report = render_report(result)
    assert report.count('<article class="witness">') == 50
    assert "Showing the first 50 of 60 counterexamples" in report
    assert len(result.witnesses) == 60
    assert len(result.model_dump(mode="json")["witnesses"]) == 60
