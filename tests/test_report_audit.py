"""Adversarial report coverage: evidence stays available while HTML remains bounded."""
import pytest
from pydantic import ValidationError

from editwitness import analyze, expand_deletions
from editwitness.exact import observe_exact
from editwitness.models import Allele, Assay, Interval, Manifest
from editwitness.report import render_report


def multisignal_manifest():
    return Manifest.model_validate({
        "schema_version": "1.1", "observation_model": "exact-local-sequence-presence-v2",
        "reference": {"name": "synthetic multiple products", "sequence": "AACGGTAGCT", "synthetic": True},
        "alleles": [{"id": "r", "edits": []}, {"id": "multi", "edits": [
            {"start": 10, "end": 10, "sequence": "AGCTACGGTTAACTTTTGCT"}]}],
        "hypotheses": [{"id": "expected", "alleles": ["multi", "multi"]},
                       {"id": "alternative", "alleles": ["multi", "r"]}],
        "expected_hypothesis": "expected",
        "assays": [{"id": "short", "left_primer": {"start": 0, "end": 3},
                    "right_primer": {"start": 7, "end": 10}}],
    })


def test_multisignal_engine_and_bounded_report(monkeypatch):
    import editwitness.report as report
    result = analyze(multisignal_manifest())
    assert any(n.code == "MULTIPLE_LOCAL_PRODUCTS" for n in result.notices)
    assert result.witnesses[0].hypothesis_id == "alternative"
    multi = next(o for o in result.allele_observations if o.allele_id == "multi")
    assert multi.signal_id is None and len(multi.signal_ids) > 1
    observed = next(h for h in result.hypothesis_observations if h.hypothesis_id == "expected")
    assert observed.signal_ids == multi.signal_ids
    monkeypatch.setattr(report, "MAX_REPORT_PRODUCTS", 1)
    monkeypatch.setattr(report, "MAX_SEQUENCE_PREVIEW_BASES", 1)
    monkeypatch.setattr(report, "MAX_REPORT_SEQUENCE_BASES", 1)
    html = render_report(result)
    assert "every product is in the JSON" in html
    assert "report preview budget reached" in html
    assert "preview: first 1" in html
    assert "final-allele sites" in html
    assert len(multi.products) > 1  # Rendering did not mutate or discard source evidence.


def test_html_witness_limit_is_declared_not_silent(demo, monkeypatch):
    import editwitness.report as report
    result = analyze(demo)
    monkeypatch.setattr(report, "MAX_REPORT_WITNESSES", 1)
    html = render_report(result)
    assert "Showing the first 1 of 2 counterexamples" in html
    assert len(result.witnesses) == 2


def test_generated_provenance_visible_in_html(demo):
    data = demo.model_dump(mode="json")
    data["deletion_scan"] = dict(start_min=195, start_max=196, end_min=215, end_max=216)
    expanded = expand_deletions(Manifest.model_validate(data))
    html = render_report(analyze(expanded))
    assert "Declared generation: 4 grid deletions" in html
    assert "not outcome frequencies" in html


def test_no_alternatives_is_not_a_positive_validation(demo):
    data = demo.model_dump(mode="json")
    data["hypotheses"] = [h for h in data["hypotheses"] if h["id"] == data["expected_hypothesis"]]
    html = render_report(analyze(Manifest.model_validate(data)))
    assert "No alternative matched within the declared model" in html
    assert "not proof of completeness" in html


def test_long_map_labels_keep_complete_accessible_identity(demo):
    data = demo.model_dump(mode="json")
    identifier = "assay_" + "x" * 70
    data["assays"][0]["id"] = identifier
    html = render_report(analyze(Manifest.model_validate(data)))
    assert f"<title>{identifier}</title>" in html
    assert "min-width:940px" in html
    assert 'tabindex="0" role="region"' in html


def test_identical_oligos_rejected_by_manifest_and_lower_level():
    data = multisignal_manifest().model_dump(mode="json")
    data["reference"]["sequence"] = "AACGGTTGTT"
    with pytest.raises(ValidationError, match="identical"):
        Manifest.model_validate(data)
    assay = Assay(id="same", left_primer=Interval(start=0, end=3),
                  right_primer=Interval(start=7, end=10))
    from editwitness.io import InputError
    with pytest.raises(InputError, match="identical"):
        observe_exact("AACGGTTGTT", Allele(id="r"), assay)
