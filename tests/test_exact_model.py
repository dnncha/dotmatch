"""Adversarial and independently calculated checks for exact-local response semantics."""
import json
import random

import pytest
from pydantic import ValidationError

from editwitness import analyze
from editwitness._version import EXACT_MODEL_VERSION
from editwitness.exact import EvidenceBudget, observe_exact
from editwitness.io import InputError
from editwitness.models import Allele, Assay, Edit, Interval, Manifest
from editwitness.observations import observe_allele


def rc(text):
    return "".join({"A": "T", "C": "G", "G": "C", "T": "A"}[b] for b in reversed(text))


def independent_reads(sequence, f, r, read_bases=None, minimum=1, maximum=None):
    """Exhaust every literal substring pair on both strands, without production helpers."""
    outputs = []
    for strand in (sequence, rc(sequence)):
        for i in range(len(strand) - len(f) + 1):
            if strand[i:i + len(f)] != f:
                continue
            for j in range(i + len(f), len(strand) - len(r) + 1):
                if strand[j:j + len(r)] != rc(r):
                    continue
                length = j + len(r) - i
                if length < minimum or (maximum is not None and length > maximum):
                    continue
                insert = strand[i + len(f):j]
                if read_bases is None:
                    outputs.append((insert,))
                else:
                    outputs.append((insert[:read_bases], rc(insert[-read_bases:])))
    return outputs


def exact_manifest(demo):
    data = demo.model_dump(mode="json")
    data.update(schema_version="1.1", observation_model=EXACT_MODEL_VERSION)
    return Manifest.model_validate(data)


def test_sequence_invariance_rescues_unchanged_site_in_broad_replacement(demo):
    assay = demo.assays[0]
    compact = demo.alleles[1]
    intended_base = compact.edits[0].sequence
    broad = Allele(id="broad", edits=(Edit(start=195, end=451,
        sequence=demo.reference.sequence[195:450] + intended_base),))
    assert observe_allele(demo.reference.sequence, broad, assay).signal_id is None
    a = observe_exact(demo.reference.sequence, compact, assay)
    b = observe_exact(demo.reference.sequence, broad, assay)
    assert a.signal_ids == b.signal_ids
    assert a.products == b.products
    assert b.status == "potentially_observable"


def test_renamed_or_reencoded_expected_state_is_not_a_counterexample(demo):
    data = demo.model_dump(mode="json")
    data["hypotheses"].append({"id": "renamed", "alleles": ["intended", "intended"]})
    for model in ("original-sites-presence-v1", EXACT_MODEL_VERSION):
        data.update(schema_version="1.1", observation_model=model)
        result = analyze(Manifest.model_validate(data))
        assert "renamed" not in {w.hypothesis_id for w in result.witnesses}
        assessment = next(h for h in result.hypotheses if h.hypothesis_id == "renamed")
        assert assessment.same_local_genomic_state_as_expected
        assert any(n.code == "IDENTICAL_LOCAL_STATES_EXCLUDED" for n in result.notices)


def test_same_state_with_different_allele_ids_is_not_an_alternative(demo):
    data = exact_manifest(demo).model_dump(mode="json")
    duplicate = dict(data["alleles"][1], id="alias")
    data["alleles"].append(duplicate)
    data["hypotheses"].append({"id": "aliased", "alleles": ["alias", "intended"]})
    result = analyze(Manifest.model_validate(data))
    assert "aliased" not in {w.hypothesis_id for w in result.witnesses}


@pytest.mark.parametrize("read_bases", [None, 1, 3, 20])
def test_randomized_independent_both_strand_oracle(read_bases):
    rng = random.Random(9274)
    checked = 0
    for _ in range(150):
        reference = "".join(rng.choice("ACGT") for _ in range(28))
        f, r = reference[2:5], rc(reference[22:25])
        if f == r:
            continue
        start = rng.randint(0, 28)
        end = rng.randint(start, 28)
        insertion = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 8)))
        if reference[start:end] == insertion:
            continue
        edited = reference[:start] + insertion + reference[end:]
        minimum, maximum = rng.choice([1, 6, 10]), rng.choice([None, 12, 40])
        assay = Assay(id="test", left_primer=Interval(start=2, end=5),
            right_primer=Interval(start=22, end=25),
            readout="full_insert" if read_bases is None else "paired_end",
            read_bases=read_bases, min_product_bp=minimum, max_product_bp=maximum)
        allele = Allele(id="edited", edits=(Edit(start=start, end=end, sequence=insertion),))
        result = observe_exact(reference, allele, assay)
        expected = independent_reads(edited, f, r, read_bases, minimum, maximum)
        assert {p.reads for p in result.products} == set(expected)
        assert len(result.products) == len(expected)
        assert len(result.signal_ids) == len(set(expected))
        checked += 1
    assert checked > 100


def test_new_products_both_orientations_and_multi_signal_union():
    reference = "AACGGTAGCT"
    assay = Assay(id="a", left_primer=Interval(start=0, end=3), right_primer=Interval(start=7, end=10))
    # Insert a reversed entire product and another forward product.
    addition = rc(reference) + "AACTTTTGCT"
    allele = Allele(id="expanded", edits=(Edit(start=10, end=10, sequence=addition),))
    result = observe_exact(reference, allele, assay)
    expected = independent_reads(reference + addition, "AAC", "AGC")
    assert len(result.products) == len(expected) > 1
    assert {p.orientation for p in result.products} == {"forward", "reverse"}
    assert {p.reads for p in result.products} == set(expected)
    assert result.signal_id is None  # Must not mean "no signal" to the engine.
    assert result.signal_ids


def test_repeat_product_budget_fails_without_silently_choosing_first_product():
    reference = "AC" * 150
    assay = Assay(id="a", left_primer=Interval(start=0, end=2), right_primer=Interval(start=298, end=300))
    with pytest.raises(InputError, match="product limit"):
        observe_exact(reference, Allele(id="reference"), assay)
    bounded = assay.model_copy(update={"max_product_bp": 1})
    assert observe_exact(reference, Allele(id="reference"), bounded).status == "outside_product_bounds"


def test_global_evidence_budget_fails_closed():
    budget = EvidenceBudget(products=20_000)
    with pytest.raises(InputError, match="No truncated analysis"):
        budget.consume(1, 0)
    assert budget.products == 20_000
    with pytest.raises(InputError, match="budget"):
        EvidenceBudget(bases=20_000_000).consume(0, 1)


def test_hit_budget_is_bounded(monkeypatch):
    import editwitness.exact as exact
    monkeypatch.setattr(exact, "MAX_SITE_HITS", 2)
    with pytest.raises(InputError, match="hit limit"):
        exact._hits("AAAA", "A")
    with pytest.raises(InputError, match="empty primer"):
        exact._hits("AAAA", "")


def test_exact_empty_insert_is_observable_not_dropout():
    reference = "AACGGTAGCT"
    assay = Assay(id="a", left_primer=Interval(start=0, end=3), right_primer=Interval(start=7, end=10))
    allele = Allele(id="empty", edits=(Edit(start=3, end=7),))
    for bases in (None, 5):
        configured = assay.model_copy(update={"readout": "full_insert" if bases is None else "paired_end", "read_bases": bases})
        result = observe_exact(reference, allele, configured)
        assert result.signal_ids
        assert result.reads == (("",) if bases is None else ("", ""))


def test_public_engine_revalidates_model_copy_and_construct(demo):
    corrupt = demo.model_copy(update={"expected_hypothesis": "missing"})
    with pytest.raises(ValidationError):
        analyze(corrupt)
    bad_assay = demo.assays[0].model_copy(update={"cost_units": -1})
    corrupt = demo.model_copy(update={"assays": (bad_assay,)})
    with pytest.raises(ValidationError):
        analyze(corrupt)
    corrupt = Manifest.model_construct(**{**demo.model_dump(), "coordinate_system": "1-based"})
    with pytest.raises(ValidationError):
        analyze(corrupt)


def test_legacy_manifest_meaning_does_not_change_implicitly(demo):
    data = demo.model_dump(mode="json")
    data.pop("observation_model", None)
    data["schema_version"] = "1.0"
    assert analyze(Manifest.model_validate(data)).model_version == "original-sites-presence-v1"
    data["observation_model"] = EXACT_MODEL_VERSION
    with pytest.raises(ValidationError, match="schema_version 1.1"):
        Manifest.model_validate(data)


def test_exact_evidence_keeps_unobserved_allele_definition(demo):
    result = analyze(exact_manifest(demo))
    loss = next(a for a in result.allele_evidence if a.allele_id == "window_deleted")
    assert loss.sequence_length == 0
    assert loss.edits[0].start == 0 and loss.edits[0].end == len(demo.reference.sequence)
    assert len(loss.sequence_sha256) == 64
    assert result.plan.unresolved_hypotheses == ("hidden_window_deletion",)
    assert json.loads(result.model_dump_json())["model_version"] == EXACT_MODEL_VERSION
