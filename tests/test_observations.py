"""Differential tests against an independent labeled-base reconstruction oracle."""
import random

import pytest

from editwitness.models import Allele, Assay, Edit, Interval
from editwitness.observations import observe_allele
from editwitness.sequence import apply_edits, reverse_complement

REFERENCE = "ACGTGCATAGCTACGATG"
ASSAY = Assay(id="test", left_primer=Interval(start=3, end=5), right_primer=Interval(start=13, end=15))


def oracle(reference, edits, assay):
    # This oracle does not use production interval predicates, coordinate mapping or reconstruction.
    bases = [(i, base) for i, base in enumerate(reference)]
    for edit in reversed(edits):
        bases[edit.start:edit.end] = [(None, b) for b in edit.sequence]
    locations = {tag: i for i, (tag, _) in enumerate(bases) if tag is not None}
    for primer in (assay.left_primer, assay.right_primer):
        ids = list(range(primer.start, primer.end))
        if any(i not in locations for i in ids):
            return "original_binding_site_disrupted", None, ()
        indices = [locations[i] for i in ids]
        if indices != list(range(indices[0], indices[0] + len(ids))):
            return "original_binding_site_disrupted", None, ()
    l = locations[assay.left_primer.start]
    r = locations[assay.right_primer.end - 1] + 1
    length = r - l
    if length < assay.min_product_bp or (assay.max_product_bp is not None and length > assay.max_product_bp):
        return "outside_product_bounds", length, ()
    lo = locations[assay.left_primer.end - 1] + 1
    hi = locations[assay.right_primer.start]
    insert = "".join(base for _, base in bases[lo:hi])
    if assay.readout == "paired_end":
        k = assay.read_bases
        reverse = "".join({"A": "T", "C": "G", "G": "C", "T": "A"}[base] for base in insert[-k:][::-1])
        reads = (insert[:k], reverse)
    else:
        reads = (insert,)
    return "potentially_observable", length, reads


EVENTS = [(start, end, "") for start in range(18) for end in range(start + 1, 19)]
EVENTS += [(p, p, "TT") for p in range(19)]
EVENTS += [(p, p + 1, "A" if REFERENCE[p] != "A" else "C") for p in range(18)]
EVENTS += [(start, end, "GGG") for start in range(0, 18, 2) for end in range(start + 1, 19, 3)]


@pytest.mark.parametrize("start,end,sequence", EVENTS)
@pytest.mark.parametrize("mode", ["full_insert", "paired_end"])
def test_every_small_event_against_independent_oracle(start, end, sequence, mode):
    edit = Edit(start=start, end=end, sequence=sequence)
    allele = Allele(id="edited", edits=(edit,))
    assay = ASSAY if mode == "full_insert" else ASSAY.model_copy(update={"readout": mode, "read_bases": 3})
    result = observe_allele(REFERENCE, allele, assay)
    expected = oracle(REFERENCE, allele.edits, assay)
    assert (result.status, result.product_length, result.reads) == expected


def test_multiple_edits_randomized_against_oracle():
    rng = random.Random(7901)
    reference = "".join(rng.choice("ACGT") for _ in range(100))
    assay = Assay(id="a", left_primer=Interval(start=10, end=20), right_primer=Interval(start=70, end=80))
    for _ in range(400):
        starts = sorted(rng.sample(range(0, 95, 5), 5))
        edits = tuple(Edit(start=p, end=p + rng.randrange(1, 4), sequence="A" * rng.randrange(5)) for p in starts)
        allele = Allele(id="a", edits=edits)
        result = observe_allele(reference, allele, assay)
        assert (result.status, result.product_length, result.reads) == oracle(reference, edits, assay)


@pytest.mark.parametrize("size,min_bp,max_bp,status", [
    (2, 10, 10, "potentially_observable"),
    (2, 11, 20, "outside_product_bounds"),
    (2, 1, 9, "outside_product_bounds"),
])
def test_product_bounds_are_inclusive(size, min_bp, max_bp, status):
    assay = ASSAY.model_copy(update={"min_product_bp": min_bp, "max_product_bp": max_bp})
    allele = Allele(id="del", edits=(Edit(start=8, end=8+size),))
    assert observe_allele(REFERENCE, allele, assay).status == status


def test_primer_sequences_are_not_counted_as_observed_insert():
    result = observe_allele(REFERENCE, Allele(id="ref"), ASSAY)
    assert result.reads == (REFERENCE[5:13],)


def test_delete_entire_insert_is_empty_signal_not_no_product():
    result = observe_allele(REFERENCE, Allele(id="del", edits=(Edit(start=5, end=13),)), ASSAY)
    assert result.reads == ("",)
    assert result.signal_id is not None
    assert result.product_length == 4


def test_reconstruction_and_reverse_complement():
    edits = (Edit(start=1, end=2, sequence="GG"), Edit(start=5, end=7), Edit(start=8, end=8, sequence="T"))
    assert apply_edits("ACGTACGT", edits) == "AGGGTATT"
    assert reverse_complement("ACGTGA") == "TCACGT"
