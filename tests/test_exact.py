"""An independent, deliberately slow, substring-pair oracle for exact local products."""
import random

import pytest

from editwitness.exact import ProductBudget, observe_exact, signal_ids
from editwitness.io import InputError
from editwitness.models import Allele, Assay, Edit, Interval
from editwitness.observations import observe_allele


def rc(value):
    return ''.join({'A':'T', 'T':'A', 'C':'G', 'G':'C'}[base] for base in value[::-1])


def oracle(ref, allele, assay):
    # Apply replacements from the end using list slicing; no production sequence helpers.
    edited = list(ref)
    for event in reversed(allele.edits):
        edited[event.start:event.end] = list(event.sequence)
    edited = ''.join(edited)
    forward = ref[assay.left_primer.start:assay.left_primer.end]
    reverse = rc(ref[assay.right_primer.start:assay.right_primer.end])
    outputs = []
    for left_word, right_word, orientation in [(forward, rc(reverse), 'forward'),
                                               (reverse, rc(forward), 'reverse')]:
        for start in range(len(edited)+1):
            if edited[start:start+len(left_word)] != left_word:
                continue
            for right in range(start+len(left_word), len(edited)+1):
                if edited[right:right+len(right_word)] != right_word:
                    continue
                end = right+len(right_word)
                if end-start < assay.min_product_bp:
                    continue
                if assay.max_product_bp is not None and end-start > assay.max_product_bp:
                    continue
                insert = edited[start+len(left_word):right]
                if orientation == 'reverse':
                    insert = rc(insert)
                k = assay.read_bases
                reads = (insert,) if k is None else (insert[:k], rc(insert[-k:]))
                outputs.append((start, end, orientation, end-start, reads))
    return sorted(outputs)


@pytest.mark.parametrize('seed', range(60))
@pytest.mark.parametrize('mode', ['full_insert', 'paired_end'])
def test_exact_products_against_independent_oracle(seed, mode):
    rng = random.Random(seed)
    ref = ''.join(rng.choice('ACGT') for _ in range(32))
    edits = tuple(Edit(start=i, end=i+rng.randrange(1, 4),
                       sequence=''.join(rng.choice('ACGT') for _ in range(rng.randrange(5))))
                  for i in (3, 13, 23))
    allele = Allele(id='changed', edits=edits)
    assay = Assay(id='a', left_primer=Interval(start=3, end=6), right_primer=Interval(start=25, end=28),
                  readout=mode, read_bases=5 if mode == 'paired_end' else None,
                  min_product_bp=6+(seed % 3), max_product_bp=15+(seed % 20))
    actual = observe_exact(ref, allele, assay)
    assert [(p.start, p.end, p.orientation, p.product_length, p.reads) for p in actual.products] == oracle(ref, allele, assay)
    assert actual.inward_exact_pairs >= len(actual.products)
    assert actual.products_excluded_by_bounds == actual.inward_exact_pairs-len(actual.products)


def sample():
    ref = 'AACGTTATGCTACGTCCAGGACTTACCGT'
    assay = Assay(id='a', left_primer=Interval(start=0, end=6), right_primer=Interval(start=22, end=28))
    return ref, assay


def test_broad_and_minimal_replacements_have_identical_sequence_observations():
    ref, assay = sample()
    minimal = Allele(id='x', edits=(Edit(start=12, end=13, sequence='T'),))
    broad = Allele(id='y', edits=(Edit(start=0, end=20, sequence=ref[:12]+'T'+ref[13:20]),))
    assert observe_allele(ref, broad, assay).signal_id is None  # reproduced old model limitation
    a, b = observe_exact(ref, minimal, assay), observe_exact(ref, broad, assay)
    assert a.products == b.products
    assert signal_ids(a)


def test_deleted_primer_can_be_replaced_by_an_inserted_exact_site():
    ref, assay = sample()
    mutant = Allele(id='new_site', edits=(Edit(start=0, end=6), Edit(start=9, end=9, sequence=ref[:6])))
    assert observe_allele(ref, mutant, assay).signal_id is None
    result = observe_exact(ref, mutant, assay)
    assert result.products
    assert result.products[0].start == 3


def test_reverse_orientation_products_are_oriented_to_forward_oligo():
    ref, assay = sample()
    a = observe_exact(ref, Allele(id='ref'), assay)
    reversed_allele = Allele(id='r', edits=(Edit(start=0, end=len(ref), sequence=rc(ref)),))
    b = observe_exact(ref, reversed_allele, assay)
    assert {p.orientation for p in b.products} == {'reverse'}
    assert signal_ids(a) == signal_ids(b)


def test_multiple_products_are_retained_not_replaced_by_absent_signal():
    ref, assay = sample()
    allele = Allele(id='duplicate', edits=(Edit(start=10, end=10, sequence=ref[:6]),))
    result = observe_exact(ref, allele, assay)
    assert len(result.products) >= 2
    assert len(signal_ids(result)) >= 2
    assert result.signal_id is None  # singular diagnostic fields cannot represent multiple products
    assert result.status == 'potentially_observable'


def test_empty_insert_and_zero_products_are_distinct():
    ref, assay = sample()
    empty = observe_exact(ref, Allele(id='empty', edits=(Edit(start=6, end=22),)), assay)
    absent = observe_exact(ref, Allele(id='absent', edits=(Edit(start=0, end=len(ref)),)), assay)
    assert ('',) in {p.reads for p in empty.products}
    assert signal_ids(empty) != signal_ids(absent)


def test_resource_limits_refuse_rather_than_truncate():
    ref = 'A'*600 + 'CGT'
    assay = Assay(id='a', left_primer=Interval(start=0, end=1), right_primer=Interval(start=600, end=603))
    with pytest.raises(InputError, match='512'):
        observe_exact(ref, Allele(id='r'), assay)
    ref = 'A'*30 + 'C'*30
    assay = Assay(id='a', left_primer=Interval(start=0, end=1), right_primer=Interval(start=59, end=60))
    with pytest.raises(InputError, match='128'):
        observe_exact(ref, Allele(id='r'), assay)
    ref, assay = sample()
    with pytest.raises(InputError, match='budget'):
        observe_exact(ref, Allele(id='r'), assay, budget=ProductBudget(limit=1))


def test_paired_end_never_includes_unobserved_product_length():
    ref, a = sample()
    assay = a.model_copy(update={'readout':'paired_end', 'read_bases':3})
    intact = observe_exact(ref, Allele(id='r'), assay)
    deleted = observe_exact(ref, Allele(id='d', edits=(Edit(start=12, end=13),)), assay)
    assert intact.product_length != deleted.product_length
    assert signal_ids(intact) == signal_ids(deleted)
