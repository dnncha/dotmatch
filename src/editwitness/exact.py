"""Sequence-invariant, bounded enumeration of exact local heteroprimer products.

This is a deterministic measurement model, not a model of PCR efficiency. It
includes new and rescued exact sites and both inward-facing orientations. Every
eligible heteroprimer product contributes signal; single-primer products,
mismatches, off-window sites and stochastic sampling remain outside the model.
"""
from __future__ import annotations

import hashlib
import json
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from typing import Literal

from .io import InputError
from .models import Allele, AlleleObservation, Assay, Interval, ProductObservation
from .sequence import apply_edits, reverse_complement

MAX_SITE_HITS = 10_000
MAX_PRODUCTS_PER_OBSERVATION = 4_096
MAX_TOTAL_PRODUCTS = 20_000
MAX_EVIDENCE_BASES = 20_000_000


@dataclass
class EvidenceBudget:
    """One cumulative budget shared by all observations in an analysis."""
    products: int = 0
    bases: int = 0

    def consume(self, products: int, bases: int) -> None:
        if self.products + products > MAX_TOTAL_PRODUCTS or self.bases + bases > MAX_EVIDENCE_BASES:
            raise InputError(
                "exact-local evidence budget exceeded; narrow the reference, use more specific "
                "primers or split the manifest. No truncated analysis was returned."
            )
        self.products += products
        self.bases += bases


def _hits(sequence: str, query: str) -> tuple[int, ...]:
    if not query:
        raise InputError("empty primer is unsupported")
    hits: list[int] = []
    start = 0
    while (start := sequence.find(query, start)) >= 0:
        if len(hits) == MAX_SITE_HITS:
            raise InputError("exact-local primer hit limit exceeded; use more specific primers or a smaller window")
        hits.append(start)
        start += 1
    return tuple(hits)


def signal_for(reads: tuple[str, ...]) -> str:
    # Same encoding as original-sites v1. Length metadata is deliberately absent.
    encoded = json.dumps(reads, separators=(",", ":"), ensure_ascii=True).encode()
    return "seq:" + hashlib.sha256(encoded).hexdigest()


def observe_exact(
    reference: str,
    allele: Allele,
    assay: Assay,
    *,
    edited_sequence: str | None = None,
    budget: EvidenceBudget | None = None,
) -> AlleleObservation:
    sequence = apply_edits(reference, allele.edits) if edited_sequence is None else edited_sequence
    budget = EvidenceBudget() if budget is None else budget
    f = assay.left_oligo or reference[assay.left_primer.start:assay.left_primer.end]
    r = assay.right_oligo or reverse_complement(
        reference[assay.right_primer.start:assay.right_primer.end]
    )
    if f == r:
        raise InputError("identical primer oligos have ambiguous read orientation; unsupported")
    arrangements: tuple[tuple[str, str, Literal["forward", "reverse"]], ...] = (
        (f, reverse_complement(r), "forward"),
        (r, reverse_complement(f), "reverse"),
    )
    products: list[ProductObservation] = []
    any_pair = False
    for left_query, right_query, orientation in arrangements:
        left_hits, right_hits = _hits(sequence, left_query), _hits(sequence, right_query)
        for left in left_hits:
            # Non-overlapping sites; an empty primer-trimmed insert is still a signal.
            adjacent = left + len(left_query)
            if bisect_left(right_hits, adjacent) < len(right_hits):
                any_pair = True
            minimum = max(adjacent, left + assay.min_product_bp - len(right_query))
            maximum = (
                len(sequence) if assay.max_product_bp is None
                else left + assay.max_product_bp - len(right_query)
            )
            low, high = bisect_left(right_hits, minimum), bisect_right(right_hits, maximum)
            count = max(0, high - low)
            if len(products) + count > MAX_PRODUCTS_PER_OBSERVATION:
                raise InputError(
                    f"{assay.id}/{allele.id}: exact-local product limit exceeded; "
                    "no arbitrary first match or partial result was used"
                )
            for index in range(low, high):
                right = right_hits[index]
                insert_length = right - adjacent
                observed_bases = insert_length if assay.readout == "full_insert" else (
                    2 * min(insert_length, assay.read_bases or 0)
                )
                budget.consume(1, observed_bases)
                # For paired ends, never materialize the entire unsequenced gap.
                reads: tuple[str, ...]
                if assay.readout == "full_insert":
                    insert = sequence[adjacent:right]
                    reads = (insert if orientation == "forward" else reverse_complement(insert),)
                else:
                    assert assay.read_bases is not None
                    k = min(assay.read_bases, insert_length)
                    start_read = sequence[adjacent:adjacent + k]
                    end_read = reverse_complement(sequence[right - k:right])
                    reads = (start_read, end_read) if orientation == "forward" else (end_read, start_read)
                products.append(ProductObservation(
                    plus_left_site=Interval(start=left, end=adjacent),
                    plus_right_site=Interval(start=right, end=right + len(right_query)),
                    orientation=orientation,
                    product_length=right + len(right_query) - left,
                    reads=reads, signal_id=signal_for(reads),
                ))
    if not products:
        return AlleleObservation(
            allele_id=allele.id, assay_id=assay.id,
            status="outside_product_bounds" if any_pair else "no_exact_local_product",
            reason=(
                "Exact inward-facing heteroprimer sites exist, but all products fall outside the declared size bounds."
                if any_pair else
                "No exact inward-facing heteroprimer product exists in the supplied final allele sequence. "
                "Mismatched, single-primer and off-window products are not modeled."
            ),
        )
    products.sort(key=lambda p: (p.plus_left_site.start, p.plus_right_site.start, p.orientation))
    signals = {p.signal_id: p.reads for p in products}
    ids = tuple(sorted(signals))
    return AlleleObservation(
        allele_id=allele.id, assay_id=assay.id, status="potentially_observable",
        reason=f"{len(products)} exact local heteroprimer product(s); {len(ids)} distinct sequence signal(s). "
               "All eligible products are assumed detectable; this is not measured PCR performance.",
        product_length=products[0].product_length if len(products) == 1 else None,
        reads=signals[ids[0]] if len(ids) == 1 else (),
        signal_id=ids[0] if len(ids) == 1 else None,
        signal_ids=ids, products=tuple(products),
    )
