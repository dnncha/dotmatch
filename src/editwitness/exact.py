"""Exact local heteroprimer products, rematched on the reconstructed allele.

This is an explicitly idealized sequence model, not a PCR simulator. Both
F->R and R->F orientations are searched. Products are oriented to the F primer.
No sampling, primer mismatch tolerance or genome-wide specificity is inferred.
"""
from __future__ import annotations

import hashlib
import json
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from typing import Literal

from .io import InputError
from .models import Allele, AlleleObservation, Assay, ExactProduct
from .sequence import apply_edits, reverse_complement

MAX_SITE_MATCHES = 512
MAX_PRODUCTS = 128
MAX_PRODUCT_BASES = 20_000_000


@dataclass
class ProductBudget:
    """Cumulative bound on reconstructed product bases, before allocating read strings."""
    used: int = 0
    limit: int = MAX_PRODUCT_BASES

    def consume(self, amount: int) -> None:
        if amount > self.limit - self.used:
            raise InputError("exact-local product work exceeds 20-million-base budget; split the manifest")
        self.used += amount


def _hits(sequence: str, oligo: str) -> tuple[int, ...]:
    positions: list[int] = []
    cursor = 0
    while (found := sequence.find(oligo, cursor)) >= 0:
        positions.append(found)
        if len(positions) > MAX_SITE_MATCHES:
            raise InputError("primer has over 512 exact local matches; use a more specific assay or smaller window")
        cursor = found + 1
    return tuple(positions)


def observe_exact(reference: str, allele: Allele, assay: Assay, *,
                  edited_sequence: str | None = None,
                  budget: ProductBudget | None = None) -> AlleleObservation:
    sequence = apply_edits(reference, allele.edits) if edited_sequence is None else edited_sequence
    budget = ProductBudget() if budget is None else budget
    forward = reference[assay.left_primer.start:assay.left_primer.end]
    reverse = reverse_complement(reference[assay.right_primer.start:assay.right_primer.end])
    orientations: tuple[tuple[str, str, Literal["forward", "reverse"]], ...] = (
        (forward, reverse_complement(reverse), "forward"),
        (reverse, reverse_complement(forward), "reverse"),
    )
    products: list[ExactProduct] = []
    inward = excluded = 0
    for left_oligo, right_site, orientation in orientations:
        right_hits = _hits(sequence, right_site)
        for left in _hits(sequence, left_oligo):
            # Adjacent primers allow an empty insert, distinct from no product.
            min_right = left + len(left_oligo)
            first = bisect_left(right_hits, min_right)
            inward_here = len(right_hits) - first
            inward += inward_here
            lower = max(min_right, left + assay.min_product_bp - len(right_site))
            upper = len(sequence) if assay.max_product_bp is None else (
                left + assay.max_product_bp - len(right_site)
            )
            lo = bisect_left(right_hits, lower)
            hi = max(lo, bisect_right(right_hits, upper))
            retained = hi - lo
            excluded += inward_here - retained
            if len(products) + retained > MAX_PRODUCTS:
                raise InputError("over 128 exact local products for an allele/assay; no results were truncated")
            for right in right_hits[lo:hi]:
                end = right + len(right_site)
                budget.consume(end - left)
                insert = sequence[left + len(left_oligo):right]
                if orientation == "reverse":
                    insert = reverse_complement(insert)
                reads: tuple[str, ...]
                if assay.readout == "full_insert":
                    reads = (insert,)
                else:
                    assert assay.read_bases is not None
                    reads = (insert[:assay.read_bases], reverse_complement(insert[-assay.read_bases:]))
                signal = "seq:" + hashlib.sha256(
                    json.dumps(reads, separators=(",", ":"), ensure_ascii=True).encode()
                ).hexdigest()
                products.append(ExactProduct(
                    start=left, end=end, orientation=orientation, product_length=end-left,
                    reads=reads, signal_id=signal,
                ))
    products.sort(key=lambda p: (p.start, p.end, p.orientation))
    if not products:
        return AlleleObservation(
            allele_id=allele.id, assay_id=assay.id,
            status="outside_product_bounds" if inward else "no_exact_local_product",
            reason="No retained exact inward-facing heteroprimer product in the local sequence. "
                   "Mismatched, same-primer and nonlocal amplification remain unmodeled.",
            inward_exact_pairs=inward, products_excluded_by_bounds=excluded,
        )
    single = products[0] if len(products) == 1 else None
    return AlleleObservation(
        allele_id=allele.id, assay_id=assay.id, status="potentially_observable",
        reason=f"{len(products)} exact local product(s) retained. Every modeled product is assumed detectable; "
               "this is not a prediction of PCR yield, sensitivity or allele dosage.",
        product_length=single.product_length if single else None,
        reads=single.reads if single else (), signal_id=single.signal_id if single else None,
        products=tuple(products), inward_exact_pairs=inward, products_excluded_by_bounds=excluded,
    )


def signal_ids(observation: AlleleObservation) -> tuple[str, ...]:
    """Project diagnostic products onto only the declared sequence-presence readout."""
    if observation.products:
        return tuple(sorted({p.signal_id for p in observation.products}))
    return (observation.signal_id,) if observation.signal_id is not None else ()
