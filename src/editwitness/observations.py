"""Explicit observation function for pristine original-site, sequence-presence assays."""
from __future__ import annotations

import hashlib
import json

from .models import Allele, AlleleObservation, Assay
from .sequence import apply_edits, disrupts_site, map_intact_base, reverse_complement


def observe_allele(reference: str, allele: Allele, assay: Assay) -> AlleleObservation:
    disrupted = [
        label for label, site in (("left", assay.left_primer), ("right", assay.right_primer))
        if any(disrupts_site(edit, site) for edit in allele.edits)
    ]
    if disrupted:
        return AlleleObservation(
            allele_id=allele.id, assay_id=assay.id,
            status="original_binding_site_disrupted",
            reason=f"Pristine {' and '.join(disrupted)} annotated primer site disrupted. "
                   "This response model emits no product; actual PCR behavior is not predicted.",
        )
    sequence = apply_edits(reference, allele.edits)
    left = map_intact_base(assay.left_primer.start, allele.edits)
    right = map_intact_base(assay.right_primer.start, allele.edits)
    left_length = assay.left_primer.end - assay.left_primer.start
    right_length = assay.right_primer.end - assay.right_primer.start
    product_length = right + right_length - left
    if product_length < assay.min_product_bp or (
        assay.max_product_bp is not None and product_length > assay.max_product_bp
    ):
        return AlleleObservation(
            allele_id=allele.id, assay_id=assay.id, status="outside_product_bounds",
            reason="Product excluded by the declared inclusive size bounds, not by a PCR-efficiency model.",
            product_length=product_length,
        )
    insert = sequence[left + left_length:right]
    reads: tuple[str, ...]
    if assay.readout == "full_insert":
        reads = (insert,)
    else:
        assert assay.read_bases is not None  # Enforced by the assay contract.
        reads = (insert[:assay.read_bases], reverse_complement(insert[-assay.read_bases:]))
    encoded = json.dumps(reads, separators=(",", ":"), ensure_ascii=True).encode()
    return AlleleObservation(
        allele_id=allele.id, assay_id=assay.id, status="potentially_observable",
        reason="Original sites retained and product within configured bounds. "
               "Signal is sequence presence, not abundance, dosage or a probability of detection.",
        product_length=product_length, reads=reads,
        signal_id="seq:" + hashlib.sha256(encoded).hexdigest(),
    )
