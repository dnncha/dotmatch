"""Streaming deletion geometry audit; fractions are not biological frequencies."""
from __future__ import annotations

from ._version import MODEL_VERSION, __version__
from .io import InputError, digest, seal
from .models import Manifest, ScanAssayCounts, ScanExample, ScanResult


def scan_deletions(manifest: Manifest) -> ScanResult:
    grid = manifest.deletion_scan
    if grid is None:
        raise InputError("manifest has no deletion_scan configuration")
    counts = {a.id: [0, 0, 0] for a in manifest.assays}
    blind_examples: list[ScanExample] = []
    total = blind = 0
    names = ("potentially_amplifiable", "binding_site_disrupted", "outside_product_bounds")
    for start in range(grid.start_min, grid.start_max + 1, grid.step):
        for end in range(grid.end_min, grid.end_max + 1, grid.step):
            length = end - start
            if length < grid.min_length or (grid.max_length is not None and length > grid.max_length):
                continue
            total += 1
            statuses: dict[str, str] = {}
            any_amplifiable = False
            for assay in manifest.assays:
                left, right = assay.left_primer, assay.right_primer
                if (start < left.end and end > left.start) or (start < right.end and end > right.start):
                    index = 1
                else:
                    product_length = right.end - left.start - max(
                        0, min(end, right.end) - max(start, left.start)
                    )
                    outside = product_length < assay.min_product_bp or (
                        assay.max_product_bp is not None and product_length > assay.max_product_bp
                    )
                    index = 2 if outside else 0
                counts[assay.id][index] += 1
                statuses[assay.id] = names[index]
                any_amplifiable |= index == 0
            if not any_amplifiable:
                blind += 1
                if len(blind_examples) < 20:
                    blind_examples.append(ScanExample(
                        start=start, end=end, length=length, assay_statuses=statuses
                    ))
    if total == 0:
        raise InputError("scan grid contains no valid deletions after length filtering")
    return seal(ScanResult(
        package_version=__version__, model_version=MODEL_VERSION,
        manifest_sha256=digest(manifest.model_dump(mode="json")), grid=grid,
        enumerated_deletions=total, all_existing_assays_structurally_blind=blind,
        assays=tuple(ScanAssayCounts(
            assay_id=assay.id, potentially_amplifiable=counts[assay.id][0],
            binding_site_disrupted=counts[assay.id][1], outside_product_bounds=counts[assay.id][2],
        ) for assay in manifest.assays),
        blind_examples=tuple(blind_examples),
        caveat="Counts describe only the declared deletion grid on the supplied reference; "
               "they are not event probabilities or empirical sensitivity. This geometry scan "
               "does not evaluate readout equivalence, compound edits, the intended allele, "
               "new primer sites, nonspecific amplification or experimental sampling.",
    ))
