"""Bounded, explicit reference-haplotype deletion hypotheses; never event probabilities."""
from __future__ import annotations

from .io import InputError, digest
from .models import Allele, Edit, GenerationProvenance, Hypothesis, Manifest, validated_manifest
from .sequence import apply_edits

MAX_GENERATION_PAIRS = 5000


def expand_deletions(manifest: Manifest) -> Manifest:
    """Add expected-allele + reference-deletion alternatives from the declared grid.

    This deliberately models deletion of the reference haplotype, not deletion
    superimposed on the intended edit. A homozygous local expectation is required
    to avoid silently choosing a phase. Existing input alleles are not modified.
    No event is silently truncated to fit the allele/work limits.
    """
    manifest = validated_manifest(manifest)
    grid = manifest.deletion_scan
    if grid is None:
        raise InputError("expand-deletions requires a deletion_scan grid")
    if manifest.generation is not None:
        raise InputError("manifest already contains generated hypotheses; start from its original source")
    pairs = ((grid.start_max-grid.start_min)//grid.step+1) * ((grid.end_max-grid.end_min)//grid.step+1)
    if pairs > MAX_GENERATION_PAIRS:
        raise InputError("hypothesis generation is limited to 5,000 endpoint pairs; increase step or narrow ranges")
    expected = next(h for h in manifest.hypotheses if h.id == manifest.expected_hypothesis)
    sequences = {a.id: apply_edits(manifest.reference.sequence, a.edits) for a in manifest.alleles}
    if sequences[expected.alleles[0]] != sequences[expected.alleles[1]]:
        raise InputError("expansion requires a homozygous local expectation; heterozygous phasing must be explicit")
    partner = expected.alleles[0]
    known = {sequence: identifier for identifier, sequence in sequences.items()}
    genotypes = {tuple(sorted(sequences[a] for a in h.alleles)) for h in manifest.hypotheses}
    allele_ids = {a.id for a in manifest.alleles}
    hypothesis_ids = {h.id for h in manifest.hypotheses}
    alleles = list(manifest.alleles)
    hypotheses = list(manifest.hypotheses)
    considered = duplicates = 0
    reference = manifest.reference.sequence
    for start in range(grid.start_min, grid.start_max+1, grid.step):
        for end in range(grid.end_min, grid.end_max+1, grid.step):
            length = end-start
            if length < grid.min_length or (grid.max_length is not None and length > grid.max_length):
                continue
            considered += 1
            sequence = reference[:start] + reference[end:]
            genotype = tuple(sorted((sequences[partner], sequence)))
            if genotype in genotypes:
                duplicates += 1
                continue
            identifier = known.get(sequence)
            if identifier is None:
                identifier = f"del_{start}_{end}"
                if identifier in allele_ids:
                    raise InputError(f"generated allele id collides with an input id: {identifier}")
                if len(alleles) >= 128:
                    raise InputError("expanded manifest would exceed 128 alleles; narrow the grid; nothing was truncated")
                alleles.append(Allele(id=identifier, edits=(Edit(start=start, end=end),),
                                     description="Generated reference-haplotype deletion hypothesis, not an observed event."))
                allele_ids.add(identifier)
                known[sequence] = identifier
            hypothesis_id = f"expected_plus_del_{start}_{end}"
            if hypothesis_id in hypothesis_ids:
                raise InputError(f"generated hypothesis id collides with an input id: {hypothesis_id}")
            if len(hypotheses) >= 1000:
                raise InputError("expanded manifest would exceed 1,000 hypotheses; narrow the grid")
            hypotheses.append(Hypothesis(id=hypothesis_id, alleles=(partner, identifier)))
            hypothesis_ids.add(hypothesis_id)
            genotypes.add(genotype)
    if not considered:
        raise InputError("scan grid contains no valid deletions after length filtering")
    data = manifest.model_dump(mode="python")
    data.update(schema_version="1.1", alleles=tuple(alleles), hypotheses=tuple(hypotheses),
                generation=GenerationProvenance(
                    source_manifest_sha256=digest(manifest.model_dump(mode="json")), grid=grid,
                    enumerated_deletions=considered, added_alleles=len(alleles)-len(manifest.alleles),
                    duplicate_local_sequences=duplicates, paired_with_expected_allele=partner,
                    caveat="Finite reference-haplotype deletion grid paired with one expected allele. "
                           "No biological frequency, repair likelihood, phase inference or exhaustive outcome coverage. "
                           "Sequence-identical local diploid states are collapsed; no unique state is silently truncated.",
                ))
    return Manifest.model_validate(data)
