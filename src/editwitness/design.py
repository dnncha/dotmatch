"""Bounded generation of explicit diploid challenges from a declared deletion grid."""
from __future__ import annotations

import hashlib

from .io import InputError, digest
from .models import Allele, DeletionGeneration, Edit, Hypothesis, Manifest
from .sequence import apply_edits

MAX_GENERATION_BASES = 200_000_000


def expand_deletions(
    manifest: Manifest, *, fixed_allele: str | None = None, max_new_hypotheses: int = 100
) -> Manifest:
    """Add all unique grid-derived local sequence states, or fail without partial output.

    Each hypothesis combines one unchanged expected allele with a single deletion
    applied to the reference haplotype. This is not an outcome-frequency model or
    an exhaustive list of repair outcomes. Equivalent repeat-shifted deletions
    are deduplicated by resulting local sequence, retaining the first geometry.
    """
    manifest = Manifest.model_validate(manifest)
    if type(max_new_hypotheses) is not int or not 1 <= max_new_hypotheses <= 999:
        raise InputError("max_new_hypotheses must be an integer from 1 to 999")
    grid = manifest.deletion_scan
    if grid is None:
        raise InputError("expand-deletions requires an explicit deletion_scan grid")
    expected = next(h for h in manifest.hypotheses if h.id == manifest.expected_hypothesis)
    reference = manifest.reference.sequence
    sequences = {a.id: apply_edits(reference, a.edits) for a in manifest.alleles}
    if fixed_allele is None:
        if len({sequences[allele_id] for allele_id in expected.alleles}) != 1:
            raise InputError(
                "heterozygous expectation requires an explicit --fixed-allele (Python: fixed_allele). "
                "Allele order must not silently choose which haplotype to preserve."
            )
        fixed = min(expected.alleles)
    else:
        fixed = fixed_allele
    if not isinstance(fixed, str) or fixed not in expected.alleles:
        raise InputError("fixed allele must belong to the expected hypothesis")
    allele_ids = {a.id for a in manifest.alleles}
    hypothesis_ids = {h.id for h in manifest.hypotheses}
    identities = {key: hashlib.sha256(value.encode()).hexdigest() for key, value in sequences.items()}
    by_sequence: dict[str, str] = {}
    for a in manifest.alleles:
        by_sequence.setdefault(sequences[a.id], a.id)
    states = {tuple(sorted(identities[a] for a in h.alleles)) for h in manifest.hypotheses}
    new_alleles: list[Allele] = []
    new_hypotheses: list[Hypothesis] = []
    valid = duplicate = 0
    generated_bases = 0

    def available(prefix: str, used: set[str]) -> str:
        name, suffix = prefix, 1
        while name in used:
            name = f"{prefix}_{suffix}"
            suffix += 1
        used.add(name)
        return name

    for start in range(grid.start_min, grid.start_max + 1, grid.step):
        for end in range(grid.end_min, grid.end_max + 1, grid.step):
            size = end - start
            if size < grid.min_length or (grid.max_length is not None and size > grid.max_length):
                continue
            generated_bases += len(reference) - size
            if generated_bases > MAX_GENERATION_BASES:
                raise InputError(
                    "deletion generation exceeds the 200-million-base reconstruction budget; "
                    "narrow the grid or increase its step. No partial manifest was returned."
                )
            valid += 1
            sequence = reference[:start] + reference[end:]
            identity = hashlib.sha256(sequence.encode()).hexdigest()
            state = tuple(sorted((identities[fixed], identity)))
            if state in states:
                duplicate += 1
                continue
            if len(new_hypotheses) >= max_new_hypotheses or len(hypothesis_ids) >= 1000:
                raise InputError("deletion expansion exceeds the hypothesis cap; narrow the grid or increase its step. "
                                 "No silent sampling or partial manifest was returned.")
            allele_id = by_sequence.get(sequence)
            if allele_id is None:
                if len(allele_ids) >= 128:
                    raise InputError("deletion expansion exceeds 128 alleles; narrow the grid or increase its step. "
                                     "No partial manifest was returned.")
                allele_id = available(f"ewdel_{start}_{end}", allele_ids)
                new_alleles.append(Allele(id=allele_id, edits=(Edit(start=start, end=end),),
                    description=f"Generated reference-haplotype deletion [{start},{end}); not an observed event."))
                by_sequence[sequence] = allele_id
            hypothesis_id = available(f"ewchallenge_{start}_{end}", hypothesis_ids)
            new_hypotheses.append(Hypothesis(id=hypothesis_id, alleles=(fixed, allele_id),
                description=f"Fixed expected allele {fixed}; reference-haplotype deletion [{start},{end})."))
            states.add(state)
    if valid == 0:
        raise InputError("deletion grid contains no valid deletions after length filtering")
    generation = DeletionGeneration(
        input_manifest_sha256=digest(manifest.model_dump(mode="json")), grid=grid,
        fixed_allele=fixed, valid_deletions=valid, deduplicated_states=duplicate,
        added_alleles=len(new_alleles), added_hypotheses=len(new_hypotheses),
        max_new_hypotheses=max_new_hypotheses,
    )
    data = manifest.model_dump(mode="json")
    data.update(schema_version="1.1", generation=generation.model_dump(mode="json"))
    data["alleles"].extend(a.model_dump(mode="json") for a in new_alleles)
    data["hypotheses"].extend(h.model_dump(mode="json") for h in new_hypotheses)
    return Manifest.model_validate(data)
