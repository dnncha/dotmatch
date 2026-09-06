# Validation and evidence

## What has been validated

The local audit tests deterministic software behavior, contract validation,
sequence reconstruction, exact local product enumeration, state equivalence,
finite panel selection, artifact integrity and packaging. See `BUILD_STATUS.md`
and `verification.json` for actual environment-specific results. Historical
0.1.0a1 remote CI results are archived separately and are **not** evidence that
0.2.0a1 has passed remote CI.

The exact model is compared with an independently written naive substring
oracle, including both orientations, multiple products, paired-end gaps,
sequence-rescuing replacements and equivalent edit representations. Generation
is checked for deterministic full-grid enumeration, deduplication, provenance,
limits and non-mutation. The optimizer is compared against independent subset
enumeration on small cases. Tests include legacy artifacts, invalid constructed
model instances, HTML escaping and release identity/checksum defenses.

Test counts include parametrized examples and seeded oracle cases, not that many
independent biological experiments. Coverage measures executed statements and
branches; it does not imply model correctness or scientific validity.

## What is not validated

No independent scientist has approved this release's observation function. No
adjudicated wet-lab dataset establishes predictive sensitivity or specificity.
No clinical use, experimental dropout frequency, PCR thermodynamics, genome-wide
specificity, allele dosage or safety certification is supported.

## Next scientific gates

First, obtain review of both the useful and misleading cases from independent
genome-engineering scientists. Keep disagreements in the repository, especially
about exact-match assumptions, off-window outcomes and negative results.

Second, curate a small benchmark with exact local references, original primers,
edit structures, read configuration, independent truth and permission to share.
Record exclusions and negative controls before scoring. Publications motivating
the problem are not automatically usable truth sets; missing assay metadata is
a reason to exclude or qualify a case, not to invent it.

Third, evaluate both false reassurance and excessive ambiguity. The important
outcome is a better validation decision, not a higher number of warnings.
Publish the scripts and per-case evidence, and separate synthetic software
checks from biological benchmark results.
