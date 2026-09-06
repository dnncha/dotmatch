# Continuing EditWitness after the 0.2 audit

## Start from evidence, not the earlier chat

Read `BUILD_STATUS.md`, `docs/audit-0.2.md`, `docs/scientific-model.md` and
`roadmap.json`. Inspect the current GitHub branch and CI before claiming any
version is published. The project remains software-tested, not empirically
biologically validated.

The public temporary source is the independent `editwitness/public` branch in
`dnncha/dotmatch`. Never merge it into DotMatch or overwrite DotMatch main. The
preferred standalone destination is `dnncha/editwitness`; creation and PyPI
publication are not implied by source delivery. Use the checked publisher only
with explicit creation capability, verified source inventory and a correctly
authenticated maintainer. It creates fresh history and refuses existing repos.

## Implemented in 0.2

Exact reconstructed-sequence primer matching searches both heteroprimer
orientations and retains every bounded in-range product. v1 is still available
and remains the default for omitted-model manifests. New `demo`/`init` inputs
select v2 explicitly. Canonical local genotype pairs prevent aliases inflating
counterexamples. No alternatives and no expected baseline signal are explicit
non-reassuring states. Bounded generation and model comparison are real commands.

Tests independently reconstruct edited sequences and exhaustively match oligos;
they also check representation invariance, rescue, reverse products, multiple
signals, read gaps, duplicate states, invalid API copies and resource refusal.
Independent software oracles do not constitute biological ground truth.

## Next scientific work

Obtain an independent genome-engineering review of the all-products-detectable
assumption, same-primer/nonlocal exclusions, local genotype equivalence and
reference-deletion phase semantics. Acquire permitted, exact benchmark inputs
with reference sequence, primer geometry, edit structures, read configuration
and independent measurements. Record source/checksums/exclusions before scoring.
Report failures and false reassurance, not just convenient examples.

Do not add probabilistic risk, allele dosage, sequencing sensitivity, a clinical
label or automatic clone-release decisions without suitable models and evidence.
Do not silently broaden a generation grid, drop products or rename a model while
retaining its old identifier. Keep source inputs and upstream caller outputs intact.

## Engineering and release checks

Run all tests, strict mypy, schemas, source inventory, builds, installed-wheel
smoke tests and coverage. Preserve output compatibility or bump schema/model
versions with migration notes. Cross-platform CI must pass for the exact shipped
commit. A package-index release needs verified ownership and authorization.
Describe what actually happened, including any partial distribution status.
