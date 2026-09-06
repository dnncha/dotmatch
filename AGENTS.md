# Working on EditWitness

Read `BUILD_STATUS.md`, `docs/scientific-model.md`, `docs/audit-0.2.0a1.md` and
`roadmap.json` before changing semantics. This is a finite-model research tool,
not clinical software or a caller. Never manufacture empirical validation,
independent review, adoption, remote CI or public distribution status.

## Invariants

- Coordinates are zero-based, half-open on the original local reference.
- Exact v2 uses final DNA, both orientations and all bounded heteroprimer
  products. Equivalent DNA must not depend on edit notation.
- Sequence presence is not dosage, counts or biological probability.
- Paired-end equivalence cannot use product length or unsequenced bases.
- Keep unresolvable alternatives; no silent grid/evidence truncation.
- New biological assumptions need an explicit model identity and migration.
- Preserve legacy input meaning; archived integrity and exact-version replay
  are different operations.
- Validate model instances at public API boundaries, even constructed/copied ones.
- Compact output is not a full replay artifact. Preserve structured errors.

## Development checks

Run tests, strict mypy, Ruff, schema drift, source hygiene and coverage. Add a
regression test before changing a scientific decision. Keep an independent
sequence oracle; don't share the optimized pairing helper with its test oracle.
Build the wheel and sdist and test imports outside the checkout. Regenerate
`release-files.json` only after reviewing source changes, then check it.

## Release boundaries

Never merge the old staging branch into DotMatch. This is a standalone project.
The publisher creates a new repo and optionally a prerelease after exact-SHA push
CI. It never force-pushes, uploads sequence samples, changes existing visibility,
or publishes to PyPI. Continue from `docs/continuation.md` and record actual
remote URLs/checks. The current session had no authenticated write capability;
this fact is not permission to claim a publication or bypass access controls.
