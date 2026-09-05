# Working on EditWitness

Read `docs/scientific-model.md`, `docs/architecture.md`, and `BUILD_STATUS.md`
before changing the engine. Read `roadmap.json` for prioritized bounded tasks.

## Non-negotiable invariants

- An observational counterexample is not evidence that a defect occurred.
- Absence of a declared counterexample never means "safe" or "biallelic confirmed".
- Primer coordinates are local, zero-based and half-open. Reverse oligos are 5′→3′.
- All edits use original reference coordinates. No hidden coordinate conversion.
- Sequence-presence sets discard dosage and read fractions intentionally.
- Paired-end observations do not include diagnostic product length or unsequenced bases.
- Original-site eligibility is not PCR thermodynamics or sequence-aware rematching.
- Unsupported events must be rejected or explicitly labeled, never silently treated as normal.
- No network access in the scientific core, telemetry, arbitrary manifest plugins, or API keys.
- Presentation must not compute or change scientific conclusions.
- Model semantics changes require a new model version and migration documentation.

## Efficient continuation

Select the highest-priority unblocked task in `roadmap.json`, but do not manufacture
biological data or validation. Complete a narrow acceptance criterion, add a
regression test first when fixing behavior, and update task status with evidence.
Do not repeatedly rediscover the architecture or replace the project with a new scaffold.

Run:

```bash
python -m pip install -e '.[dev]'
python -m pytest -q
python scripts/generate_schemas.py --check
python scripts/check_style.py
python -m mypy src/editwitness
python -m build
```

Keep independent test oracles independent: never implement expected results by
calling the same production interval/sequence functions. Synthetic fixtures must
remain explicitly labeled and reproducible.

Never push a research alpha as a clinical or empirically validated release. Do
not publish patient-derived data in issues. Do not merge an isolated GitHub
staging branch into DotMatch: export this package into its own repository using
the documented publishing script. The projects remain distinct.
