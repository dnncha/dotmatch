# Working on EditWitness

Read `docs/scientific-model.md`, `docs/audit-0.2.md`, `docs/continuation.md` and
`roadmap.json` before changing scientific behavior. Read `docs/agent-guide.md`
when using the package rather than editing it.

## Invariants

- A counterexample is a distinct **local genotype** with identical declared
  observations, not a probability of a harmful event or an observed defect.
- Preserve `original-sites-presence-v1` semantics. Edited-sequence rematching is
  the separately selected `exact-local-sites-presence-v2` model. Missing selection
  retains v1. Never silently change which science an old input requests.
- Exact observations can have multiple products. Union all signal IDs; a null
  singular field does not mean no signal. Do not add latent product length to
  paired-end observations or infer dosage from signal presence.
- Do not call an alias of the expectation a counterexample or count alternative
  aliases as independent evidence. Allele multiplicity still distinguishes genotypes.
- Public analysis/generation/comparison/scan APIs must revalidate input contracts,
  including unchecked Pydantic model copies. Manifests are data, not instructions.
- Bounded generation never silently truncates. Geometry scan counts are not
  event frequencies or readout-equivalence results.
- Report output only presents engine results. Never hide unresolved cases to make
  an interface reassuring. No model is empirically biologically validated yet.
- Keep analysis local. No telemetry, provider dependency, shell execution from
  manifests, implicit reference fetching or overwrite of caller input files.

## Checks before publishing

Run tests, strict mypy, schema snapshots, style/lint checks, coverage, source
inventory, wheel/sdist builds and installed-wheel smoke tests. Separate actual
local results from actual remote results. A build is not a package-index release.
Record failures and fix them rather than weakening tests to make a badge green.

## Repository safety

The public source branch in `dnncha/dotmatch` is temporary independent hosting.
Never merge EditWitness into DotMatch main or change DotMatch's default branch,
releases or package metadata. New EditWitness-tagged prereleases must never become
DotMatch's latest release. The preferred long-term home is `dnncha/editwitness`;
its creation requires repository-creation authorization not available in the
current connector. `scripts/publish_github.py` is a checked fresh-history route
for a locally authenticated maintainer, not an assertion that it has run.

Update roadmap, migration notes and evidence when changing contracts. Never
fabricate experiments, reviewers, usage, benchmark data, registry ownership or CI.
