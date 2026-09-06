# Continue EditWitness, do not restart it

Read `BUILD_STATUS.md`, `roadmap.json`, `docs/scientific-model.md`, and the public
`editwitness-v0.2.0a2` release plus its linked CI run before claiming status.

## Publication boundary

The temporary public distribution host is `dnncha/dotmatch`, isolated branch
`editwitness/research-alpha-20260906`. This branch contains EditWitness source,
not a change to DotMatch main. **Never merge it into DotMatch.** The namespaced
prerelease is not marked as DotMatch's latest release. No PyPI publication is
claimed. Inspect the actual release before announcing completion or retrying it.

Repository creation was not available through the connected write actions. On
an already authenticated GitHub CLI, the clean source archive can create
`dnncha/editwitness` using `python scripts/publish_github.py --public --release`.
The helper refuses to overwrite an existing repository. Use its `--resume`
mode only for an exact reviewed-source match. Update canonical links after the
move, execute new CI, then configure trusted PyPI publishing after namespace
ownership is confirmed. Do not use account credentials from another project.

## Engineering acceptance

Run pytest, strict mypy, Ruff, schema and source inventory checks. Build and
exercise the installed wheel, `self-test`, and the extracted source distribution.
Record new test results; historical green CI is not evidence for changed code.

## Scientific next step

Obtain an independently reviewed, provenance-complete benchmark case with exact
reference, primer geometry, editing structures, read configuration and orthogonal
measurements. No existing synthetic test is experimental sensitivity evidence.
Do not introduce a fake perfect copy-number assay or infer unseen allele dosage.
Do not replace callers; implement one read-only format adapter only after its
semantics and observation uncertainty are specified with real fixtures.

## New regression invariants

Heterozygous expansion requires explicit `--fixed-allele`. Identical final DNA
must yield identical exact-local observations. Product size is not secretly
observed by paired-end reads. Budget failures do not return partial evidence.
A software self-test pass never becomes a biological safety claim.
