# Build and verification status

Version: **0.2.0a1 research alpha**. Audit dated 2026-09-06.

This is a self-audit and software verification record. Independent scientific
review, an adjudicated biological benchmark and empirical assay validation have
not been performed. Historical 0.1 evidence is retained in `docs/verification.json`
and must not be used as proof for the changed 0.2 model.

## Local checks for the audited implementation

The final covered suite passed **726 tests** on Linux/Python 3.13.5, Pydantic
2.13.4 and pytest 9.0.2. Statement and branch totals, environment and command
checks are recorded in `docs/verification-0.2.json`. Coverage is software execution
coverage, not biological coverage or a probability that the implementation is correct.

The schema snapshots and dependency-free syntax/source hygiene check passed.
The latter is **not Ruff**. Strict mypy is configured as a remote gate; it was
not available in the offline local environment and is not claimed as locally run.

Both wheel and source distribution were built offline with setuptools. The
installed-wheel check uses a fresh temporary environment outside the source tree
with preinstalled dependencies, verifies the package import location, runs the
bundled demo, analysis, full-sequence witness, HTML, model comparison, schema and
result replay. Its observed result is recorded separately from build success.

The actual generated report was rendered in Chromium at 1440 and 390 pixels.
Both had zero document-level horizontal overflow, no scripts and no page errors.
Desktop and mobile screenshots were visually inspected.

## Remote and public-distribution evidence

The intended public source is `dnncha/dotmatch:editwitness/public`. This is
independent temporary hosting, **not a change to DotMatch main**. The publication
workflow runs a nine-job OS/Python matrix (Linux, macOS, Windows; Python 3.11,
3.13, 3.14), strict typing, coverage with a 95% combined floor, source inventory,
schema and build checks, and installed-wheel smoke tests.

Configured jobs are not evidence of passing jobs. The public build writes
`dist-public/publication.json` only after every required verification job passes;
it records the tested source SHA and actual run URL. Its neighboring checksums
identify the wheel, source archive and example evidence. A receipt is not an
authenticated provenance attestation or biological validation. Check the run and
exact commit before using a mutable branch in a reproducible workflow.

A separate EditWitness-tagged GitHub prerelease may be attempted after the checks,
but is not assumed to exist. The receipt records the actual outcome. GitHub
repository permissions may prevent release/tag creation for a branch that changes
workflows. Public source and versioned download assets do not imply a standalone
repository, a PyPI release, a DOI or empirical validation.

## Measured workload

`benchmarks/2026-09-06-exact-local.json` records seven warmed in-process runs of
47 local alleles, 47 hypotheses, three assays and 78 retained products. The median
was approximately 5.5 ms in this environment. It excludes process startup,
hypothesis expansion and I/O. This is a synthetic computational measurement, not
biological accuracy, a competitor comparison or a production throughput promise.
The older geometry-only benchmark measures a different model and workload.

## Remaining distribution and scientific gates

The connector can update existing repositories but cannot create
`dnncha/editwitness`. That standalone repository and package-index publication
remain unperformed unless an independently verified later record says otherwise.
The checked `scripts/publish_github.py` can create fresh standalone history with
an already authenticated maintainer CLI; it refuses an existing target.

Next: independent review of both response models, a benchmark with exact assay
metadata and independent measurements, then a real caller adapter and facility
pilots. See the audit, migration notes, continuation guide and machine roadmap.
