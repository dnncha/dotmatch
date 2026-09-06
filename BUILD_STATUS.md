# EditWitness 0.2.0a2 — verification and distribution boundary

**6 September 2026. Research alpha; not empirically biologically validated.**

The current public release target is the namespaced GitHub prerelease
`dnncha/dotmatch: editwitness-v0.2.0a2`. EditWitness is an independent package.
DotMatch main must not be changed or receive this branch as a merge. A standalone
`dnncha/editwitness` repository and PyPI publication are separate remaining steps.

## Evidence available before the remote release gate

The local suite passed **630 tests** after the hardening changes. Schema
snapshots and Python 3.11 syntax/text checks passed. New tests cover explicit
heterozygous haplotype choice, evidence budgets and the offline software self-test.
**Executed local branch coverage:** 1,046 / 1,049 statements and 322 / 328
branches; the configured 95% combined gate passed, including CLI subprocesses.
Historical 0.1.0a1 and 0.2.0a1 records are retained under `docs/history/`; they
are not reused as proof for changed code.

Local dependency downloads were unavailable. The remote CI is responsible for
strict mypy, Ruff, clean dependency installation and the expanded operating-system
and Python matrix. A configured job is not a passed job.

## Public release gate

The temporary transport first commits the complete reviewed source to the
isolated branch. CI checks out that exact committed SHA, runs all gates, builds
the wheel and source archive, installs the wheel outside the checkout, and runs
the tests from the extracted source archive. Only successful gates allow a
namespaced prerelease. Publication must not mark it as DotMatch's latest release.

The release assets include `release-evidence.json`, with the exact source commit,
source inventory digest, executed self-test and CI run URL. That run is the
source of truth for passed or failed jobs. `SHA256SUMS` covers the wheel and source
distribution; `EVIDENCE_SHA256SUMS` covers the supplementary report and provenance.
A checksum is not a signature or a claim of biological validity.

Do not report a GitHub prerelease as published until its API record is public,
not draft, and its assets and target commit match the executed gate. Do not
report a PyPI release: none is part of this workflow.

## Remaining scientific gate

Independent genome-engineering review and an adjudicated, provenance-complete
biological benchmark remain outstanding. Software checks do not establish PCR
sensitivity, allele dosage, safety, clinical validity, novelty or real adoption.
