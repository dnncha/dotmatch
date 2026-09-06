# Build and verification status

Release: **0.1.0a1 research alpha**. Recorded 2026-09-05 UTC.

## Verified GitHub checks

[Completed verification run](https://github.com/dnncha/dotmatch/actions/runs/33998897914)

Verified source commit: `377c3bfafeb722bce3eac91d339dfec9d1028e68`.
Final delivery documentation records this already-tested code; no scientific
implementation or tests were changed after that run.

| Check | Observed result |
|---|---|
| Linux, macOS and Windows; Python 3.11 and 3.13 | All six test jobs passed. |
| Test suite | 572 passed, including the covered run. |
| Strict mypy check | Passed for `src/editwitness`. |
| Source hygiene | Passed on all six matrix jobs; this is not Ruff. |
| Schema snapshots | Passed on all six matrix jobs. |
| Wheel and source builds | Passed on all six matrix jobs. |
| Branch-aware coverage | 732/733 statements and 218/222 branches; 99.48% combined. |
| Source inventory | Verified in the delivery job before running checks. |

The initial CI attempt found a schema-registry type annotation and a Windows
locale-dependent test read. Both were corrected before the passing run. Neither
fix changes the observation model. Installed developer tools are not necessarily
executed checks: Ruff was not run and is not claimed as passed.

## Additional local checks

The implementation passed 572 local tests before the two portability fixes, plus
a full branch-aware coverage run. The local environment was Linux x86_64,
Python 3.13.5, Pydantic 2.13.4 and pytest 9.0.2. The remote covered run used
Python 3.13.15, Pydantic 2.13.5, pytest 9.1.1 and coverage 7.16.0.

The wheel was installed outside the source tree, and the CLI, bundled demo,
checksums and replay were exercised. Desktop and mobile HTML rendering was
inspected in Chromium. These are software checks, not experimental validation.

Testing includes exhaustive small-event comparisons against an independent
labeled-base sequence oracle, 400 seeded compound-edit cases, and 150 randomized
small panel-selection problems compared with independent subset enumeration.
It also covers malformed and oversized inputs, coordinate boundaries, original
primer integrity, paired-end read gaps, output protection and HTML escaping.

## Measured computational benchmark

The local streaming scanner evaluated **325,250 valid deletions** from a grid of
450,000 endpoint pairs in a median **0.135 seconds** over five warmed runs on this
environment. The workload and measurements are in
`benchmarks/2026-09-05-linux.json`. This is neither a comparison with another
package nor biological sensitivity or whole-genome throughput.

## Delivery and remaining gates

The complete source is on the isolated GitHub branch
`dnncha/dotmatch:editwitness/research-alpha-20260906`. **Do not merge this branch
into DotMatch.** DotMatch main was not changed by this delivery.

The standalone `dnncha/editwitness` repository has not been created: the connected
GitHub tool cannot create repositories. With an already authenticated local
GitHub CLI, `python scripts/publish_github.py --public` creates it from the
checked source with fresh history. The publisher refuses an existing target.
Its dry-run was checked; an actual new-repository publication was not performed.

Not completed: PyPI/Bioconda publication, package namespace reservation, hosted
documentation, a registered DOI, independent scientific review, an adjudicated
biological benchmark, or laboratory adoption. No clinical suitability, safety
certification or empirical assay accuracy is claimed.

See `docs/continuation.md` and `roadmap.json` for bounded next tasks.
