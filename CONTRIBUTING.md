# Contributing to DotMatch

DotMatch is an open scientific core for deterministic short-DNA known-target assignment. Contributions are welcome when they improve correctness, reproducibility, performance, documentation, or usability without weakening the claim boundary.

## Good First Contributions

- Improve examples, error messages, and documentation for real bioinformatics workflows.
- Add small deterministic tests for CLI behavior, file schemas, ambiguity handling, and edge cases.
- Add benchmark comparators when their semantics are documented and reproducible.
- Improve packaging for source installs, Docker, PyPI wheels, or Bioconda.
- Tighten validation, audit, or report outputs that help scientists trust a run.

## Development Setup

```bash
make
make shared
make test
make cli-test
DOTMATCH_LIB=$PWD/libdotmatch.dylib PYTHONPATH=$PWD/python python3 -m pytest python/tests
make python-package-test
make publication-ready
```

On Linux, use `DOTMATCH_LIB=$PWD/libdotmatch.so` for Python tests. The CI workflow runs native build/test, CLI fixtures, Python tests, sanitizers on Linux, coverage on Linux, and shared-library builds.

## Claim Discipline

Performance and state-of-the-art claims need evidence in the repository. A claim is ready only when:

- the exact command is documented;
- the raw CSV artifact is checked in under `benchmarks/raw/`;
- the report or figure can be regenerated from that raw artifact;
- correctness validation passes against the relevant oracle;
- comparator semantics are stated plainly;
- the matching gate script passes when one exists.

Do not broaden claims from one workflow to another. In particular, DotMatch is not currently a genome aligner, and barcode/BCL SOTA claims are blocked until their strict gates pass.

## Benchmark Artifacts

Commit reproducible evidence, not large working datasets. Keep curated CSVs under `benchmarks/raw/`, figures under `benchmarks/figures/`, and reports under `docs/benchmarks/`. Downloaded FASTQ/BCL data, generated demultiplexed FASTQs, competitor environments, and scratch outputs should stay in ignored paths.

## Pull Requests

Before opening a PR:

```bash
make test
make cli-test
make python-test
make python-package-test
make publication-ready
```

Run `make coverage` for changes that touch `src/qdalign.c`, `src/qda.c`, CLI parsing, assignment semantics, or report generation. Run the relevant claim gate when changing benchmark artifacts or wording.

## Contributor Certification

By contributing, you certify that you wrote the contribution or have the right to submit it under the Apache-2.0 license used by DotMatch. Add a Developer Certificate of Origin sign-off to commits when practical:

```bash
git commit -s -m "your message"
```
