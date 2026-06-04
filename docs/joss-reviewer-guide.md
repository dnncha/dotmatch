# JOSS Reviewer Guide

This page is a compact route through DotMatch for JOSS reviewers. It focuses on
installability, core functionality, test coverage, and evidence boundaries.

## Repository And Paper

- Repository: <https://github.com/dnncha/dotmatch>
- Paper: `paper/paper.md`
- Bibliography: `paper/paper.bib`
- Software archive DOI: <https://doi.org/10.5281/zenodo.20541629>
- License: Apache-2.0 in `LICENSE`
- Submission tracker: <https://github.com/dnncha/dotmatch/issues/43>

The paper is intentionally short. Detailed command behavior, schemas, and
benchmark evidence are documented in the repository rather than repeated in the
paper.

## Quick Install Check

From a clean checkout on Linux or macOS:

```bash
git clone https://github.com/dnncha/dotmatch.git
cd dotmatch
make

./dotmatch --version
./dotmatch dist ACGT AGGT
./dotmatch leq 1 ACGT AGGT
```

Python package smoke test from the checkout:

```bash
python3 -m pip install .
python3 -c "import dotmatch; print(dotmatch.distance('ACGT', 'AGGT'))"
```

The source build requires a C compiler, `make`, Python 3.9 or newer for the
Python package, and zlib for FASTQ.gz support.

## Core Functionality Checks

Run the small native and CLI checks first:

```bash
make test
make cli-test
```

Run the Python package and wrapper tests:

```bash
make python-test
make python-package-test
```

Expected result: all targets complete successfully. These are the same test
families exercised by the GitHub Actions `ci` workflow.

## Paper Checks

The local JOSS paper gate verifies the manuscript metadata, required sections,
bibliography keys, archive DOI, AI disclosure, and 750-1750 word range:

```bash
make joss-paper-ready
```

The GitHub Actions workflow `JOSS draft PDF` compiles `paper/paper.md` with the
Open Journals draft PDF action and uploads a `joss-paper` artifact containing
`paper.pdf`.

## Claim Verification

DotMatch keeps scientific claims evidence-bounded. Reviewers should treat
`docs/scientific-claims.md` as the claim ledger and `docs/trust-and-scope.md`
as the high-level scope statement.

Useful evidence gates:

```bash
make scientific-readiness-ready
make assay-evidence-ready
make citation-metadata-ready
make release-ready
```

`make release-ready` is intentionally broad and includes Python tests, package
checks, documentation, evidence gates, citation metadata, workflow examples,
release metadata, and distribution-record checks.

## Representative Examples

For a first run, use:

- `docs/getting-started.md`;
- `docs/tutorials/crispr-count-first-run.md`;
- `examples/crispr_guides/README.md`.

For evidence-backed scenarios, use:

- `docs/evidence-gallery/README.md`;
- `docs/benchmarks/public_crispr/README.md`;
- `docs/benchmarks/barcode_demux/README.md`;
- `docs/benchmarks/feature_barcode/README.md`;
- `docs/benchmarks/amplicon_panel/README.md`;
- `docs/benchmarks/oligo_adapter/README.md`.

## Scope Boundaries

The paper and documentation should not be read as claiming that DotMatch is:

- a genome aligner;
- a basecaller;
- a production Illumina BCL converter;
- an adapter trimmer;
- a UMI/cell quantification pipeline;
- an amplicon consensus or variant-calling tool;
- a downstream CRISPR phenotype statistics package.

DotMatch is a deterministic known-target short-DNA assignment layer. Its
outputs preserve `unique`, `ambiguous`, `none`, and `invalid` states so workflow
authors can avoid silent forced assignments.

## Reviewer Issues

Reviewer-raised software, documentation, or reproducibility requests should be
opened as repository issues when practical. After review changes are accepted,
the project should create a fresh tagged release and archive that final reviewed
state on Zenodo before posting the final version and DOI in the JOSS review
issue.
