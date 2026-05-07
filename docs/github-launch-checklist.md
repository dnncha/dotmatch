# GitHub Launch Checklist

This checklist covers repository settings and external launch work that cannot be fully captured by code tests.

## Repository Settings

- Description: `Fast exact short-DNA known-target assignment for CRISPR guides, barcodes, primers, panels, and whitelists.`
- Website: link to the hosted DotMatch site or the README until a site is deployed.
- Topics:
  - `bioinformatics`
  - `computational-biology`
  - `crispr`
  - `fastq`
  - `barcode-demultiplexing`
  - `edit-distance`
  - `sequence-analysis`
  - `genomics`
  - `c`
  - `python`
- Enable Issues, Discussions, Dependabot alerts, the committed CodeQL workflow, and private vulnerability reporting.
- Add a social preview image before public launch.

## First Public Release

- Tag `v0.1.0` only after `make test`, `make cli-test`, `make python-test`, `make python-package-test`, `make publication-ready`, `make coverage`, `make public-crispr-claim-gate`, and `make crispr-sota-gate` pass.
- Release notes should lead with the narrow known-target assignment claim, not broad alignment language.
- Use the tag-driven release workflow described in `docs/release-process.md`.
- Attach or link the raw benchmark CSVs, benchmark reports, and `CITATION.cff`.
- Confirm `codemeta.json` and `docs/methods-and-citation.md` match the release version.
- Create a Zenodo archive and add the DOI to `CITATION.cff` after the DOI exists.

## Package Channels

- GitHub source release: first launch channel.
- Docker: publish an image once tags are stable.
- PyPI: publish the sdist first; publish Linux binary wheels only after manylinux/musllinux wheel repair is in place.
- Bioconda: open a recipe once the CLI install and smoke test are stable.

## Scientific Outreach

- Submit the repository to relevant software registries after release: bio.tools, SciCrunch/RRID if appropriate, and GitHub topic collections.
- Prepare a short reproducibility note that points to `docs/scientific-claims.md`, `docs/benchmarks/`, and the exact gate commands.
- Keep the manuscript/application-note claim aligned with the strictest passing gate.

## Blocked Before Broader Claims

- Barcode SOTA: requires `make barcode-sota-gate` to pass with claim-grade fixed-length real barcode sheet metadata, repeated real-data rows, Cutadapt rows, and a second comparator. The exact hash splitter counts only for the `k=0` exact-prefix lane.
- Raw BCL SOTA: requires `make bcl-sota-gate` to pass with CBCL evidence and production demux comparator validation.
- General alignment: blocked until DotMatch has reference indexing, mapping semantics, traceback/CIGAR, and appropriate aligner benchmarks.
- Production Python package: release only for platforms whose wheels pass the clean-venv native import verifier.
