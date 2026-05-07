# Changelog

All notable user-facing changes are tracked here. Claims in release notes must stay aligned with `docs/scientific-claims.md`.

## 0.1.0 - Initial GitHub Launch

### Added

- Native C short-DNA edit-distance and threshold assignment core.
- `dotmatch` CLI for pairwise distance, FASTQ assignment, demultiplexing, BCL milestone demultiplexing, count tables, CRISPR counting, audit, unmatched-read inspection, and validation.
- Python `dotmatch` package with ctypes bindings and local/GitHub wheel builds that bundle the native core.
- Deterministic assignment statuses: `unique`, `ambiguous`, `none`, and `invalid`.
- MAGeCK-compatible CRISPR count output, QC summaries, self-contained HTML reports, and audit artifacts.
- Reproducible benchmark reports, raw CSV evidence, and strict CRISPR claim gates.
- GitHub Actions CI, release artifact workflow, publication-readiness checker, contribution guide, security policy, support policy, citation metadata, and open-core boundary documentation.

### Verified Claims

- Known-target short-DNA assignment and CRISPR guide-counting claims are supported only where `make public-crispr-claim-gate` and `make crispr-sota-gate` pass on committed evidence.
- General alignment, barcode state-of-the-art, and raw BCL/CBCL state-of-the-art claims remain blocked until their dedicated gates pass.

### Packaging Status

- Source builds, local Python package builds, and GitHub release wheel/sdist artifacts are supported by repository checks.
- PyPI manylinux/musllinux Linux wheels, Bioconda packaging, Docker image publication, and Zenodo DOI registration are post-launch distribution tasks.
