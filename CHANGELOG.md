# Changelog

All notable user-facing changes are tracked here. Public statements in release notes must stay aligned with `docs/scientific-claims.md`.

## 0.1.3 - 2026-05-22

### Added

- Assay reliability evidence now records biological units, unsupported claims,
  and minimum public evidence requirements for each public claim boundary.
- Release-gate coverage now keeps BCL and paired-combinatorial evidence lanes
  from disappearing silently from `docs/assay-evidence.json`.

### Changed

- Bioconda packaging remains a native CLI and C library package. The recipe now
  passes `PKG_VERSION` directly into the native build and adds installed-file and
  tiny native `count` smoke tests without adding Python, Workbench, or browser
  dependencies.
- Package and release metadata are aligned for the next `0.1.3` update after
  the public `0.1.2` Bioconda publication.
- Public installation docs now distinguish the native Bioconda package from the
  Python workflow layer that provides `dotmatch assay`, barcode/panel
  convenience namespaces, and Workbench-backed AssaySpec workflows.

### Fixed

- Release and packaging docs no longer describe the initial Bioconda recipe as
  pending; `bioconda/bioconda-recipes#65367` has already published DotMatch
  `0.1.2` for the current Bioconda platforms.
- BCL evidence metadata no longer lists the broader `bcl-comparison-gate` as a
  passing release evidence gate; the supported public BCL statement remains the
  narrow tiny-BCL parser milestone.

## 0.1.2 - 2026-05-18

### Packaging

- Published the first Bioconda package for DotMatch with native CLI smoke tests.
  Bioconda availability is a distribution milestone only; it does not expand
  the scientific evidence boundaries documented in `docs/scientific-claims.md`.

## 0.1.0 - Initial Release

### Added

- Native C short-DNA edit-distance and threshold assignment core.
- `dotmatch` CLI for pairwise distance, FASTQ assignment, demultiplexing, BCL milestone demultiplexing, count tables, CRISPR counting, audit, unmatched-read inspection, and validation.
- Python `dotmatch` package with ctypes bindings and local/GitHub wheel builds that bundle the native core.
- Deterministic assignment statuses: `unique`, `ambiguous`, `none`, and `invalid`.
- MAGeCK-compatible CRISPR count output, QC summaries, self-contained HTML reports, and audit artifacts.
- Reproducible benchmark reports, raw CSV evidence, and strict CRISPR validation gates.
- GitHub Actions CI, release artifact workflow, repository-readiness checker, contribution guide, security policy, support policy, and citation metadata.

### Verified Evidence

- Known-target short-DNA assignment and CRISPR guide-counting statements are supported only where `make public-crispr-evidence-gate` and `make crispr-comparison-gate` pass on committed evidence.
- General alignment, barcode comparative, and raw BCL/CBCL comparative wording should stay within the evidence boundaries documented in `docs/scientific-claims.md`.

### Packaging Status

- Source builds, local Python package builds, and GitHub release wheel/sdist artifacts are supported by repository checks.
- PyPI manylinux/musllinux Linux wheels, Bioconda packaging, Docker image distribution, and Zenodo DOI registration are separate distribution tasks.
