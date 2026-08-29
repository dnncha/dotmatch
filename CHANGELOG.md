# Changelog

All notable user-facing changes are tracked here. Public statements in release notes must stay aligned with `docs/scientific-claims.md`.

## Unreleased

### Fixed

- Registered DotMatch MultiQC search patterns through the `before_config` hook
  so an installed plugin is discoverable before MultiQC indexes input files.
- Kept panel summaries out of the generic summary match and added fixtures for
  assignment summaries and top-unmatched diagnostics.

### Changed

- Added a public-record ecosystem status ledger and updated workflow
  submission records to distinguish open, accepted, released, and installable
  states.

## 0.2.2 - 2026-07-23

### Changed

- Documented the scoped nf-core CRISPR guide-counting module and its
  maintainer-review status, with repeatable Conda-profile module tests.
- Updated the public site dependencies to patched Next.js, React, React DOM,
  and Sharp releases.

## 0.2.1 - 2026-07-17

### Fixed

- Replaced repository-relative README links with stable documentation and
  project URLs so links work when the README is rendered on PyPI.
- Corrected the package homepage, documentation, changelog, and social-card
  metadata.

### Changed

- Rewrote the README, documentation landing page, package description, and
  website around the tasks DotMatch performs.
- Reduced the documentation navigation to runnable guides, workflow pages,
  APIs, file formats, benchmarks, packaging, and citation.
- Added checks for PyPI-safe README links, missing local documentation links,
  inconsistent public product names, and internal or inflated public wording.

## 0.2.0 - 2026-07-17

- Added a compatibility-safe AssayCode Bioconda metapackage template and a deterministic two-recipe release-handoff generator while preserving `dotmatch` as the engine package.

- Added a deterministic experimental panel simulator for pre-sequencing yield, ambiguity, no-call, confusion, and FDR analysis.
### AssayCode platform

- Added the additive `assaycode` CLI and Python namespace while preserving the `dotmatch` package, executable, ABI, DOI, schemas, and citation contract.
- Added AssayScript v2 compilation for R1/R2/I1/I2 segment declarations, fixed or anchored extraction, positional jitter, orientation, per-segment matching policies, allowed-combination tables, source and library fingerprints, bounded safety findings, deterministic strategy selection, and portable JSON plans.
- Added an experimental fail-closed AssayScript runtime for synchronized multi-read FASTQs, fixed or anchored extraction, ambiguity-preserving segment calls, allowed-tuple filtering and rescue, atomic assignment/count/event outputs, and fingerprinted provenance summaries.
- Added `assaycode compile`, `assaycode inspect`, and compatibility-safe assay workflow shortcuts.
- Added an explicitly experimental calibration module with per-cycle error fitting, Phred shrinkage, selective posterior calls, likelihood-ratio abstention, joint decoding over permitted tuples, Brier score, expected calibration error, held-out FDR threshold selection, and smoothed abundance priors.
- Added `assaycode watch`, a bounded-memory JSONL monitor that emits assignment-rate confidence intervals and threshold-based sequential QC decisions.
- Added a release-blocking AssayCode readiness gate, focused tests, updated scientific claim boundaries, a rewritten paper, and Bioconda smoke tests for both command identities.

### Added
- `dotmatch.tl` submodule: scverse/scanpy-style tools (`tl.assign_features`, `tl.feature_counts`, aliases for CRISPR/feature barcodes), with in-place or copy semantics and provenance in `.uns`.
- Pure, dependency-light parsers in `dotmatch.multiqc` (`parse_sample_qc_tsv`, `parse_crispr_qc_summary_tsv`, `parse_assay_manifest_summary_tsv`) usable from any Python code (notebooks, custom reports) while strictly following the documented schemas.
- Proper MultiQC plugin registration via entry point in `pyproject.toml` + improved `DotMatchModule` with full parsing of sample QC, CRISPR QC, and assay manifest (adds sections + general stats with scientific descriptions).
- R package skeleton (`R/`, `DESCRIPTION`, `NAMESPACE`, `vignettes/dotmatch.Rmd` + reticulate wrappers) for Bioconductor / tidyverse users.
- New test coverage for tl and parsers (`test_tl_and_integrations.py`).

### Changed
- Expanded optional extras (`[anndata]`, `[polars]`, `[multiqc]`) and top-level exposure.
- nf-core module meta.yml files now document threads/cpus for better resource awareness.
- README and proposals-and-roadmap.md updated with usage examples and status.

### Performance (best-of-n)
- All three candidates from `/best-of-n performance improvement` applied after confirming orthogonality (Myers edit paths vs hamming hot-path vs qda driver batch/IO), no conflicts, full builds/tests pass, and behavior-identical oracles (72/129bp myers vs DP; hamming k within counts match scalar expectation; large-batch count pipelines produce exact TSVs).
- Multi-word Myers: patterns >64bp (up to 512) now use fast generalized bit-parallel (portable carry, no __int128 dependency) instead of DP fallback. `qdaln_edit_distance*` and leq for k>=2 benefit for long guides/primers/amplicons.
- SIMD hamming: `same_length_hamming_distance_within_k` (and seed index) accelerated with AVX2 (x86) or NEON (arm64/aarch64) + scalar fallback. Direct hamming count paths (common for CRISPR k<=1 uniform) measurably faster.
- Batching/IO: raised to 1M reads/batch + `seq_buffer` reuse + `reset_seq_buffer` helper in count feeders (both single and read_threads>1 paths). Reduces malloc churn on large .fastq.gz while keeping fixed-block fastpath for uniform-length assays. RSS tradeoff documented in code.
- All changes respect evidence culture (no new public claims without raw CSV + gate); internal perf only. Verified: make test, make cli-test, 370 pytest (native), count-agreement, synthetic+real-oracle runs, A/B on mixed short/long + threaded count.

These changes add checked interfaces for scverse, MultiQC, R, and nf-core workflows while preserving unique-only counting, explicit ambiguity handling, schema-aligned parsing, and recorded provenance.

## 0.1.5 - 2026-05-26

### Added

- Added a GuideCounter-compatible `dotmatch guide-counter count` entrypoint that
  accepts GuideCounter-style count flags and writes counts, extended-counts, and
  stats outputs.
- Added Hamming `k=2` and `k=3` fixed-window guide-counting support, with exact
  audit fields that tell users whether those larger Hamming radii are safe for a
  target library.
- Added DotMatch-vs-Bowtie 1 Hamming `k=2`/`k=3` comparator evidence for the
  Sanson/Brunello public CRISPR lane. These rows are kept separate from
  GuideCounter-compatible `k=1` claims.

### Changed

- The top-level CLI help and counting help now make the supported Hamming
  radius range and Hamming `k=2`/`k=3` audit expectations visible from
  `--help`.
- Bioconda recipe readiness now tests the Python console-script package, native
  import discovery, workflow namespace help, Hamming `k=2` CRISPR counting,
  exact Hamming `k=3` audit output, GuideCounter-compatible output files,
  barcode inference, and panel design.
- The Bioconda recipe opts into `osx-arm64` Apple Silicon builds and relies on
  host `zlib` to export the linked `libzlib` runtime package instead of
  duplicating `zlib` in run requirements.
- Release and packaging metadata are aligned for the next `0.1.5` patch
  release, while public docs still point users to the latest verified Bioconda
  package until `0.1.5` propagates.

## 0.1.4 - 2026-05-23

### Changed

- Bioconda packaging now installs the Python `dotmatch` console script so
  `dotmatch assay`, `dotmatch barcode`, and `dotmatch panel` are available from
  the Bioconda package alongside the native assignment commands.
- Top-level Python help now lists the workflow namespaces, and native help no
  longer includes package-channel policy notes.
- Python native builds now honor Conda `CPPFLAGS` and `LDFLAGS` during wheel and
  source-install builds.

### Fixed

- `dotmatch assay --help` now describes the `check`, `plan`, and `run`
  subcommands.

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
