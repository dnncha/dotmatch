# DotMatch Proposals, Performance Ideas, and Bioinformatics Adoption Roadmap

This document collects concrete ideas for making DotMatch the default choice for known-target short-DNA assignment in academic cores, CROs, biotech, and pharma. It is intentionally forward-looking and not a source of current claims—see `docs/scientific-claims.md` and the gate scripts.

All performance or comparative statements must be backed by checked artifacts under `benchmarks/raw/`, regenerated figures, and passing gates before appearing in README, docs, or releases.

## 1. Performance Enhancements (High Impact, Low Risk)

- **Default auto-threads (done)**: `--threads 0` (and omitted) now detects `sysconf(_SC_NPROCESSORS_ONLN)` (capped). Transparent fallback to 1 for diagnostic outputs. Out-of-box wins on 8–128 core analysis nodes common in bioinfo.
- **Multi-word / wide Myers bit-parallel (implemented)**: Generalized `qdaln_edit_distance_myers64` (now handles nwords>1 for patterns up to 512bp) with portable carry propagation (no __int128). Routes >64bp Levenshtein (k>=2) and general distance to fast path; oracle-matched vs DP in tests + explicit 72/129bp checks. Part of best-of-n perf round.
- **SIMD hamming / popcount acceleration (implemented)**: `same_length_hamming_distance_within_k` + hot paths use AVX2 (256-bit xor+movemask pop) or NEON (aarch64) with scalar fallback. Applied in index seeding too (prefetch). Orthogonal to Myers. Verified bit-identical counts on hamming k=0/1/2 paths.
- **Allocator / batch tuning + IO reuse (implemented)**: BATCH 1M reads, `seq_buffer` static reuse + `reset_seq_buffer` (reclaims only var-len xstrndup, keeps fixed block for common uniform case). Used in both direct count and threaded read paths. Throughput win on large gz with modest RSS. All three landed cleanly after best-of-n (N=3 candidates in isolated worktrees, evaluated for correctness/safety/no overlap, all applied).
- **IO and decompress**: Optional libdeflate for gzip (faster on Linux); parallel decompress for multi-FASTQ inputs (current is per-sample sequential read + parallel compute). Memory-map for uncompressed.
- **Index for higher k / larger windows**: Safer k=2/3 Hamming already present for guides; extend seed tables or hierarchical (e.g., 2-seed + verify) for Levenshtein k=2 on 32+ bp without exploding memory. Audit tables for k>1 already gated.
- **GPU production path (Linux)**: Move the Metal brute-force experiment to a CUDA or OpenCL/Vulkan path with proper batching, CPU fallback for N/IUPAC, and end-to-end gates on real feature-barcode + CRISPR + BCL extracts. Only promote after evidence.
- **Profiling harness**: Add `make bench-profile` (perf or Instruments) and a hotspot report script that feeds the roadmap.

Target: keep or improve the current "hundreds of k reads/sec per core" while handling 100M–1B read production lanes with < few GB RSS.

## 2. Python / Data Science Adoption (Immediate Wins)

- **pandas / polars interop (implemented + expanded)**: `targets_from_dataframe`, `results_to_dataframe`, `assign_dataframe` with column heuristics + polars support (via `_to_pandas` converter). Exposed under `dotmatch[pandas]` and `dotmatch[polars]` extras.
- **scverse / AnnData + tl (implemented + UX hardened)**: Full `dotmatch.tl` submodule (`tl.assign_features`, `tl.feature_counts`, aliases for crispr/feature). In-place or copy, adds provenance. Plus the anndata bridges. `dotmatch[anndata]`. Pure parsers in multiqc for notebook use too.
- **MultiQC (real module + pure parsers)**: Full DotMatchModule + public `parse_*_tsv` functions (no multiqc dep) for direct use in reports/notebooks. Tied exactly to schemas.md for scientific fidelity. `dotmatch[multiqc]` extra.
- **nf-core**: Modules now document threads in meta.yml + pass cpus.
- **Streaming / generator API (implemented)**: `dotmatch.stream_assign(fastq_path, targets, ...)` yields FASTQ-order row assignments without loading full read lists in the Python layer. Helpers now cover `load_targets`, `iter_fastq`, `assignment_summary`, and `write_assignments_tsv` for notebook and workflow handoff.

These let data scientists stay in notebooks without shelling out or writing glue.

## 3. Workflow & Pipeline Penetration

- **nf-core**: Current examples are local modules. Goal: contribute official `nf-core/modules` for `DOTMATCH_CRISPR_COUNT`, `DOTMATCH_ASSAY_RUN`, barcode demux, etc. Requires:
  - Full `meta.yml`, `tests/`, container directive (bioconda or multi-arch), resource labels that scale with `task.ext.args`.
  - Example full pipeline under `examples/workflows/nf-core/pipeline/`.
  - Update `docs/workflow-adoption.json` and run the gate once a public nf-core PR merges or a high-profile pipeline uses it.
- **MultiQC module**: Beyond the custom-content example, ship a small `multiqc_dotmatch` Python package (or contrib to MultiQC) that auto-discovers `sample_qc.tsv`, `crispr_qc.*`, `assay_manifest.*`, unmatched summaries, and renders nice tables/plots + links to the HTML reports. pip-installable, auto in nf-core MultiQC.
- **Snakemake / Nextflow / Galaxy polish**: Make the existing examples "copy-paste ready" with test data + CI that runs them. Add a small nf-core-like test profile.
- **HPC modules / EasyBuild / Spack**: Recipes or docs so institutional clusters can `module load dotmatch`.
- **Docker / Singularity / Apptainer**: Official multi-arch images on GHCR (already in Dockerfile + release?); document for cloud batch (AWS, GCP, Azure Life Sciences).

## 4. Assay & Domain Features

- **Combinatorial / dual / pooled barcodes** (Perturb-seq style): First-class support for left+right independent or linked assignment with joint ambiguity rules and collision simulation in `panel` commands.
- **Quality-aware & position-specific models**: Extend `--max-correction-qual` rescue to full matrices; allow per-base quality masks or "do not correct if Q< in this window".
- **Longer / variable windows + anchors**: Better support for primer+amplicon starts, adapter prefixes with partial overlap, sliding or multi-window per read.
- **BCL / basecalls full**: Expand the classic lane-1 milestone to multi-lane, CBCL (NovaSeq), index-hopping diagnostics, and direct comparison to bcl-convert outputs where semantics overlap.
- **Panel design / audit enhancements**: k>2 enumeration (with cost controls), GC/ homopolymer / synthesis-constraint models that feed directly into `panel check` safety, plate randomizers that also optimize for ambiguity spheres.
- **Unmatched / failure mode taxonomy**: Expand "autopsy" to general assays (not just barcodes); ML-assisted but auditable clustering of high-count unmatched for novel biology or contamination.

Keep scope: stop at assignment + QC + audit. Leave UMI dedup, cell calling, hit stats, variant interpretation to established tools.

## 5. Language & Interop Bindings

- **R / Bioconductor**: Thin wrapper package (in `R/`) using reticulate for the Python API + S4 friendly outputs. Includes Matcher wrapper, assign_features_df, etc. Vignette with examples. Ready for Bioc/CRAN. (skeleton + polish done).
- **Julia**: `DotMatch.jl` calling the C lib via `ccall` (fast) or PythonCall.
- **C++ / Rust headers + crates**: The existing `include/qdalign.h` + libs are already usable; add pkg-config, CMake, or cargo example.

## 6. UX, Reports, Trust & Commercial Readiness

- **Rich reports**: Self-contained HTML + embedded interactive tables (sortable, filterable unmatched + per-target error spectra). PDF export for lab notebooks / audits.
- **Workbench (desktop)**: Continue Tauri work; add "load assay spec + run locally + review QC" flows, one-click "export to nf-core / MultiQC / AnnData".
- **Error messages & wet-lab empathy**: Common failure modes (offset by 1, RC barcode in pool, homopolymer synthesis slips, N in target vs read) should produce actionable suggestions + links to autopsy / panel design.
- **Provenance & reproducibility**: Every summary.json / report already contains command, versions, commit-ish, phase times, read_threads. Make it machine-verifiable against a signed manifest for regulated work (GLP?).
- **Support / SLAs**: Clear "community" vs "enterprise support" path (even if just a paid triage contract or training offering) for core facilities that need guaranteed response.
- **Licensing clarity**: Apache-2.0 is excellent for penetration. Consider additional patent grant language if any future IP appears.

## 7. Evidence, Claims, Marketing & Distribution

- **More public lanes**: Add 1–2 additional large public CRISPR (beyond Sanson/Brunello), a 10x feature-barcode or guide-capture, a real GBS or RAD-seq barcode set, a clinical amplicon panel. Each needs: fetch script, expected, comparator (MAGeCK/guide-counter/cutadapt/bowtie where fair), Edlib/scan validation, repeated runs, gate.
- **Competitor refresh**: Re-bench against latest guide-counter, MAGeCK, cutadapt, any new "ultra-fast demux" tools. Keep semantics documented.
- **Case studies**: 1–2 anonymized or public "we replaced X with DotMatch in production and got Y% faster + better QC visibility" from cores or companies.
- **JOSS / paper / blog**: Keep the JOSS process honest; add a methods-focused preprint or blog series on the indexing + ambiguity design.
- **Citation flywheel**: Prioritize official nf-core modules, a released MultiQC module, Galaxy and Snakemake wrappers, generated `methods.md` / `CITATION.bib` / `software_versions.yml` artifacts, copyable citation UI in Workbench and HTML reports, documented external pilots, and a public used-by page. Track the executable checklist in `docs/citation-flywheel.md`.
- **Distribution**: Ensure every release has:
  - Bioconda (already), PyPI wheels (linux repaired + mac), source sdist.
  - Containers published and smoke-tested in the release gate.
  - Updated `docs/distribution-release.json`.
- **Discoverability**: SEO-friendly docs, schema.org for tools, bio.tools registration, GA4GH? (stretch), prominent "used by" once real users appear.

## 8. Process & Contribution

- Keep the evidence discipline: no wording creep.
- Add "good first issue" labels around pandas polish, doc examples, small workflow fixtures.
- Automate more of the release gates (some already in `make *`).
- Consider a small "DotMatch Users" Slack/Discord or GitHub discussion for feedback without support burden.

## Prioritization Sketch (as of 2026)

1. Land pandas + AnnData + tl + MultiQC (full parsers + entrypoint) + nf-core + R (skeleton + polish: more functions, Matcher wrapper, expanded vignette) (done).
2. Auto-threads UX/perf (done).
3. Perf round: multi-word Myers (portable) + SIMD hamming (AVX2/NEON) + 1M batch/seq_buffer reuse (all 3 best-of-n candidates evaluated in isolated worktrees for correctness, safety, no overlap; all applied after verification; see "Performance Enhancements" section). (done).
4. nf-core module upstreaming (local + upstream/ self-contained PR payload + playbook ready; external PR + workflow-adoption.json update pending to flip gate). More public evidence lanes.
4. One additional large public evidence lane + gate (to safely broaden claims).
5. R package + vignette (done; in tree + docs).
6. nf-core upstream PR + record in workflow-adoption.json (local+upstream payload + playbook ready in examples/workflows/nf-core/upstream/).
7. Production GPU lane (or decide not to for v1 scope).
8. BCL expansion if demand appears.

Run `make repository-ready assay-evidence-ready scientific-readiness-ready workflow-adoption-status` (and the specific comparison gates) after material changes.

## How to Propose / Implement

- Open issue with "proposal:" or PR with small self-contained change + test + (if perf/claim) raw data + gate update.
- For large items, add a design note under `docs/superpowers/` or a design/ subdir.
- Evidence first: if it affects numbers or "faster/safer", the CSV + script + gate must land in same or preceding PR.

This roadmap is aspirational. The current strength—deterministic, auditable, fast-enough assignment with explicit ambiguity and strong validation gates—is already a solid foundation for industry use in controlled-target assays. Enhancements should amplify that without diluting the contract.

## References / Related

- `docs/benchmarks/*/README.md`
- `docs/scientific-claims.md`
- `docs/assay-evidence.json`
- `scripts/check_*.py`
- `CONTRIBUTING.md`
