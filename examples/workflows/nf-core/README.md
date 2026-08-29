# nf-core Module Preparation (for upstream contribution)

This directory contains **nf-core-style module candidates** for DotMatch, ready for
contribution to the official [nf-core/modules](https://github.com/nf-core/modules) repository.

## Current Modules

- `dotmatch/count`: Runs the general `dotmatch count` command for fixed-window
  known-target counting across guides, barcodes, primers, or panels. Emits counts,
  `summary.json`, `sample_qc.tsv`, and per-read assignments for downstream QC.

- `dotmatch/demux`: Runs `dotmatch demux` for inline barcode demultiplexing.
  Emits split FASTQ outputs, demux summary JSON, assignments, and versions.

- `dotmatch/audit`: Runs `dotmatch audit` before correction-based assignment.
  Emits audit summaries, target-safety tables, and collision-pair diagnostics.

- `dotmatch/panel_check`: Runs `dotmatch panel check` for barcode-panel safety
  certification. Emits the panel certificate and safety/collision tables.

- `dotmatch/crispr_count`: Runs the native `dotmatch crispr-count` (or `count`) for
  MAGeCK-compatible CRISPR guide counting with full QC (`sample_qc.tsv`, summary).
  Supports Hamming/Levenshtein, k=0-3, auto offset detection, etc.

- `dotmatch/assay_run`: Runs a full `dotmatch assay run` (using AssaySpec TOML).
  Emits the complete set of workflow artifacts: reports, manifests, QC (including
  CRISPR-specific), counts, and versions. Ideal for complex or multi-step assays.

All module candidates:
- Use `process_medium` label and flexible containers (Docker + Singularity).
- Forward `task.cpus` / `task.ext.cpus` to `--threads` (leverages DotMatch's
  auto CPU detection where the subcommand supports it).
- Pass through `task.ext.args` for full CLI flexibility.
- Emit `versions.yml` for nf-core version tracking.
- Include `when` directive and nf-test stubs using shared fixtures.

## nf-core Upstream Contribution Prep

To turn these into reviewed nf-core modules:

1. **Fork & PR to nf-core/modules**
   - Copy `modules/local/dotmatch/*` (without the `local/` prefix) into
     `modules/nf-core/dotmatch/<subtool>/` in your fork of nf-core/modules.
   - Follow the [nf-core module template](https://nf-co.re/docs/contributing/modules).
   - Run `nf-core modules create` if starting fresh, then adapt.

2. **Lint & Test**
   - `nf-core modules lint dotmatch/crispr_count --dir /path/to/modules`
   - Run `nf-test` on the module tests (adapt paths in `tests/main.nf.test`).
   - Add to the main nf-core CI (they have automated linting + Docker builds).

3. **Prepare for Upstream Review**
   - Use exact bioconda/singularity container hashes from a released version
     (replace the `0.1.9--<bioconda_build>` template after the accepted
     Bioconda release is visible in BioContainers).
   - The upstream tree already includes maintainers, license metadata, test
     cases, and self-contained `tests/data/` fixtures.
   - Add more nf-test cases (different k, metrics, full vs stub runs) if needed.
   - Support additional common params via `task.ext` (e.g. `--auto-offset`,
     `--max-correction-qual`, `--format`).
   - Ensure all outputs match the [DotMatch schemas](https://github.com/dnncha/dotmatch/blob/main/docs/schemas.md)
     for MultiQC custom content + downstream tools.
   - Add a `dotmatch` top-level "umbrella" module if desired for the raw CLI.

4. **Pipeline Usage Example**
   In an nf-core pipeline:
   ```nextflow
   include { DOTMATCH_CRISPR_COUNT } from '../modules/nf-core/dotmatch/crispr_count/main'
   // ...
   DOTMATCH_CRISPR_COUNT (
       ch_reads_with_meta,
       ch_library,
       params.guide_start,
       params.guide_length,
       params.k,
       params.metric
   )
   ch_counts = DOTMATCH_CRISPR_COUNT.out.counts
   ```

5. **Evidence & Claims**
   - Any new nf-core usage should be recorded in `docs/workflow-adoption.json`
     and pass `make workflow-adoption-status`.
   - Keep module behavior aligned with DotMatch's [scientific-claims.md] and
     evidence gates (no broadening of supported claims without new data).

## Local Testing (in this repo)

From the repo root:
```bash
# Requires nextflow + nf-test in PATH
cd examples/workflows/nf-core/modules/local/dotmatch/count
nf-test test tests/main.nf.test
```

## End-to-End Pipeline Example

A minimal working pipeline is included under `pipeline/`:

```bash
cd examples/workflows/nf-core/pipeline
nextflow run main.nf --outdir results
```

It exercises the CRISPR and AssaySpec modules against the shared fixtures and
produces the expected artifacts. This gives pipeline maintainers a concrete
usage example and a basis for an upstream nf-core subworkflow or module example.

See [Workflow Submission Pack](../../../docs/workflow-submissions.md) for the
external submission checklist and tracking rules.

**Status**: The `upstream/` tree is the self-contained DotMatch 0.1.9 payload
for an nf-core/modules PR after exact BioContainers tags replace the local
templates.

The `modules/local/` versions are the maintained source for this repo's internal checks (check_workflow_examples.py) and examples.

Current state after prep:
- nf-core style container conditional + resource labels (cpus/memory/time)
- when directive + full ext.args + threads support (auto-detect)
- versions.yml + nf-test (including stub)
- meta.yml with authors/maintainers/license
- Self-contained tests/data/ in upstream/
- Local module candidates for `count`, `demux`, `audit`, `panel_check`,
  `crispr_count`, and `assay_run`
- Pipeline demo in `pipeline/`
- Comprehensive contribution guide

The main remaining work is the actual PR + linting in the nf-core org (see upstream/README.md for exact copy steps). Once merged, record in docs/workflow-adoption.json to unlock "ready" status.

Before opening an upstream PR:
- Pin a specific DotMatch version/container that has passed the release and
  distribution gates.
- Align any changes with current evidence boundaries in this repository.
- Update the DotMatch changelog and `docs/workflow-adoption.json` only after
  the external PR is merged or an official module page is public.
