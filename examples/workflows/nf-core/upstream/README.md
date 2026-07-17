# nf-core/modules submission payload

This directory is a self-contained candidate payload for DotMatch modules in
`nf-core/modules`. It is prepared against DotMatch 0.1.9 and uses the immutable
BioContainers build `0.1.9--py311h13f8228_0`, which is available from Quay and
the Galaxy Singularity depot.

## Included modules

- `dotmatch/count`
- `dotmatch/demux`
- `dotmatch/audit`
- `dotmatch/panel_check`
- `dotmatch/crispr_count`
- `dotmatch/assay_run`

Each module includes `main.nf`, `meta.yml`, an nf-test definition, and the small
fixtures required by its test. The payload preserves DotMatch's `unique`,
`ambiguous`, `none`, and `invalid` assignment outcomes and exposes
`task.ext.args` for module-specific options.

## Local verification

From the DotMatch repository root, run:

```bash
make workflow-examples-ready
make workflow-integration-test
make reviewer-readiness-ready
```

The integration test requires the workflow tools named by the target. Check
that `nf-test`, Nextflow, Snakemake, Planemo, and MultiQC are available before
interpreting a missing-command failure as a module failure.

## Upstream verification

Copy `modules/nf-core/dotmatch/` into a current `nf-core/modules` checkout and
run the repository's formatter, module lint, and nf-test commands for each
module. Keep the exact container tag unless a later DotMatch release has passed
the distribution checks and is present in both Quay and the Galaxy Singularity
depot.

After an upstream pull request is accepted, add its public URL to
`docs/workflow-adoption.json` and run `make workflow-adoption-status`.
