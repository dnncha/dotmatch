# nf-core/modules submission payload

This directory is a self-contained candidate payload for DotMatch modules in
`nf-core/modules`. It is prepared against DotMatch 0.6.3 and uses the immutable
GHCR manifest `sha256:c43dd55c5c58d4b689af8e76a4d54e19973af52f133f0f77781f6949736708e0`. The release workflow verified that manifest
for both `linux/amd64` and `linux/arm64`; module pins use `docker://` for
direct Singularity/Apptainer pulls.

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
module. Keep the verified 0.6.3 GHCR manifest digest pin. Replace it only with a
reviewed public Bioconda/BioContainers build after its installed workflow checks
pass.

After an upstream pull request is accepted, add its public URL to
`docs/workflow-adoption.json` and run `make workflow-adoption-status`.
