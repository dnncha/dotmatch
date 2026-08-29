# DotMatch nf-core/modules Upstream Payload

This directory contains the nf-core/modules candidate tree prepared for
DotMatch 0.1.9. It can be copied into a fork of
https://github.com/nf-core/modules after the container placeholders are replaced
with exact BioContainers tags for the submitted release.

## Module Tree

- `modules/nf-core/dotmatch/count/main.nf`
- `modules/nf-core/dotmatch/count/meta.yml`
- `modules/nf-core/dotmatch/count/tests/main.nf.test`
- `modules/nf-core/dotmatch/demux/main.nf`
- `modules/nf-core/dotmatch/demux/meta.yml`
- `modules/nf-core/dotmatch/demux/tests/main.nf.test`
- `modules/nf-core/dotmatch/audit/main.nf`
- `modules/nf-core/dotmatch/audit/meta.yml`
- `modules/nf-core/dotmatch/audit/tests/main.nf.test`
- `modules/nf-core/dotmatch/panel_check/main.nf`
- `modules/nf-core/dotmatch/panel_check/meta.yml`
- `modules/nf-core/dotmatch/panel_check/tests/main.nf.test`
- `modules/nf-core/dotmatch/crispr_count/main.nf`
- `modules/nf-core/dotmatch/crispr_count/meta.yml`
- `modules/nf-core/dotmatch/crispr_count/tests/main.nf.test`
- `modules/nf-core/dotmatch/assay_run/main.nf`
- `modules/nf-core/dotmatch/assay_run/meta.yml`
- `modules/nf-core/dotmatch/assay_run/tests/main.nf.test`

## Container Pinning

The module files currently use the DotMatch 0.1.9 container tag template:

- Docker: `biocontainers/dotmatch:0.1.9--<bioconda_build>`
- Singularity: `https://depot.galaxyproject.org/singularity/dotmatch:0.1.9--<bioconda_build>`

Before opening an nf-core/modules PR, replace `<bioconda_build>` with the exact
build suffix published by the accepted Bioconda 0.1.9 package and mirrored by
BioContainers. Do not submit the placeholder.

## Fixtures

Each module includes a small `tests/data/` directory so the nf-test candidates
can run in an nf-core/modules checkout without a sibling DotMatch checkout.
The CRISPR and AssaySpec fixtures are copied from
`examples/workflows/fixtures/` and keep the same expected count and sample-QC
contract.

## Pre-PR Checklist

1. Clone or update a fork of `nf-core/modules`.
2. Copy `modules/nf-core/dotmatch/*` from this directory into
   `modules/nf-core/dotmatch/` in that fork.
3. Replace the 0.1.9 container tag templates with exact BioContainers tags.
4. Run the local DotMatch gates:

   ```bash
   make workflow-examples-ready
   make workflow-integration-test
   make reviewer-readiness-ready
   ```

5. In the nf-core/modules checkout, run nf-core lint and nf-test for each
   submitted module.
6. Open the PR with the module scope, output contract, fixture coverage, and
   container tag recorded in the description.

Record `docs/workflow-adoption.json` only after the external PR is merged or an
official module page is public. Use the merged PR or module page as
`adoption_url`, and use upstream CI, lint, or module documentation as
`evidence_url`.
