# Workflow Submission Records

DotMatch keeps runnable examples in `examples/workflows/`. This page links the
corresponding public upstream records and names the remaining acceptance gates.
The records were checked on 2026-08-28. An open pull request is not an accepted
or released integration.

See `docs/ecosystem-status.md` for package-manager and registry status.

## Shared preflight

```bash
make workflow-examples-ready
make workflow-integration-test
make release-ready
```

The integration test covers the local nf-test modules and minimal Nextflow
pipeline, the Snakemake workflow, Galaxy wrapper lint, and both MultiQC custom
content and plugin discovery.

## nf-core modules

- Public record: [nf-core/modules #12156](https://github.com/nf-core/modules/pull/12156)
- State: open, not accepted or released
- Checked head: `77e849b86cae10557a8b17f9c86fa87f6833ece2`
- Scope: `modules/nf-core/dotmatch/crispr_count/`
- DotMatch release: 0.2.2
- Container: `quay.io/biocontainers/dotmatch:0.2.2--py311h13f8228_1`

Additional local module payloads remain under
`examples/workflows/nf-core/upstream/modules/nf-core/dotmatch/`. Keep separate
subtools as separately reviewed changes. Preserve `unique`, `ambiguous`,
`none`, and `invalid` semantics, expose `task.ext.args`, and record the exact
DotMatch version in `versions.yml`.

Acceptance gate: merge into `nf-core/modules`, followed by verification of the
published module record. Only then add an `nf_core_module` entry to
`docs/workflow-adoption.json`.

## Galaxy / IUC

- Public record: [galaxyproject/tools-iuc #8336](https://github.com/galaxyproject/tools-iuc/pull/8336)
- State: open, not accepted or published to a ToolShed
- Checked head: `b547330d0b3f9fb38368eac3a94fd84098d51031`
- Upstream checks: wrapper lint, containerized Planemo tests, and result
  aggregation pass at the checked head
- Scope: CRISPR count, demultiplexing, panel check, and assay-run wrappers

Source wrappers and fixtures are in `examples/workflows/galaxy/`. Keep outputs
as plain TSV, JSON, FASTQ, and HTML datasets so Galaxy histories expose both
human reports and workflow-readable files.

Acceptance gate: IUC merge, then ToolShed publication. These are separate
states; add `galaxy_toolshed` to `docs/workflow-adoption.json` only after the
public ToolShed record is verified.

## Snakemake wrappers

- Public record: [snakemake/snakemake-wrappers #5825](https://github.com/snakemake/snakemake-wrappers/pull/5825)
- State: open, not accepted or released
- Checked head: `a72d2bb8bdeb97d8ac506cbc249925f969888b6f`
- Local check: Black passes for `bio/dotmatch/crispr-count/wrapper.py`
- Upstream checks: the Code quality and Tests runs require upstream maintainer
  approval at the checked head

The repository workflow remains under `examples/workflows/snakemake/`. Keep
`metric`, `k`, `ambiguity_policy`, and ambiguous-output handling explicit.

Acceptance gate: an upstream maintainer must approve the workflow runs; after
checks pass, obtain review and merge, then verify the released wrapper record
before adding a
`snakemake_workflow` entry to `docs/workflow-adoption.json`.

## MultiQC

- Public record: [MultiQC/MultiQC #3629](https://github.com/MultiQC/MultiQC/pull/3629)
- State: open, not accepted or released
- Checked head: `166f94ce70f2bc1fdbc94f460a4c857511bf1416`
- Upstream checks: eight checks pass; one additional Python 3.9 check is
  recorded as cancelled

DotMatch also exposes its parser as a package plugin. The source tree now uses
a `before_config` hook so search patterns are present before MultiQC indexes
input files. The next-release source candidate discovers `sample_qc.tsv`,
`summary.json`, `panel_summary.json`, `crispr_qc.summary.tsv`,
`assay_manifest.summary.tsv`, and `top_unmatched.tsv` without a custom config.
This fix is not part of the public 0.2.2 package.

Acceptance gates are separate: merge and release of the upstream MultiQC
module, or a DotMatch release containing the independently packaged plugin fix.
Record only the state that has actually been published.

## bio.tools and WorkflowHub

`docs/registries/biotools.yml` is draft metadata. No accepted DotMatch record
was found in bio.tools on 2026-08-28. No exact DotMatch workflow record was
found in WorkflowHub. Neither surface is recorded as submitted or accepted.
