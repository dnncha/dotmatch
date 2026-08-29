# Workflow Submission Pack

DotMatch keeps runnable local workflow examples in `examples/workflows/`. This
page is the external handoff checklist for turning those examples into reviewed
community integrations. Do not mark `docs/workflow-adoption.json` as ready until
one of these submissions is accepted or released outside this repository.

## Current Submission Priority

1. **nf-core modules**: highest leverage after PyPI/Bioconda because accepted
   modules become reusable in nf-core pipelines and inherit nf-core update and
   container automation.
2. **MultiQC module**: makes DotMatch outputs recognizable in existing pipeline
   reports without per-project custom content.
3. **Galaxy / IUC wrappers**: reaches core facilities and wet-lab teams through
   ToolShed-reviewed tools.
4. **Snakemake handoff**: lower review overhead, strong for lab
   workflow adoption and QC report visibility.

## Shared Preflight

Before submitting any external integration:

```bash
make workflow-examples-ready
make workflow-integration-test
make reviewer-readiness-ready
make release-ready
```

The integration test runs:

- nf-test for `dotmatch/count`, `dotmatch/demux`, `dotmatch/audit`,
  `dotmatch/panel_check`, `dotmatch/crispr_count`, and `dotmatch/assay_run`;
- the minimal Nextflow pipeline in `examples/workflows/nf-core/pipeline`;
- the Snakemake workflow in `examples/workflows/snakemake`;
- Planemo lint for the Galaxy CRISPR count, demux, and panel-check wrappers;
- MultiQC custom-content and registered plugin smoke tests.

## nf-core Modules

Source payload:

- `examples/workflows/nf-core/upstream/modules/nf-core/dotmatch/count/`
- `examples/workflows/nf-core/upstream/modules/nf-core/dotmatch/demux/`
- `examples/workflows/nf-core/upstream/modules/nf-core/dotmatch/audit/`
- `examples/workflows/nf-core/upstream/modules/nf-core/dotmatch/panel_check/`
- `examples/workflows/nf-core/upstream/modules/nf-core/dotmatch/crispr_count/`
- `examples/workflows/nf-core/upstream/modules/nf-core/dotmatch/assay_run/`

External target:

- repository: `https://github.com/nf-core/modules`
- target path: `modules/nf-core/dotmatch/<subtool>/`

Reviewer notes:

- The upstream payload is prepared for DotMatch 0.1.9. Replace
  `0.1.9--<bioconda_build>` container templates with exact BioContainers tags
  before opening an nf-core/modules PR.
- Pin the container to a public DotMatch Bioconda/BioContainers release only
  after `make distribution-channels` verifies that release.
- Preserve DotMatch's `unique`, `ambiguous`, `none`, and `invalid` assignment
  semantics in module docs.
- Keep `task.ext.args` available for command-specific options.
- Use the tiny fixtures already included in each upstream module test directory;
  they are copied into the upstream tree so nf-test can run without a sibling
  DotMatch checkout.

Adoption record:

- Add a `nf_core_module` entry to `docs/workflow-adoption.json` only after the
  nf-core PR is merged or released.
- Use the merged PR or nf-core module page as `adoption_url`, and CI/lint
  evidence as `evidence_url`.

## Galaxy / IUC

Source payload:

- `examples/workflows/galaxy/dotmatch_crispr_count.xml`
- `examples/workflows/galaxy/dotmatch_demux.xml`
- `examples/workflows/galaxy/dotmatch_panel_check.xml`
- `examples/workflows/galaxy/test-data/`

External target:

- repository: `https://github.com/galaxyproject/tools-iuc`
- suggested path: `tools/dotmatch/`

Reviewer notes:

- Use the Bioconda package requirement once the intended DotMatch version is
  visible on Anaconda.
- Keep wrapper outputs plain TSV/JSON/FASTQ/HTML so Galaxy histories expose both
  the human report and workflow-readable tables.
- Run Planemo lint and tests in the IUC checkout before opening the PR.

Adoption record:

- Add a `galaxy_toolshed` entry only after IUC acceptance or ToolShed
  publication.

## Snakemake

Source payload:

- `examples/workflows/snakemake/Snakefile`
- `examples/workflows/snakemake/config.json`
- `examples/workflows/fixtures/`

External targets:

- Snakemake wrapper repository if converted to a reusable wrapper.
- Lab pipeline templates if maintained as a complete workflow example.

Reviewer notes:

- Keep `metric`, `k`, `ambiguity_policy`, and `ambiguous` explicit in config.
- Use `assay_k` separately from standalone CRISPR counting when a production
  assay needs stricter preflight semantics than the raw count rule.
- Keep absolute paths in generated AssaySpec TOML files.

Adoption record:

- Add `snakemake_workflow` or `independent_workflow` only for a public external
  repository or accepted wrapper.

## MultiQC

Source payload:

- custom content: `examples/workflows/multiqc/multiqc_config.yaml`
- registered plugin: `python/dotmatch/multiqc.py`
- fixtures: `examples/workflows/multiqc/data/` and
  `examples/workflows/fixtures/assay_out/`

External target:

- MultiQC plugin packaging or upstream module discussion after public workflow
  maintainers confirm the report shape.

Reviewer notes:

- The plugin parses `summary.json`, `sample_qc.tsv`, top-unmatched tables, and
  panel-safety summaries.
- Keep sample assignment rate, ambiguity rate, unmatched rate, and panel-safety
  status visible in general stats or module sections.

Adoption record:

- Add `multiqc_plugin` after a plugin package or upstream MultiQC integration is
  published outside this repository.
