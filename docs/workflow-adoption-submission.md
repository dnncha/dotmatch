# Workflow Adoption Submission Dossier

This dossier turns the in-repository workflow examples into a concrete adoption
submission plan. It is not external adoption, not an upstream nf-core module,
not a Galaxy ToolShed release, and not a MultiQC plugin claim. Public adoption
is recorded only in `docs/workflow-adoption.json` after stable external links
exist.

## Required Checks

Run the local workflow gate before any upstream submission:

```bash
make workflow-examples-ready
```

Run the post-adoption gate only after a real external integration has landed
and `docs/workflow-adoption.json` has been updated:

```bash
make workflow-adoption-status
```

## Submission Targets

- nf-core/modules: propose `examples/workflows/nf-core/modules/local/dotmatch/crispr_count` as a `dotmatch/crispr_count` module candidate after a released DotMatch package or container is available.
- Galaxy ToolShed: package `examples/workflows/galaxy/dotmatch_crispr_count.xml` with Planemo tests and a released Bioconda package or container.
- MultiQC: keep the current custom-content recipe usable, or propose a small MultiQC plugin once multiple public DotMatch workflows emit `sample_qc.tsv`.
- Snakemake: publish the checked CRISPR counting workflow or add DotMatch to an existing CRISPR screen workflow repository.
- Nextflow: publish the checked DSL2 CRISPR counting workflow or add DotMatch to an existing screen-analysis pipeline.

## Pre-Submission Checklist

- Release DotMatch through at least one stable install channel referenced by the target repository.
- Keep all help text aligned with `docs/scientific-claims.md`; do not use genome-aligner, adapter-trimmer, full Perturb-seq, or production demultiplexing wording.
- Use `examples/workflows/fixtures/` for small FASTQ and guide-library fixtures that exercise unique, ambiguous, unmatched, and invalid outcomes.
- Keep expected MAGeCK counts and sample-QC outputs alongside the fixtures so wrapper tests have stable assertions.
- For nf-core, adapt `examples/workflows/nf-core/modules/local/dotmatch/crispr_count/tests/main.nf.test`, keep `versions.yml`, and lint against the target repository conventions.
- For Galaxy, keep the embedded Planemo test and `examples/workflows/galaxy/test-data/` fixtures current, then run `planemo lint` and `planemo test` before ToolShed publication.
- For MultiQC, verify the `sample_qc.tsv` columns against `docs/schemas.md`.
- For Snakemake and Nextflow, keep ambiguity policy, metric, guide start, guide length, and output paths explicit.

## Manifest Record Template

After an external integration lands, add one entry to `docs/workflow-adoption.json`:

```json
{
  "id": "nfcore_dotmatch_crispr_count",
  "type": "nf_core_module",
  "name": "nf-core dotmatch/crispr_count module",
  "status": "accepted",
  "adoption_url": "https://example.org/stable-adoption-page",
  "evidence_url": "https://example.org/review-or-release-record",
  "validation_notes": "Accepted upstream module with tests, version reporting, and scoped help text.",
  "recorded_date": "YYYY-MM-DD"
}
```

Use a unique `id`, set `status` to `accepted`, `released`, or `published`, and
replace every placeholder URL with stable non-example `https://` adoption and
evidence links. `recorded_date` must use `YYYY-MM-DD`.

Do not set workflow adoption status to `ready` until stable external HTTPS adoption and evidence URLs exist and `make workflow-adoption-status` passes.
