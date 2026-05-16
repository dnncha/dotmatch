# MultiQC Custom-Content Example

This example shows how to include DotMatch `sample_qc.tsv` and
`assay_manifest.summary.tsv` outputs in a MultiQC report using MultiQC custom
content. It is an integration pattern for workflow reports, not a benchmark
result and not a core MultiQC module.

Run from the repository root:

```bash
multiqc examples/workflows/multiqc/data \
  -c examples/workflows/multiqc/multiqc_config.yaml \
  -o examples/workflows/multiqc/output
```

The example data files are small fixtures with the same public schemas
documented in `docs/schemas.md`. In a real workflow, point MultiQC at the
directory containing DotMatch `sample_qc.tsv`, `*.sample_qc.tsv`, or
`*assay_manifest.summary.tsv` files and keep the same config.

The report will include a `DotMatch Sample QC` custom-content table with
assignment rate, exact/rescue rates, ambiguous/no-match rates, target coverage,
library sparsity, dominance, and candidate-verification totals. It will also
include a `DotMatch Assay Manifest` table that links the workflow run back to
the primary `assay_report.html` and `assay_manifest.json` artifacts.
