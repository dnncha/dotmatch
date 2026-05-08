# MultiQC Custom-Content Example

This example shows how to include DotMatch `sample_qc.tsv` output in a MultiQC
report using MultiQC custom content. It is an integration pattern for workflow
reports, not a benchmark claim and not a replacement for a future core MultiQC
module.

Run from the repository root:

```bash
multiqc examples/workflows/multiqc/data \
  -c examples/workflows/multiqc/multiqc_config.yaml \
  -o examples/workflows/multiqc/output
```

The example data file is a small `sample_qc.tsv` fixture with the same public
schema documented in `docs/schemas.md`. In a real workflow, point MultiQC at the
directory containing DotMatch `sample_qc.tsv` or `*.sample_qc.tsv` files and keep
the same config.

The report will include a `DotMatch Sample QC` custom-content table with
assignment rate, exact/rescue rates, ambiguous/no-match rates, target coverage,
library sparsity, dominance, and candidate-verification totals.
