# MultiQC Example

DotMatch ships a native MultiQC module in `dotmatch.multiqc`. The current source
tree registers its search patterns before MultiQC indexes input files. Install
the source candidate to check plugin discovery without a custom config:

```bash
pip install ".[multiqc]"
multiqc examples/workflows/multiqc/data --module dotmatch \
  -o examples/workflows/multiqc/output
```

The public 0.2.2 wheel contains the module entry point but not the early search
registration hook. Direct discovery with MultiQC 1.35 therefore remains a
next-release fix. For 0.2.2, use the custom-content configuration below from a
matching source checkout.

The native parser handles `sample_qc.tsv`, `crispr_qc.summary.tsv`,
`assay_manifest.summary.tsv`, `summary.json`, `panel_summary.json`, and
`top_unmatched.tsv` artifacts.

This directory also keeps a custom-content config for workflow environments
that cannot install the plugin yet. It is an integration pattern for workflow
reports, not a benchmark result.

Run from the repository root:

```bash
multiqc examples/workflows/multiqc/data \
  -c examples/workflows/multiqc/multiqc_config.yaml \
  -o examples/workflows/multiqc/output
```

The example data files are small fixtures with the same public schemas
documented in `docs/schemas.md`. In a real workflow, point MultiQC at the
directory containing DotMatch `sample_qc.tsv`, `*.sample_qc.tsv`,
`crispr_qc.summary.tsv`, or `*assay_manifest.summary.tsv` files and keep the
same config.

The report will include a `DotMatch Sample QC` custom-content table with
assignment rate, exact/rescue rates, ambiguous/no-match rates, target coverage,
library sparsity, dominance, and candidate-verification totals. It will also
include a `DotMatch Assay Manifest` table that links the workflow run back to
the primary `assay_report.html` and `assay_manifest.json` artifacts.

For the native parser implementation and dependency-free parser helpers, see
`python/dotmatch/multiqc.py`. The custom-content approach above remains
available when the plugin cannot be installed.
