# Nextflow CRISPR Counting Example

This example wraps `dotmatch crispr-count` as a small Nextflow DSL2 workflow.
It is intended as an integration pattern for labs that already run Nextflow, not
as a separate benchmark claim.

From the repository root, create the small public CRISPR fixture:

```bash
python3 scripts/fetch_mageck_demo.py --small --out examples/crispr_guides/data
```

Run the workflow:

```bash
nextflow run examples/workflows/nextflow/main.nf \
  -c examples/workflows/nextflow/nextflow.config
```

Outputs are published under `examples/workflows/nextflow/output/`:

- `counts.mageck.tsv`: MAGeCK-compatible count matrix;
- `summary.json`: DotMatch assignment and QC summary;
- `sample_qc.tsv`: sample-level QC table that can be consumed by the
  MultiQC custom-content example.

The default config uses the same Yusa/MAGeCK fixture paths as
`examples/crispr_guides/README.md`: `ERR376998.fastq.gz`,
`ERR376999.fastq.gz`, `yusa_library.csv`, `guide_start=23`, and
`guide_length=19`.
