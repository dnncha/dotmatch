# Snakemake CRISPR Counting Example

This example wraps `dotmatch crispr-count` as a small Snakemake workflow. It is
intended as an integration pattern for labs that already run Snakemake, not as a
separate benchmark claim.

From the repository root, create the small public CRISPR fixture:

```bash
python3 scripts/fetch_mageck_demo.py --small --out examples/crispr_guides/data
```

Run the workflow:

```bash
snakemake \
  -s examples/workflows/snakemake/Snakefile \
  --configfile examples/workflows/snakemake/config.json \
  --cores 1
```

Outputs are written under `examples/workflows/snakemake/output/`:

- `samples.tsv`: sample sheet generated from `config.json`;
- `counts.mageck.tsv`: MAGeCK-compatible count matrix;
- `summary.json`: DotMatch assignment and QC summary;
- `sample_qc.tsv`: sample-level QC table that can be consumed by the
  MultiQC custom-content example.

The default config uses the same Yusa/MAGeCK fixture paths as
`examples/crispr_guides/README.md`: `ERR376998.fastq.gz`,
`ERR376999.fastq.gz`, `yusa_library.csv`, `guide_start=23`, and
`guide_length=19`.
