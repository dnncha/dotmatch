# Snakemake CRISPR Counting Example

This example includes two small Snakemake rules. `dotmatch_crispr_count` runs
the native CRISPR counting command. `dotmatch_assay_run` runs an AssaySpec and
emits the assay report and manifest summary. It is intended as an integration
pattern for labs that already run Snakemake, not as benchmark evidence.

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

Rules that execute DotMatch declare
`examples/workflows/snakemake/envs/dotmatch.yaml`, which pins
`dotmatch=0.1.2` from Bioconda. Add `--use-conda` on supported Bioconda
platforms (`linux-64` or `osx-64`) when you want Snakemake to create the rule
environment; otherwise install DotMatch on `PATH` before running the workflow.

Outputs are written under `examples/workflows/snakemake/output/`:

- `samples.tsv`: sample sheet generated from `config.json`;
- `counts.mageck.tsv`: MAGeCK-compatible count matrix;
- `summary.json`: DotMatch assignment and QC summary;
- `sample_qc.tsv`: sample-level QC table that can be consumed by the
  MultiQC custom-content example;
- `assay/crispr_qc.html`, `assay/crispr_qc.json`,
  `assay/crispr_qc.summary.tsv`: CRISPR guide-count QC report, structured
  report, and workflow summary;
- `assay/assay_report.html`: primary human-readable AssaySpec report;
- `assay/assay_manifest.json`: full run provenance and command manifest;
- `assay/assay_manifest.summary.tsv`: manifest summary for MultiQC custom
  content.

The default config uses the same Yusa/MAGeCK fixture paths as
`examples/crispr_guides/README.md`: `ERR376998.fastq.gz`,
`ERR376999.fastq.gz`, `yusa_library.csv`, `guide_start=23`, and
`guide_length=19`.
