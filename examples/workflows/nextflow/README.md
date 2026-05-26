# Nextflow CRISPR Counting Example

This example includes two small Nextflow DSL2 processes. `DOTMATCH_CRISPR_COUNT`
runs the native CRISPR counting command. `DOTMATCH_ASSAY_RUN` runs an AssaySpec
and emits the assay report and manifest summary. It is intended as an integration
pattern for labs that already run Nextflow, not as benchmark evidence.

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
  MultiQC custom-content example;
- `crispr_qc.html`, `crispr_qc.json`, `crispr_qc.summary.tsv`: CRISPR
  guide-count QC report, structured report, and workflow summary;
- `assay_report.html`: primary human-readable AssaySpec report;
- `assay_manifest.json`: full run provenance and command manifest;
- `assay_manifest.summary.tsv`: manifest summary for MultiQC custom content.

The default config uses the same Yusa/MAGeCK fixture paths as
`examples/crispr_guides/README.md`: `ERR376998.fastq.gz`,
`ERR376999.fastq.gz`, `yusa_library.csv`, `guide_start=23`, and
`guide_length=19`.

`nextflow.config` keeps `metric`, `k`, `ambiguity_policy`, and `ambiguous`
explicit so production pipelines do not inherit changed defaults silently. For
fixed-length guides where you want larger mismatch radii, run an exact library
audit first and then switch the params to `metric = 'hamming'` with `k = 2` or
`k = 3` only when the audit reports the corresponding Hamming radius as safe.
