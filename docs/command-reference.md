# Command reference

This page is a compact map of the installed `dotmatch` command surface. Use
`dotmatch <command> --help` for option defaults and full argument details.

## Core Utilities

```bash
dotmatch --version
dotmatch citation
dotmatch dist SEQ1 SEQ2
dotmatch leq K SEQ1 SEQ2
```

`dist` prints the global edit distance between two short DNA strings. `leq`
prints `true` when the edit distance is less than or equal to `K`.

## Fixed-Window Assignment

```bash
dotmatch count --targets targets.tsv --reads sample.fastq.gz \
  --sample-label sample --target-start 23 --target-length 20 \
  --k 1 --metric hamming --out counts.tsv

dotmatch demux --barcodes barcodes.tsv --reads pooled.fastq.gz \
  --barcode-start 0 --barcode-length 8 --k 1 --out-dir demuxed/
```

Use `count` for guide, barcode, primer, adapter, feature-tag, and panel-target
counting from one fixed read window. Use `demux` when uniquely assigned inline
barcodes should be written to split FASTQ files.

Current correction-radius support:

- Hamming: `k=0..3` for fixed-length windows.
- Levenshtein: `k=0..2` for short insertion/deletion rescue.
- Metal GPU: experimental Darwin Metal path for eligible Hamming `k<=1`
  workloads only; CPU remains the assignment authority unless validation passes.

Before production Hamming `k=2` or `k=3`, run exact audit:

```bash
dotmatch audit --targets guides.tsv --k 3 --audit-mode exact --out-dir audit/
```

Proceed only when `safe_at_hamming_k2` or `safe_at_hamming_k3` is true for the
radius you plan to use.

## Pre-extracted Feature Matrices

```bash
dotmatch feature matrix \
  --observations feature_observations.tsv \
  --targets feature_library.tsv \
  --id-column observation_id \
  --cell-column cell_barcode \
  --sequence-column feature_seq \
  --metric hamming --k 1 --ambiguity-policy radius \
  --out-dir feature_matrix/
```

Use `feature matrix` when an upstream workflow has already made a headered
observation table with an explicit cell identifier and a feature sequence
window. It writes a deterministic sparse Matrix Market matrix with cells on
rows and features on columns, plus feature/cell axes, long-form counts,
per-observation assignments, per-cell QC, and a JSON summary.

Only unique assignments add matrix counts. The command does not pair FASTQ
reads, correct cell barcodes, deduplicate UMIs, or call cells; retain those
upstream decisions and provenance with the input table.

## AssaySpec Workflows

```bash
dotmatch assay new crispr --library guides.csv --reads-dir fastqs/ --out screen/
dotmatch assay start screen/assay.toml
dotmatch assay check assay.toml
dotmatch assay optimize assay.toml
dotmatch assay plan assay.toml
dotmatch assay run assay.toml
dotmatch assay autopsy assay.toml --out-dir autopsy/
```

`start` is the production entrypoint: it checks the spec, runs assignment, and
prints the reliability verdict. `optimize` writes an advisory CPU/GPU backend
recommendation; it does not change the count authority.

Supported templates include `crispr`, `feature-barcode`,
`inline-barcode-count`, `inline-barcode-demux`, `amplicon-panel`,
`oligo-adapter`, and `pair-count`.

## Barcode Workflows

```bash
dotmatch barcode infer --barcodes barcodes.tsv --reads pooled.fastq.gz --out offset_scan.tsv
dotmatch barcode audit --barcodes barcodes.tsv --k 1 --out-dir barcode_audit/
dotmatch barcode demux --barcodes barcodes.tsv --reads pooled.fastq.gz --out-dir demuxed/
dotmatch barcode count --barcodes barcodes.tsv --reads pooled.fastq.gz --out counts.tsv
dotmatch barcode autopsy --barcodes barcodes.tsv --reads pooled.fastq.gz --out-dir autopsy/
dotmatch barcode report --out-dir autopsy/
```

Use `barcode autopsy` when a barcode run has high unmatched, ambiguous, invalid,
or low-quality rescue rates. It combines offset inference, barcode audit,
demultiplexing/counting diagnostics, unmatched-window inspection, and a report.

## Barcode Panel Design

```bash
dotmatch panel check barcodes.tsv --k 1 --metric hamming --out-dir panel_check/
dotmatch panel design --n 96 --length 16 --preset illumina-inline-96 --out-dir panel/
dotmatch panel optimize vendor_barcodes.tsv --n 24 --out-dir optimized/
dotmatch panel simulate barcodes.tsv --reads 1000000 --out-dir simulation/
dotmatch panel layout barcodes.tsv --plate 96 --out plate_layout.tsv
dotmatch panel export barcodes.tsv --format illumina-samplesheet --out-dir sheets/
dotmatch panel compare old.tsv new.tsv --out-dir panel_compare/
dotmatch panel design-dual --samples 384 --i7-count 384 --i5-count 384 \
  --i7-length 10 --i5-length 10 --unique-dual --out-dir dual_panel/
```

Panel certificates use the same `unique`, `ambiguous`, `none`, and `invalid`
assignment semantics as counting and demultiplexing.

## CRISPR Workflows

```bash
dotmatch crispr new --library guides.csv --reads-dir fastqs/ --out screen/
dotmatch crispr start screen/assay.toml
dotmatch crispr qc --counts counts.mageck.tsv --sample-qc sample_qc.tsv --out crispr_qc.json
dotmatch crispr-count --library guides.csv --samples samples.tsv \
  --guide-start 23 --guide-length 20 --k 1 --out counts.mageck.tsv
```

The `crispr` namespace wraps AssaySpec helpers for guide-count workflows.
`crispr-count` is the direct single-command MAGeCK-compatible count writer.

## Compatibility Entrypoints

```bash
dotmatch guide-counter count --input sample.fastq.gz --samples sample \
  --library guides.tsv --output guide_counts
dotmatch guide-counter-count ...
dotmatch guide-count ...
```

These GuideCounter-compatible entrypoints preserve familiar input/output shapes
while delegating assignment to DotMatch's deterministic CPU counting engine.

## Diagnostics

```bash
dotmatch inspect-unmatched --targets targets.tsv --reads sample.fastq.gz \
  --target-start 23 --target-length 20 --k 1 --top 50 --out top_unmatched.tsv

dotmatch validate --targets targets.tsv --reads sample.fastq.gz \
  --target-start 23 --target-length 20 --k 1 --oracle scan
```

Use `inspect-unmatched` to diagnose frequent unassigned windows. Use `validate`
to compare indexed assignment with an exhaustive scan or Edlib oracle when
checking correctness-sensitive workflows.
