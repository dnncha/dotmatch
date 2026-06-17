# DotMatch to scverse for Perturb-seq and Feature Barcodes

This tutorial shows the intended handoff from DotMatch assignment artifacts to
AnnData/scverse objects. Use the CLI for FASTQ-scale assignment, then load the
small, stable TSV outputs into Python.

## 1. Count guide or feature-barcode reads

```bash
dotmatch count \
  --targets guides.tsv \
  --reads guide_capture_R2.fastq.gz \
  --sample-label guide_capture \
  --target-start 63 \
  --target-length 19 \
  --k 1 \
  --metric hamming \
  --ambiguity-policy radius \
  --ambiguous discard \
  --out guide_counts.tsv \
  --summary guide_summary.json \
  --sample-qc guide_sample_qc.tsv \
  --assignments guide_assignments.tsv
```

For TotalSeq/CITE-seq-style feature barcodes, use the feature-barcode table as
`--targets` and set `--target-start` / `--target-length` to the antibody or
feature barcode window.

## 2. Load counts into AnnData

```python
import dotmatch

guide_adata = dotmatch.counts_tsv_to_anndata("guide_counts.tsv")
guide_adata.uns["dotmatch_summary_json"] = "guide_summary.json"
guide_adata.uns["dotmatch_sample_qc_tsv"] = "guide_sample_qc.tsv"
```

The count matrix contains uniquely assigned targets only. Ambiguous reads are
reported in `summary.json`, `sample_qc.tsv`, and `assignments.tsv`; they are not
silently assigned to a guide or feature.

## 3. Attach per-read assignments when cell barcodes are available

If your assignment table includes a cell barcode column, convert it to an
AnnData observation-level table:

```python
import dotmatch

assign_adata = dotmatch.assignments_to_anndata(
    "guide_assignments.tsv",
    cell_col="cell_barcode",
    target_col="target_id",
)
```

For custom pipelines, join DotMatch assignments to cell barcodes before this
step. Keep `assignment_status` so downstream filtering can distinguish unique,
ambiguous, unmatched, and invalid reads.

## 4. Use scanpy-style helpers

```python
import scanpy as sc
import dotmatch.tl as dm_tl

library = [
    {"id": "guide_A", "sequence": "ACGTACGTACGTACGTACG"},
    {"id": "guide_B", "sequence": "TGCATGCATGCATGCATGC"},
]

dm_tl.assign_features(
    adata,
    seq_col="guide_sequence",
    library=library,
    k=1,
    metric="hamming",
)

feature_adata = dm_tl.feature_counts(
    adata,
    seq_col="guide_sequence",
    library=library,
    k=1,
    metric="hamming",
)
```

Use this path for notebook-scale inspection and prototypes. For production
FASTQ processing, prefer `dotmatch count` so the exact command, assignment
engine, ambiguity policy, and QC summaries are written as reproducible files.

## 5. Recommended scverse metadata

Store DotMatch provenance in `.uns`:

```python
adata.uns["dotmatch"] = {
    "summary": "guide_summary.json",
    "sample_qc": "guide_sample_qc.tsv",
    "assignments": "guide_assignments.tsv",
    "ambiguity_policy": "radius",
    "ambiguous_reads_counted": False,
}
```

For Perturb-seq analysis, keep guide assignment QC next to standard scRNA-seq
QC. A high ambiguous or unmatched rate usually means the guide window, barcode
library, or correction radius should be checked before interpreting guide-level
effects.
