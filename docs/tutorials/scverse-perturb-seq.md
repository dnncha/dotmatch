# DotMatch and scverse for Perturb-seq and Feature Barcodes

This tutorial covers the handoff from DotMatch feature assignments to
AnnData/scverse objects. Use `dotmatch feature matrix` when another workflow
has already produced one row per observation with an explicit cell identifier
and an extracted feature sequence.

## 1. Build a cell-by-feature matrix

The input table must be headered TSV (or CSV) and include a cell identifier and
the feature sequence window. A minimal table looks like this:

```text
observation_id	cell_barcode	feature_seq
read_001	AAACCTGAGAAACCAT	ACGTACGTACGTACGTACG
read_002	AAACCTGAGAAACCAT	ACGTACGTACGTACGTACG
read_003	AAACCTGAGCTAACAA	TGCATGCATGCATGCATGC
```

Run the matrix command with the input column names explicitly:

```bash
dotmatch feature matrix \
  --observations feature_observations.tsv \
  --targets feature_library.tsv \
  --id-column observation_id \
  --cell-column cell_barcode \
  --sequence-column feature_seq \
  --metric hamming \
  --k 1 \
  --ambiguity-policy radius \
  --out-dir feature_matrix/
```

`feature_matrix/matrix.mtx` is a sparse **cells × features** Matrix Market
matrix. Its row order is recorded in `barcodes.tsv`; its column order and
target sequences are recorded in `features.tsv`. `cell_feature_counts.tsv` is
the same unique-assignment result in long TSV form. `assignments.tsv`,
`cell_qc.tsv`, and `summary.json` retain the full outcome and run settings.

Only `unique` assignments add counts. `ambiguous`, `none`, and `invalid`
observations remain in the diagnostic artifacts instead of being forced into a
feature.

This command does not extract reads from FASTQ, pair read sides, correct cell
barcodes, deduplicate UMIs, or call cells. Perform those steps in the upstream
workflow before writing the observation table, and retain their provenance next
to the DotMatch output directory.

## 2. Load the matrix into AnnData

```python
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy.io import mmread

run = Path("feature_matrix")
cells = pd.read_csv(run / "barcodes.tsv", sep="\t")
features = pd.read_csv(run / "features.tsv", sep="\t")

feature_adata = ad.AnnData(
    X=mmread(run / "matrix.mtx").tocsr(),
    obs=cells.set_index("cell_barcode"),
    var=features.set_index("target_id"),
)
feature_adata.uns["dotmatch"] = {
    "summary": str(run / "summary.json"),
    "assignments": str(run / "assignments.tsv"),
    "cell_qc": str(run / "cell_qc.tsv"),
    "ambiguity_policy": "radius",
    "ambiguous_observations_counted": False,
}
```

## 3. Attach assignment rows when needed

For a smaller assignment table, `assignments_to_anndata` can build the same
kind of count matrix. Keep the `status` column and request unique-only counts:

```python
import pandas as pd
import dotmatch

assignments = pd.read_csv("feature_matrix/assignments.tsv", sep="\t")
feature_adata = dotmatch.assignments_to_anndata(
    assignments,
    cell_col="cell_barcode",
    feature_col="target_id",
    status_col="status",
    count_unique_only=True,
)
```

## 4. Use scanpy-style helpers for notebook-scale work

```python
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
    cell_col="cell_barcode",
    library=library,
    k=1,
    metric="hamming",
)
```

Use the command-line matrix writer for reproducible table-to-matrix runs. Use
the `dotmatch.tl` helpers for notebook-scale inspection where the observations
are already in AnnData.

## 5. Review per-cell assignment QC

For Perturb-seq and feature-barcode analysis, review `cell_qc.tsv` alongside
standard scRNA-seq QC. A high ambiguous or unmatched rate can indicate an
incorrect feature window, target library, orientation, or correction radius.
The feature matrix alone does not establish cell identity or UMI-collapsed
molecule counts.
