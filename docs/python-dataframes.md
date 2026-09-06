# Dataframes and sparse count matrices

These APIs connect known-target assignment to a researcher's existing tables.
The changes described here are **unreleased source changes after 0.5.0** until
listed under a subsequent release in the changelog.

Install the optional stack used by the workflow: `dotmatch[pandas]`,
`dotmatch[polars]` or `dotmatch[anndata]`. Basic native matching does not import
these scientific packages merely because they happen to be installed.

## Assign named dataframe columns

```python
import pandas as pd
import dotmatch

reads = pd.DataFrame({"read_id": ["0001", "0002"], "sequence": ["acgt", "TTTT"]})
targets = pd.DataFrame({"sequence": ["ACGT", "TTTT"], "target_id": ["g1", "g2"]})
assigned = dotmatch.assign_dataframe(reads, targets, k=0, metric="exact")
assert assigned["target_name"].tolist() == ["g1", "g2"]
```

Column order is not part of the experiment: recognized column names are used
before positional conventions. Unusual names can be specified with
`read_seq_col`, `read_id_col`, `target_seq_col` and `target_id_col`. A string
selects a label; an integer selects an existing integer label first, otherwise
its position in the original dataframe. Ambiguous named columns require an
explicit choice. Polars uses its own column access, without a PyArrow round trip.

High-level dataframe/FASTQ inputs uppercase sequences. Missing values are
errors, not the strings `NAN` or `NONE`. Duplicate sequences with different IDs
remain distinct candidates; duplicate IDs are rejected. Explicit `read_ids` and
`target_names` must cover the inputs exactly. For literal case-sensitive byte
matching, use the low-level `Matcher` or `assign` APIs.

## Import raw counts without rounding

```python
import dotmatch

adata = dotmatch.counts_tsv_to_anndata(
    "screen.counts.tsv",
    sample_cols=["treated", "control"],
)
assert adata.obs_names.tolist() == ["treated", "control"]
adata.write_h5ad("screen.counts.h5ad")
```

Rows in `adata.X` are samples; columns are guide/target IDs. The matrix and its
`counts` layer are sparse CSR with int64 counts. Gene and sequence metadata are
preserved in `var`. Text such as `NA` and `0001` stays text. Missing, negative,
non-finite and fractional counts fail explicitly. Counts above int64's maximum
are rejected rather than converted to an imprecise floating-point value.

For native detailed DotMatch tables, only `*_count_total` fields become samples
by default. Exact and corrected components must not become additional samples.
`sample_cols` names the source columns in the desired order; a detailed table's
`control_count_total` column becomes the sample name `control`. An explicitly
selected subset does not exempt other sample counts from input validation.

## Build a cell-by-feature matrix from labelled observations

```python
import dotmatch

adata = dotmatch.assignments_to_anndata(
    "assignments.tsv",
    cell_col="cell_barcode",
    feature_col="target_id",
    status_col="status",
    cell_names=["cell_A", "cell_B", "cell_C"],
    feature_names=["guide_1", "guide_2"],
    include_ambiguous_per_cell=True,
)
```

Cell labels must come from an explicit column, never guessed from a read ID.
Every observed cell remains on the axis, even when it has only ambiguous,
unmatched or invalid observations. Optional axes preserve supplied order and
include zero-count cells/features; omitting an observed cell or a uniquely
assigned feature is an error rather than silently dropping counts. Without
explicit axes, first-occurrence order is retained.

Only unique assignments contribute to `X`. `obs.n_observations` accounts for all
rows. The QC option adds `unique_count`, `ambiguous_count`, `unmatched_count` and
`invalid_count`. The input dataframe is not modified. An all-unassigned input
can therefore produce a valid cell-by-zero-features sparse matrix.

**The unit is observations, not necessarily molecules.** This helper does not
correct cell barcodes, deduplicate UMIs, call cells, subtract ambient guides or
infer a perturbation. Preprocess and validate those stages separately. Passing
`count_unique_only=False` is rejected because there is no implemented,
scientifically justified fractional or forced assignment policy in this adapter.

Sparse storage follows observed nonzero entries plus axis pointers, not the
full cells-times-features rectangle. A regression test constructs 50,000 cells
by 50,000 features with one nonzero, without allocating a dense matrix. This is
a software storage test, not a sequencing-throughput benchmark.

## Numeric and runtime boundaries

Thresholds must be actual nonnegative integers representable by the native C
API. Booleans, decimal values, strings and oversized integers are rejected, not
rounded or wrapped. Integer scalar types such as NumPy integers are supported.

A `Matcher` serializes native calls against `close()` on that same instance.
Use distinct instances for independent concurrent jobs. This protects index
lifetime; it is not a claim of parallel acceleration for one shared instance.

Native library/executable overrides are authoritative: a missing override
fails rather than falling back to another binary. Normal loading uses packaged
binaries or the recognized source checkout, never the working directory.

The experimental `assign_posterior` helper rejects invalid priors and preserves
ties. Its probabilities remain conditional on its simple likelihood model;
these fixes do not establish biological calibration.
