"""Explicit dataframe and sparse AnnData adapters; no optional imports at startup."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

from ._validation import identifier, named_targets, sequence
from .count_io import read_count_table
from .target_io import ID_COLUMNS, SEQ_COLUMNS


def _values(column) -> list:
    if hasattr(column, "to_list"):
        return column.to_list()
    if hasattr(column, "tolist"):
        return column.tolist()
    return list(column)


def _column(columns, explicit, aliases, *, fallback=None, name="column"):
    if explicit is not None:
        if explicit in columns:
            return explicit
        if isinstance(explicit, int) and not isinstance(explicit, bool) and 0 <= explicit < len(columns):
            return columns[explicit]
        raise ValueError(f"{name} {explicit!r} is not present")
    hits = [col for col in columns if str(col).lower() in aliases]
    if len(hits) > 1:
        raise ValueError(f"multiple possible {name}s {hits!r}; specify one explicitly")
    return hits[0] if hits else fallback


def targets_from_dataframe(df: Any, id_col=None, seq_col=None) -> list[tuple[str, str]]:
    """Read named or explicitly selected columns from pandas or Polars.

    A one-column sequence table receives target_0, target_1, ... IDs. With no
    recognized names, two-column tables retain the positional ID/sequence form.
    Missing values and duplicate IDs are rejected before matching. No PyArrow
    conversion is required for Polars inputs.
    """
    if not hasattr(df, "columns"):
        raise TypeError("targets_from_dataframe requires a dataframe with columns")
    columns = list(df.columns)
    if not columns or len(set(columns)) != len(columns):
        raise ValueError("target dataframe columns must be nonempty and unique")
    # A canonical sequence column takes precedence over legacy names such as
    # "guide", which is also a common target-ID column name.
    aliases = SEQ_COLUMNS if any(str(col).lower() in SEQ_COLUMNS for col in columns) else ("dna", "target", "guide", "barcode")
    seq = _column(columns, seq_col, aliases,
                  fallback=columns[1] if len(columns) > 1 else columns[0], name="sequence column")
    id_candidates = [col for col in columns if col != seq]
    target_id = _column(columns if id_col is not None else id_candidates, id_col, ID_COLUMNS,
                        fallback=id_candidates[0] if id_candidates else None, name="ID column")
    if target_id == seq:
        raise ValueError("ID and sequence columns must differ")
    seqs = _values(df[seq])
    ids = _values(df[target_id]) if target_id is not None else [f"target_{i}" for i in range(len(seqs))]
    if len(ids) != len(seqs):
        raise ValueError("target ID and sequence columns have different lengths")
    return named_targets(zip(ids, seqs))


def _read_sequences(reads, seq_col=None, id_col=None):
    if isinstance(reads, (str, bytes, Path)):
        raise TypeError("reads must be a sequence collection or dataframe; use stream_assign for a FASTQ path")
    if hasattr(reads, "columns"):
        columns = list(reads.columns)
        if not columns or len(set(columns)) != len(columns):
            raise ValueError("read dataframe columns must be nonempty and unique")
        seq = _column(columns, seq_col, (*SEQ_COLUMNS, "read_seq", "read_sequence", "dna"),
                      fallback=columns[0] if len(columns) == 1 else None, name="read sequence column")
        if seq is None:
            raise ValueError("no read sequence column found; specify read_seq_col")
        rid = _column(columns if id_col is not None else [col for col in columns if col != seq], id_col, ("read_id", "id", "name"), name="read ID column")
        if rid == seq:
            raise ValueError("read ID and sequence columns must differ")
        ids = _values(reads[rid]) if rid is not None else (_values(reads.index) if hasattr(reads, "index") else None)
        values = _values(reads[seq])
    else:
        if seq_col is not None or id_col is not None:
            raise ValueError("read column options require a dataframe")
        values = _values(reads)
        ids = _values(reads.index) if hasattr(reads, "index") and not callable(reads.index) else None
    return [sequence(value, f"read sequence {i}") for i, value in enumerate(values)], ids


def assign_dataframe(reads, targets, k=1, policy="radius", metric="levenshtein",
                     read_ids=None, target_names=None, *, read_seq_col=None,
                     read_id_col=None, target_seq_col=None, target_id_col=None):
    """Assign lists, Series, or dataframes and return a labelled pandas table.

    High-level sequence inputs normalize ASCII case consistently. Use the low-
    level Matcher/assign functions for literal, case-sensitive byte matching.
    Explicit label arrays must match input lengths exactly.
    """
    from . import core

    core._ensure_pandas()
    if hasattr(targets, "columns"):
        named = targets_from_dataframe(targets, target_id_col, target_seq_col)
        names, sequences = map(list, zip(*named))
    else:
        if target_seq_col is not None or target_id_col is not None:
            raise ValueError("target column options require a dataframe")
        named = core._normalize_targets(targets if isinstance(targets, (str, Path)) else _values(targets))
        names, sequences = map(list, zip(*named))
    if target_names is not None:
        names = [identifier(value, "target name") for value in target_names]
        if len(names) != len(sequences) or len(set(names)) != len(names):
            raise ValueError("target_names must contain one distinct name per target")
    read_sequences, inferred_ids = _read_sequences(reads, read_seq_col, read_id_col)
    ids = read_ids if read_ids is not None else inferred_ids
    if ids is not None:
        ids = [identifier(value, "read ID") for value in ids]
        if len(ids) != len(read_sequences):
            raise ValueError("read_ids must contain one ID per read")
    if metric == "levenshtein":
        results = core.assign(read_sequences, sequences, k=k, policy=policy)
    elif metric == "hamming":
        results = core.assign_hamming(read_sequences, sequences, k=k, policy=policy)
    elif metric == "exact":
        from ._validation import integer
        if integer(k, "k") != 0:
            raise ValueError("metric='exact' requires k=0")
        results = core.assign_exact(read_sequences, sequences, policy=policy)
    else:
        raise ValueError("metric must be 'levenshtein', 'hamming', or 'exact'")
    return core.results_to_dataframe(results, target_names=names, read_ids=ids)


def counts_tsv_to_anndata(counts_path, *, sample_cols=None,
                         var_cols=("target_id", "target_seq", "gene")):
    """Load lossless raw integer counts as a samples-by-targets sparse CSR matrix.

    Missing counts, duplicate axes, fractional counts and int64 overflow fail
    explicitly. No permissive fallback parser is used. For detailed DotMatch
    tables, only *_count_total columns become samples by default.
    """
    from . import core
    core._ensure_anndata()
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse

    table = read_count_table(counts_path, sample_cols=sample_cols, var_cols=var_cols)
    rows, cols, values = [], [], []
    for target, counts in enumerate(table.counts):
        for sample, value in enumerate(counts):
            if value > np.iinfo(np.int64).max:
                raise OverflowError(f"count for {table.target_ids[target]}/{table.sample_names[sample]} exceeds int64; no lossy conversion was made")
            if value:
                rows.append(sample)
                cols.append(target)
                values.append(value)
    matrix = sparse.csr_matrix((np.asarray(values, dtype=np.int64), (rows, cols)),
                               shape=(len(table.sample_names), len(table.target_ids)), dtype=np.int64)
    obs = pd.DataFrame(index=pd.Index(table.sample_names, dtype=str))
    var = pd.DataFrame({key: list(values) for key, values in table.metadata.items()},
                       index=pd.Index(table.target_ids, dtype=str))
    result = ad.AnnData(X=matrix, obs=obs, var=var)
    result.layers["counts"] = matrix.copy()
    result.uns["dotmatch"] = {"source": "counts_tsv_to_anndata", "raw_integer_counts": True,
                              "sample_columns": list(table.sample_columns)}
    return result


def _axis(values, name):
    result = [identifier(value, name) for value in values]
    if len(set(result)) != len(result):
        raise ValueError(f"{name} values must be unique")
    return result


def _read_assignment_table(path):
    """Reject ragged/duplicate-column TSV before pandas can infer a new index."""
    import csv
    import gzip
    import pandas as pd

    source = Path(path)
    opener = gzip.open if source.name.lower().endswith(".gz") else open
    with opener(source, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", strict=True)
        header = next(reader, None)
        if not header or any(not name or name != name.strip() or any(ord(c) < 32 or ord(c) == 127 for c in name) for name in header):
            raise ValueError("assignment columns must be nonempty without surrounding whitespace or controls")
        if len(set(header)) != len(header):
            raise ValueError("assignment columns must be unique")
        rows = []
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"assignment row {reader.line_num} has {len(row)} fields; expected {len(header)}")
            rows.append(row)
    return pd.DataFrame(rows, columns=header, dtype=object)


def assignments_to_anndata(assignments, *, cell_col="cell_barcode", feature_col="target_name",
                           status_col=None, count_unique_only=True, include_ambiguous_per_cell=False,
                           cell_names: Sequence[str] | None = None, feature_names: Sequence[str] | None = None):
    """Build sparse cell-by-feature counts from explicitly cell-labelled rows.

    This counts observations, not molecules: it does not infer cells from read
    IDs, correct cell barcodes, deduplicate UMIs, or call perturbations. All
    observed cells remain on the axis, including cells with no unique calls.
    Optional cell_names/feature_names fix axis order and retain zero-count axes.
    """
    from . import core
    core._ensure_anndata()
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse

    if count_unique_only is not True:
        raise ValueError("only unique assignments may contribute counts; count_unique_only must be True")
    if isinstance(assignments, (str, Path)):
        # Preserve NA-like and numeric-looking identifiers exactly.
        frame = _read_assignment_table(assignments)
    elif hasattr(assignments, "columns"):
        frame = assignments
    else:
        raise TypeError("assignments must be a dataframe or a TSV path")
    columns = list(frame.columns)
    if len(set(columns)) != len(columns):
        raise ValueError("assignment columns must be unique")
    if cell_col not in columns:
        raise ValueError(f"missing explicit cell column {cell_col!r}; cell barcodes are never inferred from read IDs")
    if feature_col not in columns:
        if feature_col == "target_name" and "target_id" in columns:
            feature_col = "target_id"
        else:
            raise ValueError(f"missing feature column {feature_col!r}")
    if status_col is None:
        status_col = "status_name" if "status_name" in columns else "status"
    if status_col not in columns:
        raise ValueError(f"missing assignment status column {status_col!r}")
    cells = [identifier(value, "cell barcode") for value in _values(frame[cell_col])]
    features = _values(frame[feature_col])
    aliases = {"unique": "unique", "1": "unique", "1.0": "unique", "ambiguous": "ambiguous",
               "2": "ambiguous", "2.0": "ambiguous", "none": "none", "0": "none", "0.0": "none",
               "invalid": "invalid", "-1": "invalid", "-1.0": "invalid"}
    statuses = []
    for value in _values(frame[status_col]):
        key = str(value).strip().lower()
        if key not in aliases:
            raise ValueError(f"unknown assignment status {value!r}")
        statuses.append(aliases[key])
    unique_features = [identifier(features[i], "assigned feature") for i, status in enumerate(statuses) if status == "unique"]
    cell_axis = list(dict.fromkeys(cells)) if cell_names is None else _axis(cell_names, "cell name")
    feature_axis = list(dict.fromkeys(unique_features)) if feature_names is None else _axis(feature_names, "feature name")
    cell_index = {name: i for i, name in enumerate(cell_axis)}
    feature_index = {name: i for i, name in enumerate(feature_axis)}
    if any(name not in cell_index for name in cells):
        raise ValueError("cell_names omits an observed cell; refusing to drop observations")
    if any(name not in feature_index for name in unique_features):
        raise ValueError("feature_names omits a uniquely assigned feature; refusing to drop counts")
    rows, cols = [], []
    qc = {status: np.zeros(len(cell_axis), dtype=np.int64) for status in ("unique", "ambiguous", "none", "invalid")}
    feature_iterator = iter(unique_features)
    for cell, status in zip(cells, statuses):
        row = cell_index[cell]
        qc[status][row] += 1
        if status == "unique":
            rows.append(row)
            cols.append(feature_index[next(feature_iterator)])
    matrix = sparse.csr_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)),
                               shape=(len(cell_axis), len(feature_axis)), dtype=np.int64)
    obs = pd.DataFrame(index=pd.Index(cell_axis, dtype=str))
    obs["n_observations"] = sum(qc.values())
    if include_ambiguous_per_cell:
        for status, label in (("unique", "unique"), ("ambiguous", "ambiguous"), ("none", "unmatched"), ("invalid", "invalid")):
            obs[f"{label}_count"] = qc[status]
    result = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=pd.Index(feature_axis, dtype=str)))
    result.layers["counts"] = matrix.copy()
    result.uns["dotmatch"] = {"source": "assignments_to_anndata", "unique_only": True,
                              "unit": "observations", "cell_barcodes_inferred": False,
                              "umi_deduplicated": False}
    return result
