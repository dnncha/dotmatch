"""
dotmatch.tl - scverse / scanpy-style tools for seamless AnnData integration.

This provides a `tl` (tools) submodule following the conventions of scanpy, pertpy, etc.
Functions operate on AnnData in-place where appropriate, or return new AnnData objects
for count matrices. They build on the core assign / anndata helpers while enforcing
DotMatch's scientific guarantees (unique-only assignment by default, explicit ambiguity).

Typical usage in a notebook after having cell-barcoded reads or extracted windows:

    import scanpy as sc
    import dotmatch as dm

    # adata.obs has 'feature_seq' column (extracted fixed window from R2 etc.)
    dm.tl.assign_features(
        adata,
        library="feature_barcodes.tsv",
        seq_col="feature_seq",
        k=1,
        metric="hamming",
        result_col="assigned_feature",
    )

    # Or build a dedicated feature count AnnData
    feature_adata = dm.tl.feature_counts(
        adata,
        library=library_df,
        seq_col="feature_seq",
        cell_col="cell_barcode",
    )

All operations respect the radius ambiguity policy by default and only count
`unique` assignments.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any, Sequence

from ._optional import optional_module

ad, _HAS_ANNDATA = optional_module("anndata")
pd, _HAS_PANDAS = optional_module("pandas")

# Use absolute import from the installed package name to avoid relative import
# issues when running under PYTHONPATH or pytest collection.
from dotmatch.core import (
    Matcher,
    assignments_to_anndata,
    counts_tsv_to_anndata,
    targets_from_dataframe,
)


def _ensure_anndata() -> None:
    if not (_HAS_ANNDATA and _HAS_PANDAS):
        raise ImportError(
            "anndata (and pandas) are required for dotmatch.tl. "
            "Install with: pip install 'dotmatch[anndata]'"
        )


def _load_library(library: Any) -> list[tuple[str, str]]:
    """Accept a path, DataFrame, mappings, tuples, or sequences (auto ids)."""
    if isinstance(library, (str, Path)):
        # Assume tsv/csv with id,seq or just seqs
        try:
            df = pd.read_csv(library, sep="\t")
            if len(df.columns) >= 2:
                return targets_from_dataframe(df)
            else:
                seqs = df.iloc[:, 0].astype(str).tolist()
                return [(f"target_{i}", s) for i, s in enumerate(seqs)]
        except Exception:
            # fallback single column?
            pass
    if isinstance(library, pd.DataFrame):
        return targets_from_dataframe(library)
    if isinstance(library, (list, tuple)):
        if library and isinstance(library[0], (list, tuple)) and len(library[0]) == 2:
            return list(library)
        if library and isinstance(library[0], Mapping):
            normalized: list[tuple[str, str]] = []
            for index, item in enumerate(library):
                if not isinstance(item, Mapping):
                    raise TypeError("library mappings cannot be mixed with other entry types")
                target_id = next(
                    (
                        item[key]
                        for key in ("target_id", "id", "guide_id", "feature_id", "name")
                        if key in item and str(item[key]).strip()
                    ),
                    f"target_{index}",
                )
                sequence = next(
                    (
                        item[key]
                        for key in ("target_seq", "sequence", "seq", "guide_seq", "feature_seq")
                        if key in item and str(item[key]).strip()
                    ),
                    None,
                )
                if sequence is None:
                    raise ValueError(f"library mapping at position {index} has no sequence field")
                normalized.append((str(target_id), str(sequence)))
            return normalized
        else:
            return [(f"target_{i}", str(s)) for i, s in enumerate(library)]
    raise TypeError("library must be path, DataFrame, or list of (id,seq) / seqs")


def assign_features(
    adata: Any,
    *,
    library: Any,
    seq_col: str,
    k: int = 1,
    metric: str = "hamming",
    ambiguity_policy: str = "radius",
    result_col: str = "assigned_feature",
    distance_col: str = "feature_distance",
    status_col: str = "feature_status",
    copy: bool = False,
) -> Any | None:
    """
    Assign feature / guide / barcode sequences stored in adata.obs[seq_col]
    against a known library. Adds result columns to .obs.

    This is the in-place equivalent of using Matcher + assign_dataframe,
    but with AnnData ergonomics.

    Parameters
    ----------
    adata : AnnData
        Must have the window sequences in a column of .obs.
    library : str | Path | DataFrame | list
        Target library (id + sequence). Same formats accepted by targets_from_dataframe.
    seq_col : str
        Column in adata.obs containing the extracted target window sequences.
    k, metric, ambiguity_policy : assignment parameters (passed to core).
    result_col, distance_col, status_col : names for the columns written to .obs.
    copy : bool
        If True, return a copy of adata instead of modifying in place.

    Returns
    -------
    AnnData or None
        The (modified) AnnData. None if copy=False (modified in place).

    Notes
    -----
    Only 'unique' assignments receive the target id in result_col.
    For ambiguous / none, result_col will be empty or NaN and status_col will reflect it.
    This prevents silent assignment of ambiguous reads.
    """
    _ensure_anndata()
    if copy:
        adata = adata.copy()

    if seq_col not in adata.obs.columns:
        raise KeyError(f"seq_col '{seq_col}' not found in adata.obs")

    seqs = adata.obs[seq_col].astype(str).tolist()
    targets = _load_library(library)

    if metric not in {"hamming", "levenshtein", "exact"}:
        raise ValueError("metric must be 'hamming', 'levenshtein', or 'exact'")
    if metric == "exact" and k != 0:
        raise ValueError("metric='exact' requires k=0")

    # Use efficient indexed path
    with Matcher([t[1] for t in targets]) as matcher:
        if metric == "hamming":
            results = matcher.assign_hamming(seqs, k=k, policy=ambiguity_policy)
        elif metric == "exact":
            results = matcher.assign_exact(seqs, policy=ambiguity_policy)
        else:
            results = matcher.assign(seqs, k=k, policy=ambiguity_policy)

    # Map back
    id_map = {i: t[0] for i, t in enumerate(targets)}
    assigned = []
    dists = []
    statuses = []

    for r in results:
        if r.status == 1 and r.target_index >= 0:  # UNIQUE
            assigned.append(id_map.get(r.target_index, ""))
            dists.append(r.best_distance)
            statuses.append("unique")
        else:
            assigned.append("")
            dists.append(r.best_distance if r.best_distance >= 0 else -1)
            if r.status == 2:
                statuses.append("ambiguous")
            elif r.status == 0:
                statuses.append("none")
            else:
                statuses.append("invalid")

    adata.obs[result_col] = assigned
    adata.obs[distance_col] = dists
    adata.obs[status_col] = statuses

    # Also store the full policy used for reproducibility
    adata.uns["dotmatch_tl"] = adata.uns.get("dotmatch_tl", {})
    adata.uns["dotmatch_tl"]["last_assign"] = {
        "seq_col": seq_col,
        "k": k,
        "metric": metric,
        "ambiguity_policy": ambiguity_policy,
        "library_size": len(targets),
    }

    return adata if copy else None


def feature_counts(
    adata: Any,
    *,
    library: Any,
    seq_col: str,
    cell_col: str | None = None,
    k: int = 1,
    metric: str = "hamming",
    ambiguity_policy: str = "radius",
    **kwargs,
) -> Any:
    """
    High-level convenience: from an AnnData that contains per-"read" (or per-UMI)
    feature sequences, produce a cells x features count AnnData.

    Internally uses assign + the assignments_to_anndata helper, guaranteeing
    that only uniquely assigned reads contribute to counts.

    This is ideal for building the feature count matrix for CITE-seq, hashing,
    guide capture, etc. that can then be used with scanpy or concatenated with
    the gene expression AnnData.

    Returns a new AnnData with .X = counts, .var = features, .obs = cells.
    Additional QC columns (ambiguous_count etc.) can be requested via kwargs
    passed to assignments_to_anndata.
    """
    _ensure_anndata()

    seqs = adata.obs[seq_col].astype(str).tolist()

    # Build a temporary assignments table using the core assign
    targets = _load_library(library)
    target_seqs = [t[1] for t in targets]
    target_names = [t[0] for t in targets]

    if metric not in {"hamming", "levenshtein", "exact"}:
        raise ValueError("metric must be 'hamming', 'levenshtein', or 'exact'")
    if metric == "exact" and k != 0:
        raise ValueError("metric='exact' requires k=0")

    with Matcher(target_seqs) as matcher:
        if metric == "hamming":
            results = matcher.assign_hamming(seqs, k=k, policy=ambiguity_policy)
        elif metric == "exact":
            results = matcher.assign_exact(seqs, policy=ambiguity_policy)
        else:
            results = matcher.assign(seqs, k=k, policy=ambiguity_policy)

    # Turn results into a DataFrame that assignments_to_anndata understands
    import pandas as _pd  # local to avoid top-level dep if not needed

    tmp = _pd.DataFrame(
        {
            "read_index": range(len(results)),
            "target_name": [
                target_names[r.target_index] if r.target_index >= 0 else ""
                for r in results
            ],
            "status_name": [
                {1: "unique", 2: "ambiguous", 0: "none", -1: "invalid"}.get(r.status, "invalid")
                for r in results
            ],
        }
    )

    if cell_col and cell_col in adata.obs.columns:
        tmp[cell_col] = adata.obs[cell_col].values
    else:
        # fallback: use index as cell if not provided
        tmp[cell_col or "cell"] = adata.obs.index.astype(str).values
        cell_col = cell_col or "cell"

    # Delegate to the battle-tested helper (enforces unique-only)
    return assignments_to_anndata(
        tmp,
        cell_col=cell_col,
        feature_col="target_name",
        count_unique_only=True,
        **kwargs,
    )


# Convenience aliases for common assay types
crispr_guide_assignment = assign_features
feature_barcode_assignment = assign_features
