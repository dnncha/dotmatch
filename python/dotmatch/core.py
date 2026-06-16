from __future__ import annotations

import ctypes
import os
import platform
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence, Any

try:
    import pandas as pd  # type: ignore
    _HAS_PANDAS = True
except Exception:  # noqa: BLE001
    pd = None  # type: ignore
    _HAS_PANDAS = False

try:
    import polars as pl  # type: ignore
    _HAS_POLARS = True
except Exception:  # noqa: BLE001
    pl = None  # type: ignore
    _HAS_POLARS = False

try:
    import anndata as ad  # type: ignore
    _HAS_ANNDATA = True
except Exception:  # noqa: BLE001
    ad = None  # type: ignore
    _HAS_ANNDATA = False

MATCH_INVALID = -1
MATCH_NONE = 0
MATCH_UNIQUE = 1
MATCH_AMBIGUOUS = 2


class _CMatchResult(ctypes.Structure):
    _fields_ = [
        ("target_index", ctypes.c_int),
        ("best_distance", ctypes.c_int),
        ("second_best_distance", ctypes.c_int),
        ("match_count", ctypes.c_int),
        ("status", ctypes.c_int),
    ]


class _CIndexStats(ctypes.Structure):
    _fields_ = [
        ("candidates_considered", ctypes.c_size_t),
        ("candidates_verified", ctypes.c_size_t),
    ]


@dataclass(frozen=True)
class MatchResult:
    target_index: int
    best_distance: int
    second_best_distance: int
    match_count: int
    status: int


@dataclass(frozen=True)
class AssignmentStats:
    candidates_considered: int
    candidates_verified: int


def _platform_ext() -> str:
    return "dylib" if platform.system() == "Darwin" else "so"


def _candidate_paths() -> list[Path]:
    env = os.environ.get("DOTMATCH_LIB") or os.environ.get("QUICKDNA_LIB")
    paths = [Path(env)] if env else []
    here = Path(__file__).resolve()
    ext = _platform_ext()
    names = [f"libdotmatch.{ext}", f"libqdalign.{ext}"]
    for name in names:
        paths.extend(
            [
                here.parent / name,
                here.parents[2] / name,
                Path.cwd() / name,
            ]
        )
    return paths


def _load_lib() -> ctypes.CDLL:
    for path in _candidate_paths():
        if path.exists():
            lib = ctypes.CDLL(str(path))
            lib.qdaln_alphabet_policy.argtypes = []
            lib.qdaln_alphabet_policy.restype = ctypes.c_char_p
            lib.qdaln_edit_distance.argtypes = [
                ctypes.c_char_p,
                ctypes.c_size_t,
                ctypes.c_char_p,
                ctypes.c_size_t,
            ]
            lib.qdaln_edit_distance.restype = ctypes.c_int
            lib.qdaln_edit_distance_leq.argtypes = [
                ctypes.c_char_p,
                ctypes.c_size_t,
                ctypes.c_char_p,
                ctypes.c_size_t,
                ctypes.c_int,
            ]
            lib.qdaln_edit_distance_leq.restype = ctypes.c_int
            lib.qdaln_match_many.argtypes = [
                ctypes.POINTER(ctypes.c_char_p),
                ctypes.POINTER(ctypes.c_size_t),
                ctypes.c_size_t,
                ctypes.POINTER(ctypes.c_char_p),
                ctypes.POINTER(ctypes.c_size_t),
                ctypes.c_size_t,
                ctypes.c_int,
                ctypes.POINTER(_CMatchResult),
            ]
            lib.qdaln_match_many.restype = ctypes.c_int
            lib.qdaln_index_build.argtypes = [
                ctypes.POINTER(ctypes.c_char_p),
                ctypes.POINTER(ctypes.c_size_t),
                ctypes.c_size_t,
            ]
            lib.qdaln_index_build.restype = ctypes.c_void_p
            lib.qdaln_index_free.argtypes = [ctypes.c_void_p]
            lib.qdaln_index_free.restype = None
            lib.qdaln_index_assign_stats.argtypes = [
                ctypes.c_void_p,
                ctypes.POINTER(ctypes.c_char_p),
                ctypes.POINTER(ctypes.c_size_t),
                ctypes.c_size_t,
                ctypes.c_int,
                ctypes.POINTER(_CMatchResult),
                ctypes.POINTER(_CIndexStats),
            ]
            lib.qdaln_index_assign_stats.restype = ctypes.c_int
            return lib
    searched = ", ".join(str(p) for p in _candidate_paths())
    raise RuntimeError(f"could not find DotMatch native library; searched: {searched}")


_LIB = _load_lib()


def alphabet_policy() -> str:
    policy = _LIB.qdaln_alphabet_policy()
    if policy is None:
        raise RuntimeError("DotMatch native library returned no alphabet policy")
    return policy.decode("ascii")


def _as_bytes(seq: str | bytes) -> bytes:
    if isinstance(seq, bytes):
        return seq
    if isinstance(seq, str):
        return seq.encode("ascii")
    raise TypeError("sequence must be str or bytes")


def distance(a: str | bytes, b: str | bytes) -> int:
    aa = _as_bytes(a)
    bb = _as_bytes(b)
    result = int(_LIB.qdaln_edit_distance(aa, len(aa), bb, len(bb)))
    if result < 0:
        raise ValueError("invalid sequence input")
    return result


def distance_leq(a: str | bytes, b: str | bytes, k: int) -> bool:
    aa = _as_bytes(a)
    bb = _as_bytes(b)
    result = int(_LIB.qdaln_edit_distance_leq(aa, len(aa), bb, len(bb), int(k)))
    if result < 0:
        raise ValueError("invalid sequence input")
    return bool(result)


def _array_inputs(seqs: Sequence[str | bytes]) -> tuple[list[bytes], ctypes.Array, ctypes.Array]:
    encoded = [_as_bytes(s) for s in seqs]
    ptrs = (ctypes.c_char_p * len(encoded))()
    lens = (ctypes.c_size_t * len(encoded))()
    for i, seq in enumerate(encoded):
        ptrs[i] = seq
        lens[i] = len(seq)
    return encoded, ptrs, lens


def _results_to_python(results: ctypes.Array) -> list[MatchResult]:
    return [
        MatchResult(
            target_index=r.target_index,
            best_distance=r.best_distance,
            second_best_distance=r.second_best_distance,
            match_count=r.match_count,
            status=r.status,
        )
        for r in results
    ]


def _normalize_policy(policy: str) -> str:
    if policy not in {"radius", "best"}:
        raise ValueError("policy must be 'radius' or 'best'")
    return policy


def _apply_policy(results: list[MatchResult], policy: str) -> list[MatchResult]:
    policy = _normalize_policy(policy)
    if policy == "best":
        return results
    return [
        MatchResult(
            target_index=r.target_index,
            best_distance=r.best_distance,
            second_best_distance=r.second_best_distance,
            match_count=r.match_count,
            status=MATCH_AMBIGUOUS,
        )
        if r.status == MATCH_UNIQUE and r.match_count > 1
        else r
        for r in results
    ]


def assign(
    reads: Sequence[str | bytes],
    barcodes: Sequence[str | bytes],
    k: int = 1,
    policy: str = "radius",
) -> list[MatchResult]:
    if k < 0:
        raise ValueError("k must be non-negative")
    _normalize_policy(policy)
    _read_bytes, read_ptrs, read_lens = _array_inputs(reads)
    _target_bytes, target_ptrs, target_lens = _array_inputs(barcodes)
    results = (_CMatchResult * len(reads))()

    rc = int(
        _LIB.qdaln_match_many(
            read_ptrs,
            read_lens,
            len(reads),
            target_ptrs,
            target_lens,
            len(barcodes),
            int(k),
            results,
        )
    )
    if rc != 0:
        raise ValueError("invalid batch assignment input")

    return _apply_policy(_results_to_python(results), policy)


class Matcher:
    def __init__(self, barcodes: Sequence[str | bytes]):
        self._closed = False
        self._target_bytes, target_ptrs, target_lens = _array_inputs(barcodes)
        self._index = _LIB.qdaln_index_build(target_ptrs, target_lens, len(barcodes))
        if not self._index:
            raise ValueError("invalid barcode input")

    def close(self) -> None:
        if not self._closed:
            _LIB.qdaln_index_free(self._index)
            self._index = None
            self._closed = True

    def __enter__(self) -> "Matcher":
        if self._closed:
            raise ValueError("matcher is closed")
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def assign(self, reads: Sequence[str | bytes], k: int = 1, policy: str = "radius") -> list[MatchResult]:
        results, _stats = self.assign_with_stats(reads, k=k, policy=policy)
        return results

    def assign_with_stats(
        self,
        reads: Sequence[str | bytes],
        k: int = 1,
        policy: str = "radius",
    ) -> tuple[list[MatchResult], AssignmentStats]:
        if self._closed:
            raise ValueError("matcher is closed")
        if k < 0:
            raise ValueError("k must be non-negative")
        _normalize_policy(policy)

        _read_bytes, read_ptrs, read_lens = _array_inputs(reads)
        results = (_CMatchResult * len(reads))()
        stats = _CIndexStats()
        rc = int(
            _LIB.qdaln_index_assign_stats(
                self._index,
                read_ptrs,
                read_lens,
                len(reads),
                int(k),
                results,
                ctypes.byref(stats),
            )
        )
        if rc != 0:
            raise ValueError("invalid indexed assignment input")

        return _apply_policy(_results_to_python(results), policy), AssignmentStats(
            candidates_considered=int(stats.candidates_considered),
            candidates_verified=int(stats.candidates_verified),
        )


def _ensure_pandas() -> None:
    if not _HAS_PANDAS:
        raise ImportError(
            "pandas is required; install with 'pip install \"dotmatch[pandas]\"' (or add pandas to your environment)"
        )


def _to_pandas(df: Any) -> Any:
    """Convert polars or other to pandas for internal use; return as-is if already pandas."""
    if _HAS_POLARS and isinstance(df, pl.DataFrame):
        return df.to_pandas()
    if _HAS_PANDAS and hasattr(df, "to_pandas"):
        # polars lazy etc, but simple
        try:
            return df.to_pandas()
        except Exception:
            pass
    return df

def targets_from_dataframe(
    df: Any,
    id_col: str | int | None = None,
    seq_col: str | int | None = None,
) -> list[tuple[str, str]]:
    """Extract (id, sequence) pairs from a pandas or polars DataFrame or Series-like.

    If id_col/seq_col omitted, assumes first two columns (or column named 'id'/'seq'/'sequence'/'barcode').
    Returns list of (id, seq) suitable for Matcher or assign.
    Supports polars DataFrames (converted internally).
    """
    if _HAS_PANDAS:
        _ensure_pandas()
    else:
        if not (_HAS_POLARS and (isinstance(df, pl.DataFrame) if pl else False or hasattr(df, "to_pandas"))):
            raise ImportError("pandas (or polars) required for DataFrame input to targets_from_dataframe")
    data = _to_pandas(df)
    if _HAS_PANDAS:
        cols = list(data.columns)
    else:
        cols = list(getattr(data, "columns", []))
    if id_col is None:
        for cand in ("id", "target_id", "name", "guide_id", "barcode_id", cols[0] if cols else None):
            if cand in cols:
                id_col = cand
                break
        else:
            id_col = cols[0] if cols else 0
    if seq_col is None:
        for cand in ("seq", "sequence", "dna", "target", "guide", "barcode", cols[1] if len(cols) > 1 else None):
            if cand in cols:
                seq_col = cand
                break
        else:
            seq_col = cols[1] if len(cols) > 1 else 1
    ids = data[id_col].astype(str).tolist()
    seqs = data[seq_col].astype(str).tolist()
    return list(zip(ids, seqs))


def results_to_dataframe(
    results: Sequence[MatchResult],
    target_names: Sequence[str] | None = None,
    read_ids: Sequence[str] | None = None,
) -> Any:
    """Convert MatchResult list to pandas DataFrame if pandas available.

    Adds human columns: status_name, target_name (if names provided), etc.
    Returns DataFrame or raises if no pandas.
    """
    _ensure_pandas()
    status_map = {
        MATCH_UNIQUE: "unique",
        MATCH_AMBIGUOUS: "ambiguous",
        MATCH_NONE: "none",
        MATCH_INVALID: "invalid",
    }
    rows = []
    for i, r in enumerate(results):
        row = {
            "read_index": i,
            "target_index": r.target_index,
            "best_distance": r.best_distance,
            "second_best_distance": r.second_best_distance,
            "match_count": r.match_count,
            "status": r.status,
            "status_name": status_map.get(r.status, str(r.status)),
        }
        if target_names is not None and 0 <= r.target_index < len(target_names):
            row["target_name"] = target_names[r.target_index]
        if read_ids is not None and i < len(read_ids):
            row["read_id"] = read_ids[i]
        rows.append(row)
    return pd.DataFrame(rows)


def assign_dataframe(
    reads: Any,
    targets: Any,
    k: int = 1,
    policy: str = "radius",
    read_ids: Sequence[str] | None = None,
    target_names: Sequence[str] | None = None,
) -> Any:
    """High-level: assign using pandas/polars Series/DataFrame inputs, return pandas DataFrame of results.

    reads/targets can be list, Series of seqs, or DataFrame (will use seq col heuristics).
    Polars inputs are converted internally; result is pandas DataFrame (call .to_polars() if desired).
    """
    if _HAS_PANDAS:
        _ensure_pandas()
    targets = _to_pandas(targets)
    reads = _to_pandas(reads)
    # normalize targets to seq list
    if _HAS_PANDAS and hasattr(targets, "iloc"):
        if len(getattr(targets, "shape", (0,))) > 1 and targets.shape[1] > 0:  # df
            tseqs = targets.iloc[:, 1].astype(str).tolist() if targets.shape[1] > 1 else targets.iloc[:, 0].astype(str).tolist()
            tnames = targets.iloc[:, 0].astype(str).tolist() if target_names is None and targets.shape[1] > 0 else (target_names or None)
        else:
            tseqs = targets.astype(str).tolist()
            tnames = target_names
    else:
        tseqs = [str(x) for x in targets]
        tnames = target_names
    # reads
    if _HAS_PANDAS and hasattr(reads, "iloc"):
        rseqs = reads.astype(str).tolist() if hasattr(reads, "astype") else [str(x) for x in reads]
        rids = read_ids or (reads.index.astype(str).tolist() if hasattr(reads, "index") else None)
    else:
        rseqs = [str(x) for x in reads]
        rids = read_ids
    res = assign(rseqs, tseqs, k=k, policy=policy)
    return results_to_dataframe(res, target_names=tnames, read_ids=rids)


def _ensure_anndata() -> None:
    if not _HAS_ANNDATA:
        raise ImportError(
            "anndata is required; install with 'pip install \"dotmatch[anndata]\"' "
            "(pulls in pandas too). For full scanpy workflows, also install scanpy."
        )
    _ensure_pandas()


def counts_tsv_to_anndata(
    counts_path: str | Path,
    *,
    sample_cols: list[str] | None = None,
    var_cols: list[str] = ("target_id", "target_seq", "gene"),
) -> Any:
    """Load a DotMatch counts TSV (mageck or dotmatch format) into an AnnData object.

    The resulting AnnData has:
    - X: counts (cells/samples x features) as sparse or dense
    - var: feature metadata (target_id, seq, gene if present)
    - obs: samples/cells
    Useful bridge after running `dotmatch count --format mageck ...` or the Python CLI.
    """
    _ensure_anndata()
    _ensure_pandas()
    path = Path(counts_path)
    # Reuse the robust parser from cli if possible, else simple pandas read
    try:
        # Try to use internal if exposed, else fall back
        from . import cli as _cli  # type: ignore
        counts_dict = _cli._read_crispr_count_matrix(str(path))  # may be private
        # Reconstruct simple matrix
        guides = counts_dict.get("guides", [])
        samples = counts_dict.get("samples", [])
        data = []
        for g in guides:
            row = [g["counts"].get(s, 0) for s in samples]
            data.append(row)
        X = pd.DataFrame(data, index=[g["id"] for g in guides], columns=samples).T
        var_df = pd.DataFrame(
            [{"target_id": g["id"], "gene": g.get("gene", "")} for g in guides]
        ).set_index("target_id")
    except Exception:
        # Simple fallback: read tsv, assume first cols are id/gene, rest numeric samples
        df = pd.read_csv(path, sep="\t")
        # Try to detect id col
        id_col = None
        for c in ("target_id", "guide", "id", "feature_id", df.columns[0]):
            if c in df.columns:
                id_col = c
                break
        numeric = df.select_dtypes(include="number")
        if numeric.shape[1] == 0:
            numeric = df.iloc[:, 2:] if df.shape[1] > 2 else df.iloc[:, 1:]
        X = numeric.T
        X.columns = df[id_col].astype(str).values if id_col else [f"f{i}" for i in range(X.shape[1])]
        var_df = pd.DataFrame(index=X.columns)
        for c in var_cols:
            if c in df.columns and c != id_col:
                var_df[c] = df.set_index(id_col)[c].reindex(X.columns).values if id_col else None
    adata = ad.AnnData(X=X.values if hasattr(X, "values") else X)
    adata.var = var_df if len(var_df) else pd.DataFrame(index=X.columns)
    adata.obs = pd.DataFrame(index=X.index.astype(str))
    adata.var_names = X.columns.astype(str)
    adata.obs_names = X.index.astype(str)
    return adata


def assignments_to_anndata(
    assignments: Any,
    *,
    cell_col: str = "cell_barcode",
    feature_col: str = "target_name",
    count_unique_only: bool = True,
    include_ambiguous_per_cell: bool = False,
) -> Any:
    """Build a cells x features count AnnData from a per-read assignments table (pandas/polars DF or path to assignments.tsv).

    This is useful for 10x-style feature barcode or guide capture / perturb-seq where you have
    pre-extracted cell barcodes + the feature/guide window sequence, ran assignment (via CLI `dotmatch count --assignments` or Python `assign_dataframe`),
    and now want a count matrix in AnnData form for scanpy/pertpy/etc. downstream analysis.

    By default (count_unique_only=True) only status==unique reads contribute to the count matrix.
    This preserves DotMatch's core scientific contract: ambiguous reads are never silently assigned.

    If your assignments came from the CLI, the read_id often encodes cell info (parse or supply cell_col).
    """
    _ensure_anndata()
    _ensure_pandas()
    if isinstance(assignments, (str, Path)):
        df = pd.read_csv(assignments, sep="\t")
    else:
        df = _to_pandas(assignments)
    if cell_col not in df.columns:
        # heuristic: try common names or assume first col or read_id contains it
        for cand in ("cell_barcode", "cell", "barcode", "CB", "cell_id", "barcode_id"):
            if cand in df.columns:
                cell_col = cand
                break
        else:
            if "read_id" in df.columns:
                # naive; real pipelines pre-extract or use proper CB from 10x R1
                df[cell_col] = df["read_id"].astype(str).str.split("_").str[0]
            else:
                raise ValueError(f"Could not find cell column '{cell_col}'; pass cell_col= or pre-populate the DF")
    if feature_col not in df.columns:
        for cand in ("target_name", "target_id", "feature", "guide", "id", "target"):
            if cand in df.columns:
                feature_col = cand
                break
        else:
            feature_col = "target_index"  # fallback

    # Always compute per-cell stats for QC/accuracy visibility
    if "status_name" in df.columns:
        unique_mask = df["status_name"].isin(["unique"])
        ambig_mask = df["status_name"].isin(["ambiguous"])
    else:
        unique_mask = df.get("status", 0) == 1
        ambig_mask = df.get("status", 0) == 2

    if count_unique_only:
        df_unique = df[unique_mask]
    else:
        df_unique = df[unique_mask]  # still only uniques for matrix

    # aggregate unique counts
    grp = df_unique.groupby([cell_col, feature_col]).size().reset_index(name="count")
    pivot = grp.pivot(index=cell_col, columns=feature_col, values="count").fillna(0).astype(int)

    adata = ad.AnnData(X=pivot.values)
    adata.obs = pd.DataFrame(index=pivot.index.astype(str))
    adata.var = pd.DataFrame(index=pivot.columns.astype(str))
    adata.var_names = pivot.columns.astype(str)
    adata.obs_names = pivot.index.astype(str)
    adata.layers["counts"] = adata.X.copy()

    if include_ambiguous_per_cell:
        ambig_counts = df[ambig_mask].groupby(cell_col).size()
        adata.obs["ambiguous_count"] = adata.obs_names.map(ambig_counts).fillna(0).astype(int)
        adata.obs["unique_count"] = adata.obs_names.map(
            df_unique.groupby(cell_col).size()
        ).fillna(0).astype(int)

    # Reproducibility / scientific metadata (policy, k etc. can be joined from summary if user passes more context)
    adata.uns["dotmatch"] = {"source": "assignments_to_anndata", "unique_only": count_unique_only}
    return adata
