from __future__ import annotations

import ctypes
import csv
import gzip
import math
import os
import platform
import threading
from functools import wraps
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence, TextIO

from ._optional import optional_module
from ._validation import integer, named_targets
from .fastq_io import FastqRecord, iter_fastq

pd, _HAS_PANDAS = optional_module("pandas")

pl, _HAS_POLARS = optional_module("polars")

ad, _HAS_ANNDATA = optional_module("anndata")

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


@dataclass(frozen=True)
class PosteriorAssignment:
    target_index: int
    posterior: float
    second_posterior: float
    status: int
    posteriors: tuple[float, ...]


@dataclass(frozen=True)
class StreamAssignment:
    read_id: str
    observed_seq: str
    target_index: int
    target_name: str
    target_seq: str
    best_distance: int
    second_best_distance: int
    match_count: int
    status: int
    status_name: str


def _platform_ext() -> str:
    return "dylib" if platform.system() == "Darwin" else "so"


def _candidate_paths() -> list[Path]:
    # Overrides are explicit, never a hint followed by an unrelated fallback.
    env = os.environ.get("DOTMATCH_LIB") or os.environ.get("QUICKDNA_LIB")
    if env:
        return [Path(env).expanduser().resolve()]
    here = Path(__file__).resolve()
    names = [f"libdotmatch.{_platform_ext()}", f"libqdalign.{_platform_ext()}"]
    paths = [here.parent / name for name in names]
    # Only the known source layout permits a source-tree build. Never search CWD
    # or arbitrary ancestors of an installed/vendored package for executable code.
    if here.parent.name == "dotmatch" and here.parent.parent.name == "python":
        paths.extend(here.parents[2] / name for name in names)
    return paths


def _load_lib() -> ctypes.CDLL:
    for path in _candidate_paths():
        if path.is_file():
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
            lib.qdaln_index_assign_status_stats.argtypes = [
                ctypes.c_void_p,
                ctypes.POINTER(ctypes.c_char_p),
                ctypes.POINTER(ctypes.c_size_t),
                ctypes.c_size_t,
                ctypes.c_int,
                ctypes.POINTER(_CMatchResult),
                ctypes.POINTER(_CIndexStats),
            ]
            lib.qdaln_index_assign_status_stats.restype = ctypes.c_int
            lib.qdaln_index_assign_hamming_stats.argtypes = [
                ctypes.c_void_p,
                ctypes.POINTER(ctypes.c_char_p),
                ctypes.POINTER(ctypes.c_size_t),
                ctypes.c_size_t,
                ctypes.c_int,
                ctypes.POINTER(_CMatchResult),
                ctypes.POINTER(_CIndexStats),
            ]
            lib.qdaln_index_assign_hamming_stats.restype = ctypes.c_int
            lib.qdaln_index_lookup_exact_many_stats.argtypes = [
                ctypes.c_void_p,
                ctypes.POINTER(ctypes.c_char_p),
                ctypes.POINTER(ctypes.c_size_t),
                ctypes.c_size_t,
                ctypes.POINTER(_CMatchResult),
                ctypes.POINTER(_CIndexStats),
            ]
            lib.qdaln_index_lookup_exact_many_stats.restype = ctypes.c_int
            lib.qdaln_index_lookup_exact_ascii_many_stats.argtypes = [
                ctypes.c_void_p,
                ctypes.POINTER(ctypes.c_char_p),
                ctypes.POINTER(ctypes.c_size_t),
                ctypes.c_size_t,
                ctypes.POINTER(_CMatchResult),
                ctypes.POINTER(_CIndexStats),
            ]
            lib.qdaln_index_lookup_exact_ascii_many_stats.restype = ctypes.c_int
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


def status_name(status: int) -> str:
    """Return the stable text name for a DotMatch assignment status."""
    return _STATUS_NAMES.get(status, f"unknown:{status}")


_STATUS_NAMES = {
        MATCH_INVALID: "invalid",
        MATCH_NONE: "none",
        MATCH_UNIQUE: "unique",
        MATCH_AMBIGUOUS: "ambiguous",
}


def _open_text(path: str | Path, mode: str = "rt") -> TextIO:
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def _looks_like_target_header(cols: Sequence[str]) -> bool:
    normalized = {str(col).strip().lower() for col in cols[:3]}
    return bool(normalized & {"target_id", "guide_id", "barcode_id", "id", "name"}) and bool(
        normalized & {"target_seq", "guide_seq", "barcode_seq", "sequence", "seq"}
    )


def load_targets(path: str | Path) -> list[tuple[str, str]]:
    """Read plain/gzipped CSV or TSV with named or positional target columns.

    Duplicate sequences are retained; empty or duplicate IDs are rejected.
    """
    from .target_io import read_target_table

    return [(row.target_id, row.sequence) for row in read_target_table(path)]


def _normalize_targets(targets: Any) -> list[tuple[str, str]]:
    if isinstance(targets, (str, Path)):
        return load_targets(targets)
    if hasattr(targets, "columns"):
        return targets_from_dataframe(targets)
    rows = []
    for i, item in enumerate(targets):
        if isinstance(item, (tuple, list)) and len(item) >= 2:
            rows.append((item[0], item[1]))
        else:
            rows.append((f"target_{i}", item))
    return named_targets(rows)


def _extract_window(seq: str, start: int, length: int) -> str | None:
    if start < 0:
        raise ValueError("target_start must be non-negative")
    if length <= 0:
        raise ValueError("target_length must be positive")
    end = start + length
    if end > len(seq):
        return None
    return seq[start:end]


def _chunks(items: Iterable[FastqRecord], size: int) -> Iterator[list[FastqRecord]]:
    if size <= 0:
        raise ValueError("batch_size must be positive")
    batch: list[FastqRecord] = []
    for item in items:
        batch.append(item)
        if len(batch) == size:
            yield batch
            batch = []
    if batch:
        yield batch


def distance(a: str | bytes, b: str | bytes) -> int:
    aa = _as_bytes(a)
    bb = _as_bytes(b)
    result = int(_LIB.qdaln_edit_distance(aa, len(aa), bb, len(bb)))
    if result < 0:
        raise ValueError("invalid sequence input")
    return result


def distance_leq(a: str | bytes, b: str | bytes, k: int) -> bool:
    k = integer(k, "k")
    aa = _as_bytes(a)
    bb = _as_bytes(b)
    result = int(_LIB.qdaln_edit_distance_leq(aa, len(aa), bb, len(bb), int(k)))
    if result < 0:
        raise ValueError("invalid sequence input")
    return bool(result)


def _array_inputs(seqs: Sequence[str | bytes]) -> tuple[list[bytes], ctypes.Array, ctypes.Array]:
    if isinstance(seqs, (str, bytes, bytearray)) or not hasattr(seqs, "__len__"):
        raise TypeError("expected a sized sequence of sequences; wrap a single sequence in a list")
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
    applied: list[MatchResult] | None = None
    for i, r in enumerate(results):
        if r.status == MATCH_UNIQUE and r.match_count > 1:
            if applied is None:
                applied = results[:i]
            applied.append(
                MatchResult(
                    target_index=r.target_index,
                    best_distance=r.best_distance,
                    second_best_distance=r.second_best_distance,
                    match_count=r.match_count,
                    status=MATCH_AMBIGUOUS,
                )
            )
        elif applied is not None:
            applied.append(r)
    return applied if applied is not None else results


def assign(
    reads: Sequence[str | bytes],
    barcodes: Sequence[str | bytes],
    k: int = 1,
    policy: str = "radius",
) -> list[MatchResult]:
    k = integer(k, "k")
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


def assign_hamming(
    reads: Sequence[str | bytes],
    barcodes: Sequence[str | bytes],
    k: int = 1,
    policy: str = "radius",
) -> list[MatchResult]:
    """Assign fixed-length reads by Hamming distance using the native Hamming index."""
    with Matcher(barcodes) as matcher:
        return matcher.assign_hamming(reads, k=k, policy=policy)


def assign_exact(
    reads: Sequence[str | bytes],
    barcodes: Sequence[str | bytes],
    *,
    ascii_fold: bool = False,
    policy: str = "radius",
) -> list[MatchResult]:
    """Assign exact fixed windows using the native exact lookup table."""
    with Matcher(barcodes) as matcher:
        return matcher.assign_exact(reads, ascii_fold=ascii_fold, policy=policy)


def _phred33_probability(ch: int) -> float:
    q = ch - 33
    if q < 0 or q > 93:
        raise ValueError("quality string must use Phred+33 ASCII characters 33–126")
    return 10.0 ** (-q / 10.0)


def assign_posterior(
    read: str | bytes,
    targets: Sequence[str | bytes],
    quality: str | bytes,
    *,
    min_posterior: float = 0.95,
    priors: Sequence[float] | None = None,
) -> PosteriorAssignment:
    """Assign one fixed-window read using a simple Phred likelihood model.

    This is an auditable posterior call helper, not a replacement for the fast
    indexed batch matcher. Bases compare literally; ambiguity symbols are not
    expanded as wildcards.
    """
    read_b = _as_bytes(read)
    qual_b = _as_bytes(quality)
    target_b = [_as_bytes(t) for t in targets]
    if not read_b:
        raise ValueError("posterior assignment requires a nonempty read")
    if not target_b:
        raise ValueError("targets must not be empty")
    if len(read_b) != len(qual_b):
        raise ValueError("read and quality must have the same length")
    if any(len(t) != len(read_b) for t in target_b):
        return PosteriorAssignment(-1, 0.0, 0.0, MATCH_INVALID, tuple())
    if not (0.0 <= min_posterior <= 1.0):
        raise ValueError("min_posterior must be between 0 and 1")

    if priors is None:
        log_priors = [-math.log(len(target_b))] * len(target_b)
    else:
        if len(priors) != len(target_b):
            raise ValueError("priors must have one entry per target")
        values = [float(p) for p in priors]
        if any(not math.isfinite(p) or p < 0.0 for p in values) or not any(p > 0.0 for p in values):
            raise ValueError("priors must be finite, non-negative with positive total mass")
        # Normalize in log space below. Summing large finite priors first can
        # overflow; dividing tiny priors by that sum can underflow to zero.
        log_priors = [math.log(p) if p > 0.0 else -math.inf for p in values]

    log_likelihoods: list[float] = []
    for target, log_prior in zip(target_b, log_priors):
        ll = log_prior
        for rb, tb, qb in zip(read_b, target, qual_b):
            p_error = min(max(_phred33_probability(qb), 1e-12), 1.0 - 1e-12)
            ll += math.log1p(-p_error) if rb == tb else math.log(p_error / 3.0)
        log_likelihoods.append(ll)

    max_ll = max(log_likelihoods)
    weights = [math.exp(ll - max_ll) for ll in log_likelihoods]
    total_weight = sum(weights)
    posteriors = tuple(w / total_weight for w in weights)
    order = sorted(range(len(posteriors)), key=lambda i: posteriors[i], reverse=True)
    best = order[0]
    second = posteriors[order[1]] if len(order) > 1 else 0.0
    status = MATCH_UNIQUE if posteriors[best] >= min_posterior and posteriors[best] > second else MATCH_AMBIGUOUS
    return PosteriorAssignment(best, posteriors[best], second, status, posteriors)


def _synchronized(method):
    @wraps(method)
    def call(self, *args, **kwargs):
        # ctypes releases the GIL: close() must never free an in-flight index.
        with self._lock:
            return method(self, *args, **kwargs)
    return call


class Matcher:
    """Reusable native index with synchronized lifetime management.

    Calls on one Matcher are serialized, including close(). Use an independent
    Matcher per worker when parallel assignment is required.
    """
    def __init__(self, barcodes: Sequence[str | bytes]):
        self._lock = threading.RLock()
        self._closed = True
        self._index = None
        self._target_bytes, target_ptrs, target_lens = _array_inputs(barcodes)
        self._index = _LIB.qdaln_index_build(target_ptrs, target_lens, len(barcodes))
        if not self._index:
            raise ValueError("invalid barcode input")
        self._closed = False

    @_synchronized
    def close(self) -> None:
        if not self._closed:
            _LIB.qdaln_index_free(self._index)
            self._index = None
            self._closed = True

    @_synchronized
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
        return self._assign_with_stats_func(_LIB.qdaln_index_assign_stats, reads, k=k, policy=policy)

    def assign_hamming(self, reads: Sequence[str | bytes], k: int = 1, policy: str = "radius") -> list[MatchResult]:
        """Assign fixed-length reads by Hamming distance with the native Hamming kernel.

        This is the fast path for one-mismatch guide or barcode workflows where
        insertions and deletions are intentionally out of scope.
        """
        results, _stats = self.assign_hamming_with_stats(reads, k=k, policy=policy)
        return results

    def assign_hamming_with_stats(
        self,
        reads: Sequence[str | bytes],
        k: int = 1,
        policy: str = "radius",
    ) -> tuple[list[MatchResult], AssignmentStats]:
        k = integer(k, "k", maximum=3)
        if k == 0:
            return self.assign_exact_with_stats(reads, policy=policy)
        return self._assign_with_stats_func(_LIB.qdaln_index_assign_hamming_stats, reads, k=k, policy=policy)

    def assign_exact(
        self,
        reads: Sequence[str | bytes],
        *,
        ascii_fold: bool = False,
        policy: str = "radius",
    ) -> list[MatchResult]:
        """Assign exact fixed windows using the native exact lookup table."""
        results, _stats = self.assign_exact_with_stats(reads, ascii_fold=ascii_fold, policy=policy)
        return results

    @_synchronized
    def assign_exact_with_stats(
        self,
        reads: Sequence[str | bytes],
        *,
        ascii_fold: bool = False,
        policy: str = "radius",
    ) -> tuple[list[MatchResult], AssignmentStats]:
        if self._closed:
            raise ValueError("matcher is closed")
        _normalize_policy(policy)
        _read_bytes, read_ptrs, read_lens = _array_inputs(reads)
        results = (_CMatchResult * len(reads))()
        stats = _CIndexStats()
        func = _LIB.qdaln_index_lookup_exact_ascii_many_stats if ascii_fold else _LIB.qdaln_index_lookup_exact_many_stats
        rc = int(
            func(
                self._index,
                read_ptrs,
                read_lens,
                len(reads),
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

    def assign_status_with_stats(
        self,
        reads: Sequence[str | bytes],
        k: int = 1,
        policy: str = "radius",
    ) -> tuple[list[MatchResult], AssignmentStats]:
        """Assign with an early-stop native path when only status/best target matters.

        For ambiguous calls, ``match_count`` and ``second_best_distance`` may be
        lower-bound values. Use ``assign_with_stats`` when exact ambiguity counts
        are required.
        """
        return self._assign_with_stats_func(_LIB.qdaln_index_assign_status_stats, reads, k=k, policy=policy)

    @_synchronized
    def _assign_with_stats_func(
        self,
        func: Any,
        reads: Sequence[str | bytes],
        *,
        k: int,
        policy: str,
    ) -> tuple[list[MatchResult], AssignmentStats]:
        if self._closed:
            raise ValueError("matcher is closed")
        k = integer(k, "k")
        _normalize_policy(policy)

        _read_bytes, read_ptrs, read_lens = _array_inputs(reads)
        results = (_CMatchResult * len(reads))()
        stats = _CIndexStats()
        rc = int(
            func(
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


def stream_assign(
    fastq_path: str | Path,
    targets: Any,
    *,
    target_start: int = 0,
    target_length: int | None = None,
    k: int = 1,
    metric: str = "levenshtein",
    policy: str = "radius",
    batch_size: int = 4096,
    include_invalid: bool = True,
) -> Iterator[StreamAssignment]:
    """Stream fixed-window FASTQ assignments without loading reads into memory.

    ``targets`` may be a path, a DataFrame accepted by ``targets_from_dataframe``,
    a list of ``(id, seq)`` pairs, or a list of sequences. Only ``unique`` rows
    carry a target id/sequence; ambiguous, none, and invalid windows remain
    explicit rather than being forced into a target.
    """
    k = integer(k, "k", maximum=3 if metric == "hamming" else 2147483647)
    target_start = integer(target_start, "target_start")
    batch_size = integer(batch_size, "batch_size", minimum=1)
    if target_length is not None:
        target_length = integer(target_length, "target_length", minimum=1)
    _normalize_policy(policy)
    if metric not in {"levenshtein", "hamming", "exact"}:
        raise ValueError("metric must be 'levenshtein', 'hamming', or 'exact'")
    if metric == "exact" and k != 0:
        raise ValueError("metric='exact' requires k=0")
    normalized_targets = _normalize_targets(targets)
    target_names = [target_id for target_id, _seq in normalized_targets]
    target_seqs = [seq for _target_id, seq in normalized_targets]
    if target_length is None:
        lengths = {len(seq) for seq in target_seqs}
        if len(lengths) != 1:
            raise ValueError("target_length is required when targets have mixed lengths")
        target_length = lengths.pop()
    _extract_window("", target_start, target_length)

    with Matcher(target_seqs) as matcher:
        for batch in _chunks(iter_fastq(fastq_path), batch_size):
            observed: list[str] = []
            valid_slots: list[int] = []
            output: list[StreamAssignment | None] = [None] * len(batch)
            for slot, record in enumerate(batch):
                window = _extract_window(record.seq, target_start, target_length)
                if window is None:
                    if include_invalid:
                        output[slot] = StreamAssignment(
                            read_id=record.read_id,
                            observed_seq="",
                            target_index=-1,
                            target_name="",
                            target_seq="",
                            best_distance=-1,
                            second_best_distance=-1,
                            match_count=0,
                            status=MATCH_INVALID,
                            status_name=status_name(MATCH_INVALID),
                        )
                    continue
                observed.append(window)
                valid_slots.append(slot)

            if observed:
                if metric == "levenshtein":
                    results = matcher.assign(observed, k=k, policy=policy)
                elif metric == "hamming":
                    results = matcher.assign_hamming(observed, k=k, policy=policy)
                elif metric == "exact":
                    results = matcher.assign_exact(observed, policy=policy)
                for slot, window, result in zip(
                    valid_slots,
                    observed,
                    results,
                ):
                    record = batch[slot]
                    if result.status == MATCH_UNIQUE and 0 <= result.target_index < len(target_names):
                        target_name = target_names[result.target_index]
                        target_seq = target_seqs[result.target_index]
                    else:
                        target_name = ""
                        target_seq = ""
                    output[slot] = StreamAssignment(
                        read_id=record.read_id,
                        observed_seq=window,
                        target_index=result.target_index,
                        target_name=target_name,
                        target_seq=target_seq,
                        best_distance=result.best_distance,
                        second_best_distance=result.second_best_distance,
                        match_count=result.match_count,
                        status=result.status,
                        status_name=status_name(result.status),
                )
            for assignment in output:
                if assignment is not None:
                    yield assignment


def assignment_summary(assignments: Iterable[StreamAssignment]) -> dict[str, int | float]:
    """Summarize streamed assignments into counts and core rates."""
    total_reads = 0
    assigned_unique = 0
    assigned_exact = 0
    assigned_corrected = 0
    ambiguous = 0
    unmatched = 0
    invalid = 0
    for assignment in assignments:
        total_reads += 1
        if assignment.status == MATCH_UNIQUE:
            assigned_unique += 1
            if assignment.best_distance == 0:
                assigned_exact += 1
            else:
                assigned_corrected += 1
        elif assignment.status == MATCH_AMBIGUOUS:
            ambiguous += 1
        elif assignment.status == MATCH_NONE:
            unmatched += 1
        else:
            invalid += 1
    return _summary_from_counts(
        total_reads,
        assigned_unique,
        assigned_exact,
        assigned_corrected,
        ambiguous,
        unmatched,
        invalid,
    )


def _empty_assignment_summary() -> dict[str, int | float]:
    return {
        "total_reads": 0,
        "assigned_unique": 0,
        "assigned_exact": 0,
        "assigned_corrected": 0,
        "ambiguous": 0,
        "unmatched": 0,
        "invalid": 0,
    }


def _add_assignment_to_summary(summary: dict[str, int | float], assignment: StreamAssignment) -> None:
    summary["total_reads"] = int(summary["total_reads"]) + 1
    if assignment.status == MATCH_UNIQUE:
        summary["assigned_unique"] = int(summary["assigned_unique"]) + 1
        if assignment.best_distance == 0:
            summary["assigned_exact"] = int(summary["assigned_exact"]) + 1
        else:
            summary["assigned_corrected"] = int(summary["assigned_corrected"]) + 1
    elif assignment.status == MATCH_AMBIGUOUS:
        summary["ambiguous"] = int(summary["ambiguous"]) + 1
    elif assignment.status == MATCH_NONE:
        summary["unmatched"] = int(summary["unmatched"]) + 1
    else:
        summary["invalid"] = int(summary["invalid"]) + 1


def _summary_from_counts(
    total_reads: int,
    assigned_unique: int,
    assigned_exact: int,
    assigned_corrected: int,
    ambiguous: int,
    unmatched: int,
    invalid: int,
) -> dict[str, int | float]:
    return _finish_assignment_summary(
        {
            "total_reads": total_reads,
            "assigned_unique": assigned_unique,
            "assigned_exact": assigned_exact,
            "assigned_corrected": assigned_corrected,
            "ambiguous": ambiguous,
            "unmatched": unmatched,
            "invalid": invalid,
        }
    )


def _finish_assignment_summary(summary: dict[str, int | float]) -> dict[str, int | float]:
    total = int(summary["total_reads"])
    summary["assignment_rate"] = int(summary["assigned_unique"]) / total if total else 0.0
    summary["ambiguous_rate"] = int(summary["ambiguous"]) / total if total else 0.0
    summary["no_match_rate"] = int(summary["unmatched"]) / total if total else 0.0
    summary["invalid_rate"] = int(summary["invalid"]) / total if total else 0.0
    return summary


def write_assignments_tsv(assignments: Iterable[StreamAssignment], path: str | Path) -> dict[str, int | float]:
    """Write streamed assignments to TSV and return ``assignment_summary``."""
    columns = [
        "read_id",
        "observed_seq",
        "target_id",
        "target_seq",
        "distance",
        "status",
        "match_count",
        "second_best_distance",
    ]
    total_reads = 0
    assigned_unique = 0
    assigned_exact = 0
    assigned_corrected = 0
    ambiguous = 0
    unmatched = 0
    invalid = 0
    with _open_text(path, "wt") as fh:
        write = fh.write
        write("\t".join(columns) + "\n")
        for row in assignments:
            write(
                f"{row.read_id}\t{row.observed_seq}\t{row.target_name}\t{row.target_seq}\t"
                f"{row.best_distance}\t{row.status_name}\t{row.match_count}\t{row.second_best_distance}\n"
            )
            total_reads += 1
            if row.status == MATCH_UNIQUE:
                assigned_unique += 1
                if row.best_distance == 0:
                    assigned_exact += 1
                else:
                    assigned_corrected += 1
            elif row.status == MATCH_AMBIGUOUS:
                ambiguous += 1
            elif row.status == MATCH_NONE:
                unmatched += 1
            else:
                invalid += 1
    return _summary_from_counts(
        total_reads,
        assigned_unique,
        assigned_exact,
        assigned_corrected,
        ambiguous,
        unmatched,
        invalid,
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

def targets_from_dataframe(df: Any, id_col=None, seq_col=None) -> list[tuple[str, str]]:
    """Extract validated (ID, sequence) pairs; see dotmatch.dataframes."""
    from .dataframes import targets_from_dataframe as convert
    return convert(df, id_col=id_col, seq_col=seq_col)


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
    if read_ids is not None and len(read_ids) != len(results):
        raise ValueError("read_ids must contain one ID per result")
    if target_names is not None and any(r.status == MATCH_UNIQUE and not 0 <= r.target_index < len(target_names) for r in results):
        raise ValueError("target_names does not cover every uniquely assigned target index")
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
        if target_names is not None:
            # A candidate index is diagnostic information, not an assignment.
            row["target_name"] = (
                target_names[r.target_index]
                if r.status == MATCH_UNIQUE and 0 <= r.target_index < len(target_names)
                else ""
            )
        if read_ids is not None and i < len(read_ids):
            row["read_id"] = read_ids[i]
        rows.append(row)
    columns = ["read_index", "target_index", "best_distance", "second_best_distance", "match_count", "status", "status_name"]
    if target_names is not None:
        columns.append("target_name")
    if read_ids is not None:
        columns.append("read_id")
    return pd.DataFrame(rows, columns=columns)


def assign_dataframe(reads: Any, targets: Any, k: int = 1, policy: str = "radius",
                     metric: str = "levenshtein", read_ids=None, target_names=None, *,
                     read_seq_col=None, read_id_col=None, target_seq_col=None, target_id_col=None):
    """Assign lists, Series or named dataframes without lossy string coercion."""
    from .dataframes import assign_dataframe as convert
    return convert(reads, targets, k, policy, metric, read_ids, target_names,
                   read_seq_col=read_seq_col, read_id_col=read_id_col,
                   target_seq_col=target_seq_col, target_id_col=target_id_col)


def _ensure_anndata() -> None:
    if not _HAS_ANNDATA:
        raise ImportError(
            "anndata is required; install with 'pip install \"dotmatch[anndata]\"' "
            "(pulls in pandas too). For full scanpy workflows, also install scanpy."
        )
    _ensure_pandas()


def counts_tsv_to_anndata(counts_path: str | Path, *, sample_cols=None,
                         var_cols=("target_id", "target_seq", "gene")):
    """Load raw integer counts as sparse samples-by-targets AnnData.

    sample_cols selects source column names in the requested order. Detailed
    DotMatch output uses only *_count_total columns by default. Missing values,
    fractional/negative counts and overflowing integers raise explicit errors.
    """
    from .dataframes import counts_tsv_to_anndata as convert
    return convert(counts_path, sample_cols=sample_cols, var_cols=var_cols)


def assignments_to_anndata(assignments: Any, *, cell_col="cell_barcode", feature_col="target_name",
                           status_col=None, count_unique_only=True, include_ambiguous_per_cell=False,
                           cell_names=None, feature_names=None):
    """Count uniquely assigned observations in a sparse cell-by-feature matrix.

    Cell labels must be explicit; they are never inferred from read IDs. Keep
    all observed cells, including those with no unique assignments. Optional
    cell_names/feature_names fix order and retain zero-count cells and features.
    No UMI deduplication or biological cell/perturbation calling is performed.
    """
    from .dataframes import assignments_to_anndata as convert
    return convert(assignments, cell_col=cell_col, feature_col=feature_col,
                   status_col=status_col, count_unique_only=count_unique_only,
                   include_ambiguous_per_cell=include_ambiguous_per_cell,
                   cell_names=cell_names, feature_names=feature_names)
