from __future__ import annotations

import ctypes
import os
import platform
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

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


def assign(reads: Sequence[str | bytes], barcodes: Sequence[str | bytes], k: int = 1) -> list[MatchResult]:
    if k < 0:
        raise ValueError("k must be non-negative")
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

    return _results_to_python(results)


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

    def assign(self, reads: Sequence[str | bytes], k: int = 1) -> list[MatchResult]:
        results, _stats = self.assign_with_stats(reads, k=k)
        return results

    def assign_with_stats(self, reads: Sequence[str | bytes], k: int = 1) -> tuple[list[MatchResult], AssignmentStats]:
        if self._closed:
            raise ValueError("matcher is closed")
        if k < 0:
            raise ValueError("k must be non-negative")

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

        return _results_to_python(results), AssignmentStats(
            candidates_considered=int(stats.candidates_considered),
            candidates_verified=int(stats.candidates_verified),
        )
