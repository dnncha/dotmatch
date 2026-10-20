"""Private ctypes ABI declarations and trusted native-library discovery.

Keep signatures here synchronized with include/qdalign.h. Matching semantics
and file parsing belong in their own modules.
"""

from __future__ import annotations

import ctypes
import os
import platform
from pathlib import Path


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


LIB = _load_lib()
