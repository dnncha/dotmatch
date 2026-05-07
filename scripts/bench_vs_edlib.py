#!/usr/bin/env python3
"""Compare DotMatch against Edlib for global edit-distance workloads.

Usage:
  make shared
  python3 -m pip install edlib
  python3 scripts/bench_vs_edlib.py

The exact-distance benchmark compares qdaln_edit_distance against
edlib.align(..., mode="NW", task="distance", k=-1).

The threshold benchmark compares qdaln_edit_distance_leq against Edlib's
bounded global alignment, treating editDistance >= 0 as distance <= k and
editDistance == -1 as distance > k.
"""

from __future__ import annotations

import ctypes
import os
import platform
import random
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

try:
    import edlib  # type: ignore
except ImportError as exc:  # pragma: no cover
    raise SystemExit("Install edlib first: python3 -m pip install edlib") from exc

ROOT = Path(__file__).resolve().parents[1]
LIB_NAME = "libdotmatch.dylib" if platform.system() == "Darwin" else "libdotmatch.so"
LIB_PATH = ROOT / LIB_NAME

ExactFn = Callable[[str, str], int]
ThresholdFn = Callable[[str, str, int], int]


@dataclass(frozen=True)
class BenchResult:
    kind: str
    tool: str
    length: int
    k: int | None
    err: float
    calls: int
    seconds: float
    checksum: int

    @property
    def calls_per_sec(self) -> float:
        return self.calls / self.seconds

    @property
    def ns_per_call(self) -> float:
        return 1e9 / self.calls_per_sec

    def csv_row(self) -> str:
        k_value = "" if self.k is None else str(self.k)
        return (
            f"{self.kind},{self.tool},{self.length},{k_value},{self.err:.2f},"
            f"{self.calls},{self.seconds:.6f},{self.calls_per_sec:.1f},"
            f"{self.ns_per_call:.1f},{self.checksum}"
        )


@dataclass(frozen=True)
class SpeedupResult:
    kind: str
    length: int
    k: int | None
    err: float
    speedup: float
    qda_ns: float
    edlib_ns: float


def build_shared() -> None:
    if not LIB_PATH.exists():
        subprocess.run(["make", "shared"], cwd=ROOT, check=True)


def load_lib() -> ctypes.CDLL:
    build_shared()
    lib = ctypes.CDLL(str(LIB_PATH))
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
    return lib


def rand_seq(n: int) -> str:
    return "".join(random.choice("ACGT") for _ in range(n))


def mutate(seq: str, err: float) -> str:
    out = []
    for c in seq:
        if random.random() < err:
            choices = [x for x in "ACGT" if x != c]
            out.append(random.choice(choices))
        else:
            out.append(c)
    return "".join(out)


def make_pairs(length: int, err: float, calls: int) -> list[tuple[str, str]]:
    return [(a := rand_seq(length), mutate(a, err)) for _ in range(calls)]


def qda_distance(lib: ctypes.CDLL, a: str, b: str) -> int:
    aa = a.encode()
    bb = b.encode()
    return int(lib.qdaln_edit_distance(aa, len(aa), bb, len(bb)))


def qda_leq(lib: ctypes.CDLL, a: str, b: str, k: int) -> int:
    aa = a.encode()
    bb = b.encode()
    return int(lib.qdaln_edit_distance_leq(aa, len(aa), bb, len(bb), k))


def edlib_distance(a: str, b: str) -> int:
    return int(edlib.align(a, b, mode="NW", task="distance", k=-1)["editDistance"])


def edlib_leq(a: str, b: str, k: int) -> int:
    edit_distance = int(edlib.align(a, b, mode="NW", task="distance", k=k)["editDistance"])
    return 1 if edit_distance >= 0 else 0


def bench_exact(tool: str, fn: ExactFn, length: int, err: float, pairs: list[tuple[str, str]]) -> BenchResult:
    checksum = 0
    t0 = time.perf_counter()
    for a, b in pairs:
        checksum += fn(a, b)
    dt = time.perf_counter() - t0
    return BenchResult("exact", tool, length, None, err, len(pairs), dt, checksum)


def bench_threshold(
    tool: str,
    fn: ThresholdFn,
    length: int,
    k: int,
    err: float,
    pairs: list[tuple[str, str]],
) -> BenchResult:
    checksum = 0
    t0 = time.perf_counter()
    for a, b in pairs:
        checksum += fn(a, b, k)
    dt = time.perf_counter() - t0
    return BenchResult("threshold", tool, length, k, err, len(pairs), dt, checksum)


def verify_exact(lib: ctypes.CDLL, length: int, err: float, pairs: list[tuple[str, str]]) -> None:
    for a, b in pairs[:500]:
        qd = qda_distance(lib, a, b)
        ed = edlib_distance(a, b)
        if qd != ed:
            raise AssertionError(("exact", length, err, a, b, qd, ed))


def verify_threshold(lib: ctypes.CDLL, length: int, k: int, err: float, pairs: list[tuple[str, str]]) -> None:
    for a, b in pairs[:500]:
        qd = qda_leq(lib, a, b, k)
        ed = edlib_leq(a, b, k)
        if qd != ed:
            raise AssertionError(("threshold", length, k, err, a, b, qd, ed))


def print_pair(qda: BenchResult, edlib_result: BenchResult) -> SpeedupResult:
    print(qda.csv_row())
    print(edlib_result.csv_row())
    speedup = edlib_result.seconds / qda.seconds
    k_label = "" if qda.k is None else f" k={qda.k}"
    print(
        f"# {qda.kind} len={qda.length}{k_label} err={qda.err:.0%} "
        f"speedup_vs_edlib={speedup:.2f}x "
        f"dotmatch_ns={qda.ns_per_call:.1f} edlib_ns={edlib_result.ns_per_call:.1f}"
    )
    return SpeedupResult(qda.kind, qda.length, qda.k, qda.err, speedup, qda.ns_per_call, edlib_result.ns_per_call)


def print_interpretation(results: list[SpeedupResult]) -> None:
    threshold = [r for r in results if r.kind == "threshold"]
    wins = sorted((r for r in threshold if r.speedup > 1.0), key=lambda r: r.speedup, reverse=True)
    losses = sorted((r for r in threshold if r.speedup <= 1.0), key=lambda r: r.speedup)

    print()
    print("# interpretation")
    if wins:
        print("# DotMatch threshold wins")
        for r in wins[:8]:
            print(
                f"#   len={r.length} k={r.k} err={r.err:.0%} "
                f"speedup={r.speedup:.2f}x dotmatch_ns={r.qda_ns:.1f} edlib_ns={r.edlib_ns:.1f}"
            )
    else:
        print("# DotMatch did not beat Edlib in any threshold regime measured here")

    if losses:
        print("# edlib threshold wins")
        for r in losses[:8]:
            print(
                f"#   len={r.length} k={r.k} err={r.err:.0%} "
                f"speedup={r.speedup:.2f}x dotmatch_ns={r.qda_ns:.1f} edlib_ns={r.edlib_ns:.1f}"
            )
    else:
        print("# Edlib did not beat DotMatch in any threshold regime measured here")

    print("# broad state-of-the-art claims require broader benchmarks beyond these short global-distance regimes")


def main() -> None:
    random.seed(12345)
    lib = load_lib()

    lengths = [16, 24, 32, 64]
    thresholds = [0, 1, 2, 3]
    errors = [0.0, 0.01, 0.03, 0.10]
    calls = 20_000
    speedups: list[SpeedupResult] = []

    print(f"# python={platform.python_version()} platform={platform.platform()}")
    print("# external_baseline=edlib.align(mode='NW', task='distance')")
    print("kind,tool,len,k,err,calls,seconds,calls_per_sec,ns_per_call,checksum")

    print()
    print("# exact distance")
    for length in lengths:
        for err in errors:
            pairs = make_pairs(length, err, calls)
            verify_exact(lib, length, err, pairs)
            qda = bench_exact("dotmatch", lambda a, b: qda_distance(lib, a, b), length, err, pairs)
            ed = bench_exact("edlib", edlib_distance, length, err, pairs)
            speedups.append(print_pair(qda, ed))

    print()
    print("# threshold distance <= k")
    for length in lengths:
        for k in thresholds:
            for err in errors:
                pairs = make_pairs(length, err, calls)
                verify_threshold(lib, length, k, err, pairs)
                qda = bench_threshold("dotmatch", lambda a, b, kk: qda_leq(lib, a, b, kk), length, k, err, pairs)
                ed = bench_threshold("edlib", edlib_leq, length, k, err, pairs)
                speedups.append(print_pair(qda, ed))

    print_interpretation(speedups)


if __name__ == "__main__":
    os.chdir(ROOT)
    main()
