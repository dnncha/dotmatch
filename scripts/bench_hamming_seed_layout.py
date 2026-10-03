#!/usr/bin/env python3
"""Paired native Hamming benchmark with an independent exhaustive byte oracle.

The public library workload uses simulated substitutions, not public FASTQs.
Both libraries must implement the same complete qdaln_match_result contract.
"""

from __future__ import annotations

import argparse
import csv
import ctypes as ct
import hashlib
import json
import platform
import random
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


class Match(ct.Structure):
    _fields_ = [(name, ct.c_int) for name in (
        "target_index", "best_distance", "second_best_distance", "match_count", "status"
    )]


class Stats(ct.Structure):
    _fields_ = [("candidates_considered", ct.c_size_t), ("candidates_verified", ct.c_size_t)]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sequence_digest(sequences: list[bytes]) -> str:
    h = hashlib.sha256()
    for seq in sequences:
        h.update(len(seq).to_bytes(8, "little"))
        h.update(seq)
    return h.hexdigest()


def result_array(results: ct.Array) -> np.ndarray:
    return np.frombuffer(results, dtype=np.int32).reshape(-1, 5).copy()


class Engine:
    def __init__(self, path: Path, targets: list[bytes], reads: list[bytes]):
        self.lib = ct.CDLL(str(path.resolve()))
        ptrs, lens = ct.POINTER(ct.c_char_p), ct.POINTER(ct.c_size_t)
        self.lib.qdaln_index_build.argtypes = [ptrs, lens, ct.c_size_t]
        self.lib.qdaln_index_build.restype = ct.c_void_p
        self.lib.qdaln_index_free.argtypes = [ct.c_void_p]
        self.lib.qdaln_index_assign_hamming_stats.argtypes = [
            ct.c_void_p, ptrs, lens, ct.c_size_t, ct.c_int, ct.POINTER(Match), ct.POINTER(Stats)
        ]
        self.lib.qdaln_index_assign_hamming_stats.restype = ct.c_int
        target_ptrs = (ct.c_char_p * len(targets))(*targets)
        target_lens = (ct.c_size_t * len(targets))(*(len(x) for x in targets))
        start = time.perf_counter()
        self.index = self.lib.qdaln_index_build(target_ptrs, target_lens, len(targets))
        self.build_seconds = time.perf_counter() - start
        if not self.index:
            raise RuntimeError("index construction failed")
        self.ptrs = (ct.c_char_p * len(reads))(*reads)
        self.lens = (ct.c_size_t * len(reads))(*(len(x) for x in reads))
        self.results = (Match * len(reads))()
        self.n_reads = len(reads)

    def assign(self, k: int) -> tuple[float, np.ndarray, Stats]:
        stats = Stats()
        start = time.perf_counter()
        rc = self.lib.qdaln_index_assign_hamming_stats(
            self.index, self.ptrs, self.lens, self.n_reads, k, self.results, ct.byref(stats)
        )
        elapsed = time.perf_counter() - start
        if rc != 0:
            raise RuntimeError(f"assignment failed: {rc}")
        return elapsed, result_array(self.results), stats

    def close(self) -> None:
        self.lib.qdaln_index_free(self.index)
        self.index = None


def oracle_distances(targets: list[bytes], reads: list[bytes], ids: list[int]) -> list[np.ndarray]:
    """Exhaustive literal byte comparison, with no seeds or packed DNA codes."""
    groups: dict[int, list[int]] = {}
    for i, target in enumerate(targets):
        groups.setdefault(len(target), []).append(i)
    matrices = {
        length: (np.array(indices), np.frombuffer(b"".join(targets[i] for i in indices), dtype=np.uint8)
                 .reshape(len(indices), length))
        for length, indices in groups.items() if length != 0
    }
    result = []
    for i in ids:
        read = reads[i]
        distances = np.full(len(targets), -1, dtype=np.int16)
        if len(read) == 0:
            distances[groups.get(0, [])] = 0
        elif len(read) in matrices:
            indices, matrix = matrices[len(read)]
            distances[indices] = np.count_nonzero(matrix != np.frombuffer(read, dtype=np.uint8), axis=1)
        result.append(distances)
    return result


def oracle_results(distances: list[np.ndarray], k: int) -> np.ndarray:
    results = []
    for d in distances:
        hits = np.flatnonzero((d >= 0) & (d <= k))
        if len(hits) == 0:
            results.append([-1, -1, -1, 0, 0])
            continue
        best = int(d[hits].min())
        ties = hits[d[hits] == best]
        other = d[hits][d[hits] > best]
        results.append([int(ties[0]), best, int(other.min()) if len(other) else -1,
                        len(hits), 2 if len(ties) > 1 else 1])
    return np.array(results, dtype=np.int32)


def simulated_reads(targets: list[bytes], n: int, seed: int, *, unknowns: bool = False) -> list[bytes]:
    rng = random.Random(seed)
    reads = []
    source_targets = targets if unknowns else [t for t in targets if set(t) <= set(b"ACGT")]
    for i in range(n):
        target = source_targets[rng.randrange(len(source_targets))]
        seq = bytearray(target)
        mode = i % (6 if unknowns else 5)
        if mode == 4:
            seq = bytearray(rng.choices(b"ACGT", k=len(target)))
        else:
            for pos in rng.sample(range(len(seq)), min(mode, 3, len(seq))):
                seq[pos] = rng.choice([b for b in b"ACGT" if b != seq[pos]])
        if unknowns and mode == 5 and seq:
            seq[rng.randrange(len(seq))] = ord("N")
        reads.append(bytes(seq))
    return reads


def workloads(args: argparse.Namespace):
    rng = random.Random(args.seed)
    for n, length in [(2048, 20), (65536, 20), (65536, 32), (1024, 8)]:
        dense = length == 8
        targets = [bytes(rng.choices(b"ACGT", k=length)) for _ in range(n)]
        if dense:
            targets = [b"AA" + t[2:] for t in targets]
        # Duplicate rows and literal targets exercise full metadata semantics.
        targets[-2] = targets[0]
        targets[-1] = targets[0][:-1] + b"N"
        label = "dense_shared_seed_1024_8" if dense else f"random_{n}_{length}"
        yield label, targets, simulated_reads(targets, args.reads, args.seed + n + length)
        if dense:
            yield "dense_literal_reads_1024_8", targets, simulated_reads(
                targets, min(args.reads, 3000), args.seed + n + length, unknowns=True
            )
    if args.public_library:
        with args.public_library.open(newline="") as fh:
            targets = [row["gRNA.sequence"].encode("ascii") for row in csv.DictReader(fh)]
        yield "yusa_library_simulated_reads", targets, simulated_reads(targets, args.reads, args.seed + 19)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-lib", required=True, type=Path)
    parser.add_argument("--candidate-lib", default=ROOT / "libdotmatch.so", type=Path)
    parser.add_argument("--baseline-source", required=True, type=Path)
    parser.add_argument("--candidate-source", default=ROOT / "src/qdalign.c", type=Path)
    parser.add_argument("--baseline-ref", required=True)
    parser.add_argument("--public-library", type=Path)
    parser.add_argument("--reads", type=int, default=30000)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--oracle-reads", type=int, default=256)
    parser.add_argument("--seed", type=int, default=20261002)
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks/raw/hamming_seed_layout.csv")
    args = parser.parse_args()
    if min(args.reads, args.repeats, args.oracle_reads) < 1:
        parser.error("reads, repeats and oracle-reads must be positive")
    rows = []
    provenance = {
        "command": [sys.executable, *sys.argv],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(), "processor": platform.processor(),
        "python": platform.python_version(), "numpy": np.__version__,
        "compiler": subprocess.check_output(["cc", "--version"], text=True).splitlines()[0],
        "baseline_ref": args.baseline_ref, "seed": args.seed,
        "semantics": "literal-byte, same-length Hamming k=0..3; complete best-match metadata",
        "read_generation": "ACGT reads in equal groups: exact, 1/2/3 substitutions, random decoy; dense_literal_reads adds a sixth group with 3 substitutions plus N",
        "timing_scope": "single-thread native batch assignment; includes stats; excludes index build and Python marshalling",
        "oracle": "exhaustive NumPy byte inequality/count_nonzero; all target rows, sampled reads",
        "files": {str(p): digest(p) for p in [args.baseline_source, args.candidate_source,
                  args.baseline_lib, args.candidate_lib, Path(__file__)]},
        "workloads": [],
    }
    if args.public_library:
        provenance["files"][str(args.public_library)] = digest(args.public_library)
    for name, targets, reads in workloads(args):
        engines = {"baseline": Engine(args.baseline_lib, targets, reads),
                   "candidate": Engine(args.candidate_lib, targets, reads)}
        try:
            # Exhaust every read for small dense libraries; sample large
            # libraries where a full scan would dominate the benchmark run.
            ids = (list(range(len(reads))) if len(targets) <= 1024 else
                   sorted(random.Random(args.seed).sample(range(len(reads)), min(args.oracle_reads, len(reads)))))
            distances = oracle_distances(targets, reads, ids)
            oracle_hash = hashlib.sha256(b"".join(d.tobytes() for d in distances)).hexdigest()
            provenance["workloads"].append({
                "name": name, "n_targets": len(targets), "n_reads": len(reads),
                "target_sha256": sequence_digest(targets), "read_sha256": sequence_digest(reads),
                "oracle_read_ids": ids, "oracle_distances_sha256": oracle_hash,
            })
            for k in range(4):
                expected = oracle_results(distances, k)
                # Warm both engines, validate all fields against each other and
                # sampled exhaustive results, then alternate paired run order.
                reference = engines["baseline"].assign(k)[1]
                candidate = engines["candidate"].assign(k)[1]
                oracle_errors = int(np.count_nonzero(np.any(candidate[ids] != expected, axis=1)))
                baseline_errors = int(np.count_nonzero(np.any(reference[ids] != expected, axis=1)))
                changed = np.any(candidate != reference, axis=1)
                disagreement = int(np.count_nonzero(changed))
                # The old unknown-read path used a Levenshtein verifier. Only
                # those reads may change; every candidate result must be right.
                changed_acgt = any(set(reads[i]) <= set(b"ACGT") for i in np.flatnonzero(changed))
                if oracle_errors or changed_acgt:
                    raise RuntimeError(f"{name} k={k}: oracle={oracle_errors}/{baseline_errors}, baseline={disagreement}")
                for repeat in range(args.repeats):
                    order = ("baseline", "candidate") if repeat % 2 == 0 else ("candidate", "baseline")
                    for tool in order:
                        seconds, actual, stats = engines[tool].assign(k)
                        if not np.array_equal(actual, candidate if tool == "candidate" else reference):
                            raise RuntimeError(f"nondeterministic result: {name}/{k}/{tool}")
                        rows.append({
                            "workload": name, "tool": tool, "k": k, "repeat": repeat,
                            "n_reads": len(reads), "n_targets": len(targets),
                            "length": len(targets[0]), "seconds": seconds,
                            "reads_per_second": len(reads) / seconds,
                            "index_build_seconds": engines[tool].build_seconds,
                            "candidates_considered": stats.candidates_considered,
                            "candidates_verified": stats.candidates_verified,
                            "oracle_reads": len(ids), "oracle_mismatches": oracle_errors if tool == "candidate" else baseline_errors,
                            "baseline_mismatches": disagreement,
                            "baseline_status_changes": int(np.count_nonzero(candidate[:, 4] != reference[:, 4])),
                            "baseline_unique_assignment_changes": int(np.count_nonzero(
                                ((candidate[:, 4] == 1) | (reference[:, 4] == 1)) &
                                ((candidate[:, 4] != reference[:, 4]) | (candidate[:, 0] != reference[:, 0])))),
                            "results_sha256": hashlib.sha256(actual.tobytes()).hexdigest(),
                            "targets_sha256": sequence_digest(targets), "reads_sha256": sequence_digest(reads),
                        })
                print(f"{name} k={k}: {len(reads)} results, {disagreement} corrected; "
                      f"{len(ids)} exhaustive checks, baseline errors={baseline_errors}, candidate errors={oracle_errors}", flush=True)
        finally:
            for engine in engines.values():
                engine.close()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    provenance["csv_sha256"] = digest(args.out)
    targets, reads = [b"AATCGGAT", b"AAACGGAT"], [b"ANATGGAT", b"ARATGGAT", b"AaATGGAT"]
    expected = oracle_results(oracle_distances(targets, reads, [0, 1, 2]), 3)
    regression = {"targets": [t.decode() for t in targets], "reads": [r.decode() for r in reads],
                  "k": 3, "fields": [name for name, _ in Match._fields_], "oracle": expected.tolist()}
    for tool, path in [("baseline", args.baseline_lib), ("candidate", args.candidate_lib)]:
        engine = Engine(path, targets, reads)
        try:
            result = engine.assign(3)[1]
            regression[tool] = result.tolist()
            if tool == "candidate" and not np.array_equal(result, expected):
                raise RuntimeError("unknown-read regression failed")
        finally:
            engine.close()
    provenance["unknown_read_regression"] = regression
    args.out.with_suffix(".json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
