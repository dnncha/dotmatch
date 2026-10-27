#!/usr/bin/env python3
"""Deterministic Python-API throughput harness.

This measures the whole binding path (argument marshalling, native kernel, and
Python result objects), which is what a DotMatch Python user actually pays.
"""

from __future__ import annotations

import argparse
import json
import random
import statistics
import tempfile
import time
from typing import Callable, Sequence

import dotmatch as dm


def _sequences(rng: random.Random, count: int, length: int) -> list[str]:
    return ["".join(rng.choices("ACGT", k=length)) for _ in range(count)]


def _reads(rng: random.Random, targets: Sequence[str], count: int, edit_rate: float) -> list[str]:
    reads: list[str] = []
    for _ in range(count):
        source = list(rng.choice(targets))
        if rng.random() < edit_rate and source:
            index = rng.randrange(len(source))
            source[index] = rng.choice([base for base in "ACGT" if base != source[index]])
        reads.append("".join(source))
    return reads


def _time_once(operation: Callable[[], object]) -> float:
    start = time.perf_counter()
    operation()
    return time.perf_counter() - start


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=int, default=4096)
    parser.add_argument("--reads", type=int, default=40000)
    parser.add_argument("--length", type=int, default=20)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument(
        "--workload", choices=("hamming", "levenshtein", "exact", "stream"), default="hamming"
    )
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    targets = _sequences(rng, args.targets, args.length)
    reads = _reads(rng, targets, args.reads, 0.7)
    matcher = dm.Matcher(targets)
    temporary = tempfile.TemporaryDirectory(prefix="dotmatch-bench-")
    fastq = f"{temporary.name}/reads.fastq"
    with open(fastq, "w", encoding="ascii") as handle:
        for index, read in enumerate(reads):
            handle.write(f"@r{index}\n{read}\n+\n{'I' * len(read)}\n")
    try:
        def operation() -> object:
            if args.workload == "hamming":
                return matcher.assign_hamming(reads, k=1)
            if args.workload == "levenshtein":
                return matcher.assign(reads, k=1)
            if args.workload == "exact":
                return matcher.assign_exact(reads)
            return sum(
                1
                for _assignment in dm.stream_assign(fastq, targets, metric="hamming", k=1)
            )
        operation()
        samples = sorted(_time_once(operation) for _ in range(args.repeats))
    finally:
        matcher.close()
        temporary.cleanup()
    median = statistics.median(samples)
    payload = {
        "workload": args.workload,
        "targets": args.targets,
        "reads": args.reads,
        "length": args.length,
        "repeats": args.repeats,
        "median_seconds": median,
        "min_seconds": samples[0],
        "max_seconds": samples[-1],
        "median_reads_per_second": args.reads / median if median else None,
        "samples_seconds": samples,
    }
    if args.json:
        print(json.dumps(payload, sort_keys=True))
    else:
        print(
            f"{args.workload}: {payload['median_reads_per_second']:.0f} reads/s "
            f"(median {median:.4f}s, {args.reads} reads, {args.targets} targets)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
