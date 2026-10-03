#!/usr/bin/env python3
"""Measure complete count commands on a public library and simulated FASTQs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from bench_hamming_seed_layout import Engine, ROOT, digest, simulated_reads


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-bin", required=True, type=Path)
    parser.add_argument("--candidate-bin", default=ROOT / "dotmatch", type=Path)
    parser.add_argument("--candidate-lib", default=ROOT / "libdotmatch.so", type=Path)
    parser.add_argument("--library", required=True, type=Path)
    parser.add_argument("--reads", default=30000, type=int)
    parser.add_argument("--repeats", default=5, type=int)
    parser.add_argument("--seed", default=20261002, type=int)
    parser.add_argument("--work-dir", default=ROOT / "benchmarks/work/hamming_seed_layout", type=Path)
    parser.add_argument("--out", default=ROOT / "benchmarks/raw/hamming_seed_layout_cli.csv", type=Path)
    args = parser.parse_args()
    if min(args.reads, args.repeats) < 1:
        parser.error("reads and repeats must be positive")
    with args.library.open(newline="") as handle:
        targets = [row["gRNA.sequence"].encode("ascii") for row in csv.DictReader(handle)]
    lengths = {len(t) for t in targets}
    if len(lengths) != 1:
        parser.error("this fixed-window benchmark requires a single target length")
    length = lengths.pop()
    reads = simulated_reads(targets, args.reads, args.seed + 19)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    fastq = args.work_dir / "simulated.fastq.gz"
    with fastq.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as handle:
            for i, read in enumerate(reads):
                handle.write(b"@simulated_" + str(i).encode() + b"\n" + read + b"\n+\n" + b"I" * len(read) + b"\n")
    rows = []
    provenance = {
        "command": [sys.executable, *sys.argv], "platform": platform.platform(),
        "seed": args.seed, "timing_scope": "single-thread count subprocess, including index build, gzip parsing and TSV/JSON output",
        "input_scope": "Yusa public guide library; simulated ACGT reads, five equal error groups; not public FASTQ evidence",
        "files": {str(p): digest(p) for p in [args.baseline_bin, args.candidate_bin, args.candidate_lib,
                  args.library, fastq, Path(__file__), ROOT / "scripts/bench_hamming_seed_layout.py"]},
    }
    engine = Engine(args.candidate_lib, targets, reads)
    try:
        for k in (2, 3):
            results = engine.assign(k)[1]
            for policy in ("best", "radius"):
                accepted = results[:, 4] == 1
                if policy == "radius":
                    accepted &= results[:, 3] == 1
                expected = np.bincount(results[accepted, 0], minlength=len(targets))
                baseline_hash = None
                baseline_outcomes = None
                for repeat in range(args.repeats):
                    order = ("baseline", "candidate") if repeat % 2 == 0 else ("candidate", "baseline")
                    for tool in order:
                        binary = args.baseline_bin if tool == "baseline" else args.candidate_bin
                        output = args.work_dir / f"counts-{k}-{policy}-{repeat}-{tool}.tsv"
                        summary = output.with_suffix(".json")
                        command = [str(binary.resolve()), "count", "--targets", str(args.library.resolve()),
                                   "--reads", str(fastq.resolve()), "--sample-label", "simulated",
                                   "--target-start", "0", "--target-length", str(length), "--k", str(k),
                                   "--metric", "hamming", "--ambiguity-policy", policy, "--threads", "1",
                                   "--format", "mageck", "--out", str(output.resolve()), "--summary", str(summary.resolve())]
                        start = time.perf_counter()
                        completed = subprocess.run(command, text=True, capture_output=True, check=True)
                        elapsed = time.perf_counter() - start
                        with output.open(newline="") as handle:
                            counts = np.array([int(row["simulated"]) for row in csv.DictReader(handle, delimiter="\t")])
                        errors = int(np.count_nonzero(counts != expected))
                        if errors:
                            raise RuntimeError(f"native batch/count disagreement: {k}/{policy}/{tool}: {errors}")
                        data = json.loads(summary.read_text())["samples"][0]
                        outcomes = {key: int(data[key]) for key in (
                            "total_reads", "assigned_unique", "assigned_exact", "assigned_corrected",
                            "ambiguous", "unmatched", "invalid", "candidates_considered", "candidates_verified"
                        )}
                        output_hash = digest(output)
                        if baseline_hash is None:
                            baseline_hash, baseline_outcomes = output_hash, outcomes
                        if output_hash != baseline_hash or outcomes != baseline_outcomes:
                            raise RuntimeError("baseline/candidate output or outcome disagreement")
                        rows.append({
                            "workload": "yusa_library_simulated_fastq", "tool": tool, "k": k, "policy": policy,
                            "repeat": repeat, "n_reads": len(reads), "n_targets": len(targets), "length": length,
                            "seconds": elapsed, "reads_per_second": len(reads) / elapsed,
                            "count_mismatches": errors, "counts_sha256": output_hash,
                            "outcomes_sha256": hashlib.sha256(json.dumps(outcomes, sort_keys=True).encode()).hexdigest(),
                            "exit_code": completed.returncode, "command": json.dumps(command), **outcomes,
                        })
                print(f"k={k}, {policy}: all {len(targets)} guide counts and read outcomes agree", flush=True)
    finally:
        engine.close()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    provenance["csv_sha256"] = digest(args.out)
    args.out.with_suffix(".json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
