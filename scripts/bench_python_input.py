#!/usr/bin/env python3
"""Compare Python input throughput with a pinned Git baseline, checking outputs.

Synthetic parsing measurements only; not whole-assay or native CLI timings.
Run from a source checkout with its shared library built (make shared).
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import io
import json
import platform
import statistics
import subprocess
import sys
import tempfile
import time
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))


def load_pair(name: str, baseline_ref: str):
    path = f"python/dotmatch/{name}.py"
    before = subprocess.run(["git", "show", f"{baseline_ref}:{path}"], cwd=ROOT,
                            check=True, capture_output=True).stdout
    after = (ROOT / path).read_bytes()
    module_name = f"dotmatch._input_benchmark_baseline_{name}"
    baseline = types.ModuleType(module_name)
    sys.modules[module_name] = baseline
    exec(compile(before, path, "exec"), baseline.__dict__)
    candidate = importlib.import_module(f"dotmatch.{name}")
    hashes = {"baseline": hashlib.sha256(before).hexdigest(),
              "candidate": hashlib.sha256(after).hexdigest()}
    return baseline, candidate, hashes


def measure(before, after, repeats: int):
    expected = before()  # Untimed warm-up and independent baseline output.
    if after() != expected:
        raise RuntimeError("candidate output differs from baseline")
    samples = {"baseline": [], "candidate": []}
    for repetition in range(repeats):
        calls = [("baseline", before), ("candidate", after)]
        if repetition % 2:
            calls.reverse()
        for name, call in calls:
            start = time.perf_counter()
            result = call()
            elapsed = time.perf_counter() - start
            if result != expected:
                raise RuntimeError(f"{name} output changed during measurement")
            samples[name].append(elapsed)
    medians = {name: statistics.median(values) for name, values in samples.items()}
    return {"seconds": samples, "median_seconds": medians,
            "speedup": medians["baseline"] / medians["candidate"],
            "output_sha256": expected}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-ref", required=True)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--reads", type=int, default=100000)
    parser.add_argument("--rows", type=int, default=10000)
    parser.add_argument("--samples", type=int, default=32)
    args = parser.parse_args()
    if min(args.repeats, args.reads, args.rows, args.samples) < 1:
        parser.error("workload sizes and repeats must be positive")
    baseline_ref = subprocess.run(["git", "rev-parse", "--verify", args.baseline_ref + "^{commit}"],
                                  cwd=ROOT, check=True, capture_output=True, text=True).stdout.strip()
    old_fastq, new_fastq, fastq_sources = load_pair("fastq_io", baseline_ref)
    old_counts, new_counts, count_sources = load_pair("count_io", baseline_ref)
    sequence = ("acgtn" * 30)
    quality = "I" * len(sequence)
    fastq = "".join(f"@r{i}\r\n{sequence}\r\n+r{i}\r\n{quality}\r\n" for i in range(args.reads))

    def parse_fastq(module):
        digest = hashlib.sha256()
        for record in module.iter_fastq_records(io.StringIO(fastq), "synthetic.fastq"):
            digest.update(f"{record.read_id}\0{record.seq}\0{record.qual}\n".encode())
        return digest.hexdigest()

    with tempfile.TemporaryDirectory(prefix="dotmatch-input-bench-") as directory:
        path = Path(directory) / "counts.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            handle.write("sgRNA\tGene\t" + "\t".join(f"s{i}" for i in range(args.samples)) + "\n")
            for row in range(args.rows):
                values = [str((row * 37 + column * 101) % 1000000) for column in range(args.samples)]
                if row % 1000 == 0:
                    values[0] = "9007199254740993"
                handle.write(f"g{row}\tG{row % 100}\t" + "\t".join(values) + "\n")

        def parse_counts(module):
            table = module.read_count_table(path)
            return hashlib.sha256(repr(table).encode()).hexdigest()

        report = {"baseline_commit": baseline_ref, "python": sys.version,
                  "platform": platform.platform(), "repeats": args.repeats,
                  "scope": "synthetic Python parsing including output hashing; excludes assignment, native CLI and gzip",
                  "fastq": {"records": args.reads, "read_length": len(sequence),
                            "input_sha256": hashlib.sha256(fastq.encode()).hexdigest(),
                            "source_sha256": fastq_sources,
                            **measure(lambda: parse_fastq(old_fastq), lambda: parse_fastq(new_fastq), args.repeats)},
                  "counts": {"rows": args.rows, "samples": args.samples,
                             "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                             "source_sha256": count_sources,
                             **measure(lambda: parse_counts(old_counts), lambda: parse_counts(new_counts), args.repeats)}}
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
