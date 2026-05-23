#!/usr/bin/env python3
"""Run the real public CRISPR Metal GPU benchmark lane."""

from __future__ import annotations

import argparse
import csv
import os
import platform
import subprocess
from io import StringIO
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "crispr_guides"
RAW = ROOT / "benchmarks" / "raw" / "gpu_crispr.csv"
FIELDS = [
    "tool",
    "backend",
    "path",
    "status",
    "workload",
    "total_reads",
    "packable_reads",
    "n_targets",
    "target_start",
    "target_length",
    "k",
    "input_seconds",
    "prep_seconds",
    "seconds",
    "total_seconds",
    "reads_per_sec",
    "total_reads_per_sec",
    "assigned_unique",
    "assigned_exact",
    "assigned_corrected",
    "ambiguous",
    "unmatched",
    "invalid_windows",
    "non_acgt_windows",
    "skipped_targets",
    "checksum",
    "mismatches",
    "count_delta",
    "candidate_count",
    "avg_candidates",
    "max_candidates",
    "device",
    "notes",
]


def unavailable_row(reason: str) -> dict[str, str]:
    row = {field: "" for field in FIELDS}
    row.update({
        "tool": "dotmatch_gpu_metal",
        "backend": "metal",
        "path": "unavailable",
        "status": "unavailable",
        "workload": "public_crispr_yusa_hamming",
        "total_reads": "0",
        "packable_reads": "0",
        "n_targets": "0",
        "target_start": "23",
        "target_length": "19",
        "k": "1",
        "input_seconds": "0.0",
        "prep_seconds": "0.0",
        "seconds": "0.0",
        "total_seconds": "0.0",
        "reads_per_sec": "0.0",
        "total_reads_per_sec": "0.0",
        "checksum": "0",
        "mismatches": "0",
        "count_delta": "0",
        "candidate_count": "0",
        "avg_candidates": "0.0",
        "max_candidates": "0",
        "device": platform.platform(),
        "notes": reason,
    })
    return row


def parse_rows(text: str) -> list[dict[str, str]]:
    rows = list(csv.DictReader(StringIO(text)))
    for row in rows:
        for field in FIELDS:
            row.setdefault(field, "")
    return rows


def ensure_public_data(records: str) -> None:
    data = EXAMPLE / "data"
    required = [
        data / "yusa_library.csv",
        data / "ERR376998.fastq.gz",
        data / "ERR376999.fastq.gz",
    ]
    if all(path.exists() for path in required):
        return
    cmd = [
        "python3",
        str(ROOT / "scripts" / "fetch_mageck_demo.py"),
        "--out",
        str(data),
        "--subsample",
        records,
    ]
    subprocess.run(cmd, cwd=ROOT, check=True)


def run_metal(args: argparse.Namespace) -> list[dict[str, str]]:
    metal_bin = Path(args.metal_bin)
    if platform.system() != "Darwin":
        return [unavailable_row("Metal benchmark is only available on Darwin hosts")]
    if not metal_bin.exists():
        return [unavailable_row(f"Metal CRISPR benchmark binary not found: {metal_bin}")]
    ensure_public_data(str(args.max_reads_per_sample))
    data = EXAMPLE / "data"
    cmd = [
        str(metal_bin),
        "--targets",
        str(data / "yusa_library.csv"),
        "--reads",
        str(data / "ERR376998.fastq.gz"),
        "--reads",
        str(data / "ERR376999.fastq.gz"),
        "--target-start",
        str(args.target_start),
        "--target-length",
        str(args.target_length),
        "--k",
        str(args.k),
        "--max-reads-per-sample",
        str(args.max_reads_per_sample),
    ]
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True)
    if result.returncode != 0:
        reason = result.stderr.strip() or result.stdout.strip() or f"exit {result.returncode}"
        return [unavailable_row(reason)]
    rows = parse_rows(result.stdout)
    return rows or [unavailable_row("Metal CRISPR benchmark produced no rows")]


def write_rows(rows: list[dict[str, str]], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in FIELDS})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metal-bin", default=str(ROOT / "build" / "bench_gpu_crispr_metal"))
    parser.add_argument("--out", default=str(RAW))
    parser.add_argument("--target-start", type=int, default=23)
    parser.add_argument("--target-length", type=int, default=19)
    parser.add_argument("--k", type=int, default=1)
    parser.add_argument("--max-reads-per-sample", type=int, default=int(os.environ.get("DOTMATCH_GPU_CRISPR_READS", "10000")))
    args = parser.parse_args()

    rows = run_metal(args)
    write_rows(rows, Path(args.out))
    print(Path(args.out))


if __name__ == "__main__":
    main()
