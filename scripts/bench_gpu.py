#!/usr/bin/env python3
"""Run the experimental GPU acceleration benchmark and write raw CSV evidence."""

from __future__ import annotations

import argparse
import csv
import platform
import subprocess
from io import StringIO
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "gpu_acceleration.csv"
FIELDS = [
    "tool",
    "backend",
    "status",
    "workload",
    "n_reads",
    "n_targets",
    "len",
    "k",
    "error_rate",
    "prep_seconds",
    "seconds",
    "total_seconds",
    "reads_per_sec",
    "total_reads_per_sec",
    "pairs_per_sec",
    "checksum",
    "mismatches",
    "device",
    "notes",
]


def unavailable_row(reason: str) -> dict[str, str]:
    return {
        "tool": "dotmatch_gpu_metal",
        "backend": "metal",
        "status": "unavailable",
        "workload": "synthetic_hamming",
        "n_reads": "0",
        "n_targets": "0",
        "len": "0",
        "k": "1",
        "error_rate": "0.0",
        "prep_seconds": "0.0",
        "seconds": "0.0",
        "total_seconds": "0.0",
        "reads_per_sec": "0.0",
        "total_reads_per_sec": "0.0",
        "pairs_per_sec": "0.0",
        "checksum": "0",
        "mismatches": "0",
        "device": platform.platform(),
        "notes": reason,
    }


def parse_rows(text: str) -> list[dict[str, str]]:
    rows = list(csv.DictReader(StringIO(text)))
    for row in rows:
        for field in FIELDS:
            row.setdefault(field, "")
    return rows


def run_metal(metal_bin: Path) -> list[dict[str, str]]:
    if platform.system() != "Darwin":
        return [unavailable_row("Metal benchmark is only available on Darwin hosts")]
    if not metal_bin.exists():
        return [unavailable_row(f"Metal benchmark binary not found: {metal_bin}")]
    result = subprocess.run([str(metal_bin)], cwd=ROOT, text=True, capture_output=True)
    if result.returncode != 0:
        reason = result.stderr.strip() or result.stdout.strip() or f"exit {result.returncode}"
        return [unavailable_row(reason)]
    rows = parse_rows(result.stdout)
    return rows or [unavailable_row("Metal benchmark produced no rows")]


def write_rows(rows: list[dict[str, str]], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in FIELDS})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metal-bin", default=str(ROOT / "build" / "bench_gpu_metal"))
    parser.add_argument("--out", default=str(RAW))
    args = parser.parse_args()

    rows = run_metal(Path(args.metal_bin))
    write_rows(rows, Path(args.out))
    print(Path(args.out))


if __name__ == "__main__":
    main()
