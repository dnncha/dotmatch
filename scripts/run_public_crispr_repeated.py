#!/usr/bin/env python3
"""Run repeated real-data public CRISPR benchmark rows for claim-grade summaries."""

from __future__ import annotations

import argparse
import csv
import os
import platform
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "benchmarks" / "raw" / "public_crispr_repeated.csv"


def parse_sizes(text: str) -> list[int]:
    out: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        value = int(part)
        if value <= 0:
            raise ValueError("read sizes must be positive")
        out.append(value)
    if not out:
        raise ValueError("at least one read size is required")
    return out


def metadata() -> dict[str, str]:
    return {
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python": platform.python_version(),
    }


def run_once(args: argparse.Namespace, records_per_sample: int | None, repeat: int, tmp_out: Path) -> list[dict[str, str]]:
    cmd = ["python3", str(ROOT / "scripts" / "run_public_crispr_benchmark.py"), "--out", str(tmp_out)]
    env = os.environ.copy()
    if records_per_sample is not None:
        cmd.append("--small")
        env["DOTMATCH_PUBLIC_SUBSAMPLE"] = str(records_per_sample)
    if args.run_mageck:
        cmd.append("--run-mageck")
    if args.run_guide_counter:
        cmd.append("--run-guide-counter")
    if args.run_cutadapt:
        cmd.append("--run-cutadapt")
    if args.run_bowtie2:
        cmd.append("--run-bowtie2")
    if args.dotmatch_threads > 1:
        cmd.extend(["--dotmatch-threads", str(args.dotmatch_threads)])
    subprocess.run(cmd, cwd=ROOT, env=env, check=True)
    with tmp_out.open() as fh:
        rows = list(csv.DictReader(fh))
    extra = metadata()
    for row in rows:
        row["repeat"] = str(repeat)
        row["requested_records_per_sample"] = "full" if records_per_sample is None else str(records_per_sample)
        row["benchmark_command"] = " ".join(cmd)
        row.update(extra)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--read-sizes", default=os.environ.get("DOTMATCH_PUBLIC_READ_SIZES", "10000,100000"),
            help="comma-separated FASTQ records per sample for subsampled runs")
    parser.add_argument("--repeats", type=int, default=int(os.environ.get("DOTMATCH_PUBLIC_REPEATS", "5")))
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--run-mageck", action="store_true")
    parser.add_argument("--run-guide-counter", action="store_true")
    parser.add_argument("--run-cutadapt", action="store_true")
    parser.add_argument("--run-bowtie2", action="store_true")
    parser.add_argument("--dotmatch-threads", type=int, default=int(os.environ.get("DOTMATCH_COUNT_THREADS", "1")))
    parser.add_argument("--full", action="store_true", help="also run the full public FASTQ inputs")
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")

    sizes = parse_sizes(args.read_sizes)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, str]] = []

    with tempfile.TemporaryDirectory(prefix="dotmatch-public-crispr-") as tmp_dir:
        tmp_root = Path(tmp_dir)
        for records_per_sample in sizes:
            for repeat in range(1, args.repeats + 1):
                tmp_out = tmp_root / f"public_crispr_{records_per_sample}_{repeat}.csv"
                all_rows.extend(run_once(args, records_per_sample, repeat, tmp_out))
        if args.full:
            for repeat in range(1, args.repeats + 1):
                tmp_out = tmp_root / f"public_crispr_full_{repeat}.csv"
                all_rows.extend(run_once(args, None, repeat, tmp_out))

    if not all_rows:
        raise SystemExit("no benchmark rows were produced")
    fieldnames: list[str] = []
    for row in all_rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)
    print(out_path)


if __name__ == "__main__":
    main()
