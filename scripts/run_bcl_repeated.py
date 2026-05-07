#!/usr/bin/env python3
"""Run repeated BCL demultiplexing benchmarks and combine raw rows."""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "bcl_demux_repeated.csv"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-folder", required=True)
    parser.add_argument("--sample-sheet", required=True)
    parser.add_argument("--workflow-name", default="real_bcl_repeated")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--gzip-level", type=int, default=1)
    parser.add_argument("--run-installed-competitors", action="store_true")
    parser.add_argument("--detect-competitors", action="store_true")
    parser.add_argument("--out", default=str(RAW))
    args = parser.parse_args()

    rows: list[dict[str, str]] = []
    fieldnames: list[str] | None = None
    with tempfile.TemporaryDirectory(prefix="dotmatch-bcl-repeated-") as tmp:
        tmpdir = Path(tmp)
        for repeat in range(1, args.repeats + 1):
            out = tmpdir / f"repeat_{repeat}.csv"
            cmd = [
                "python3", "scripts/bench_bcl_demux.py",
                "--run-folder", args.run_folder,
                "--sample-sheet", args.sample_sheet,
                "--workflow-name", args.workflow_name,
                "--threads", str(args.threads),
                "--gzip-level", str(args.gzip_level),
                "--out", str(out),
            ]
            if args.run_installed_competitors:
                cmd.append("--run-installed-competitors")
            if args.detect_competitors:
                cmd.append("--detect-competitors")
            subprocess.run(cmd, cwd=ROOT, check=True)
            with out.open() as fh:
                reader = csv.DictReader(fh)
                if fieldnames is None:
                    fieldnames = ["repeat", *reader.fieldnames] if reader.fieldnames else ["repeat"]
                for row in reader:
                    row = {"repeat": str(repeat), **row}
                    rows.append(row)

    if fieldnames is None:
        raise SystemExit("no benchmark rows produced")
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(out_path)


if __name__ == "__main__":
    main()
