#!/usr/bin/env python3
"""Fail unless the raw-BCL benchmark evidence is strong enough for broad comparison."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "bcl_demux.csv"


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh))


def fail(reason: str) -> None:
    raise SystemExit(f"BCL comparison gate failed: {reason}")


def repeat_count(rows: list[dict[str, str]]) -> int:
    repeats = {r.get("repeat", "").strip() for r in rows if r.get("repeat", "").strip()}
    if repeats:
        return len(repeats)
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default=str(RAW))
    parser.add_argument("--require-cbcl", action="store_true", default=True)
    parser.add_argument("--min-speedup", type=float, default=10.0)
    parser.add_argument("--min-repeats", type=int, default=5)
    parser.add_argument("--allow-tiny-demo", action="store_true")
    args = parser.parse_args()

    rows = read_rows(Path(args.csv))
    if not rows:
        fail("no benchmark rows found")
    if not args.allow_tiny_demo and any("tiny" in r.get("workflow", "").lower() for r in rows):
        fail("tiny demo rows are present; run larger real BCL/CBCL benchmarks before broad comparison")
    if any("synthetic" in r.get("workflow", "") for r in rows):
        fail("synthetic rows are present; run real BCL/CBCL benchmarks before broad comparison")
    dotmatch = [r for r in rows if r.get("tool") == "dotmatch_bcl_demux" and r.get("exit_code") == "0"]
    if not dotmatch:
        fail("no successful DotMatch BCL row")
    dotmatch_repeats = repeat_count(dotmatch)
    if dotmatch_repeats < args.min_repeats:
        fail(f"repeated DotMatch BCL evidence required: {dotmatch_repeats} < {args.min_repeats} distinct successful repeats")
    if args.require_cbcl and not any(
        r.get("tool") == "dotmatch_bcl_demux" and "cbcl" in r.get("format", "").lower() and r.get("exit_code") == "0"
        for r in rows
    ):
        fail("no successful DotMatch CBCL row")
    competitors = [r for r in rows if r.get("tool") in {"bcl-convert", "bcl2fastq", "cuda-demux"} and r.get("exit_code") == "0"]
    if not competitors:
        fail("no successful competitor rows")
    validated = [r for r in competitors if r.get("validation_mismatches") == "0" and r.get("validation_exit_code") == "0"]
    if not validated:
        fail("no competitor row has zero-mismatch validation against DotMatch output")
    validated_repeats = repeat_count(validated)
    if validated_repeats < args.min_repeats:
        fail(f"repeated validated competitor evidence required: {validated_repeats} < {args.min_repeats} distinct successful repeats")
    dot_speed = max(float(r.get("clusters_per_sec") or 0.0) for r in dotmatch)
    comp_speed = max(float(r.get("clusters_per_sec") or 0.0) for r in validated)
    required = comp_speed * args.min_speedup
    if dot_speed < required:
        fail(f"DotMatch is below {args.min_speedup:.1f}x speedup over validated competitor rows ({dot_speed:.1f} < {required:.1f})")
    print("BCL comparison gate passed")


if __name__ == "__main__":
    main()
