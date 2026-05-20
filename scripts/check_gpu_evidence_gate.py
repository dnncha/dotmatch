#!/usr/bin/env python3
"""Gate the experimental GPU acceleration evidence lane."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "gpu_acceleration.csv"
CRISPR_RAW = ROOT / "benchmarks" / "raw" / "gpu_crispr.csv"
REPORT = ROOT / "docs" / "benchmarks" / "gpu" / "README.md"


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def as_int(value: str | None) -> int:
    try:
        return int(value or 0)
    except ValueError:
        return 0


def case_key(row: dict[str, str]) -> tuple[str, str, str, str, str]:
    return (
        row.get("n_reads", ""),
        row.get("n_targets", ""),
        row.get("len", ""),
        row.get("k", ""),
        row.get("error_rate", ""),
    )


def real_case_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (
        row.get("packable_reads", ""),
        row.get("n_targets", ""),
        row.get("target_length", ""),
        row.get("k", ""),
    )


def row_gate(rows: list[dict[str, str]], failures: list[str], label: str = "gpu_acceleration.csv") -> None:
    if not rows:
        failures.append(f"{label} is empty")
        return
    unavailable = [row for row in rows if row.get("tool") == "dotmatch_gpu_metal" and row.get("status") == "unavailable"]
    gpu_rows = [row for row in rows if row.get("tool") == "dotmatch_gpu_metal" and row.get("status") == "ok"]
    cpu_rows = {
        case_key(row): row
        for row in rows
        if row.get("tool") == "dotmatch_cpu_index" and row.get("status") == "ok"
    }
    if not gpu_rows:
        if not unavailable:
            failures.append("missing successful GPU rows or explicit unavailable row")
        return
    for row in gpu_rows:
        key = case_key(row)
        if key not in cpu_rows:
            failures.append(f"missing CPU baseline for GPU case {key}")
            continue
        if as_int(row.get("mismatches")) != 0:
            failures.append(f"GPU mismatches for case {key}: {row.get('mismatches')}")
        if row.get("checksum") != cpu_rows[key].get("checksum"):
            failures.append(f"GPU checksum differs from CPU baseline for case {key}")
        if not row.get("device"):
            failures.append(f"GPU row missing device name for case {key}")


def real_row_gate(rows: list[dict[str, str]], failures: list[str]) -> None:
    if not rows:
        failures.append("gpu_crispr.csv is empty")
        return
    unavailable = [row for row in rows if row.get("tool") == "dotmatch_gpu_metal" and row.get("status") == "unavailable"]
    gpu_rows = [row for row in rows if row.get("tool") == "dotmatch_gpu_metal" and row.get("status") == "ok"]
    cpu_rows = {
        real_case_key(row): row
        for row in rows
        if row.get("tool") == "dotmatch_cpu_index" and row.get("status") == "ok"
    }
    if not gpu_rows:
        if not unavailable:
            failures.append("missing successful public CRISPR GPU rows or explicit unavailable row")
        return
    for row in gpu_rows:
        key = real_case_key(row)
        if key not in cpu_rows:
            failures.append(f"missing public CRISPR CPU baseline for GPU case {key}")
            continue
        if as_int(row.get("mismatches")) != 0:
            failures.append(f"public CRISPR GPU mismatches for case {key}: {row.get('mismatches')}")
        if as_int(row.get("count_delta")) != 0:
            failures.append(f"public CRISPR GPU count delta for case {key}: {row.get('count_delta')}")
        if row.get("checksum") != cpu_rows[key].get("checksum"):
            failures.append(f"public CRISPR GPU checksum differs from CPU baseline for case {key}")
        if as_int(row.get("packable_reads")) <= 0 or as_int(row.get("n_targets")) <= 0:
            failures.append(f"public CRISPR GPU row missing positive reads or targets for case {key}")
        if not row.get("device"):
            failures.append(f"public CRISPR GPU row missing device name for case {key}")


def report_gate(path: Path, failures: list[str]) -> None:
    if not path.exists():
        failures.append(f"missing GPU benchmark report: {path}")
        return
    text = path.read_text(encoding="utf-8")
    if "Experimental GPU Acceleration Benchmark" not in text:
        failures.append("GPU report missing expected title")
    if "not a production speed claim" not in text:
        failures.append("GPU report must keep the experimental scope boundary")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default=str(RAW))
    parser.add_argument("--crispr-csv", default=str(CRISPR_RAW))
    parser.add_argument("--report", default=str(REPORT))
    args = parser.parse_args(argv)

    failures: list[str] = []
    row_gate(read_rows(Path(args.csv)), failures)
    real_row_gate(read_rows(Path(args.crispr_csv)), failures)
    report_gate(Path(args.report), failures)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        print("GPU EVIDENCE: FAIL")
        return 1
    print("GPU EVIDENCE: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
