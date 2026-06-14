#!/usr/bin/env python3
"""Gate the experimental GPU acceleration evidence lane."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "gpu_acceleration.csv"
CRISPR_RAW = ROOT / "benchmarks" / "raw" / "gpu_crispr.csv"
PRODUCTION_CRISPR_RAW = ROOT / "benchmarks" / "raw" / "crispr_cpu_metal.csv"
REPORT = ROOT / "docs" / "benchmarks" / "gpu" / "README.md"
PRODUCTION_CRISPR_REPORT = ROOT / "docs" / "benchmarks" / "gpu" / "production_crispr_cpu_metal.md"


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


def real_case_key(row: dict[str, str]) -> tuple[str, str, str, str, str, str, str]:
    return (
        row.get("workload", ""),
        row.get("total_reads", ""),
        row.get("packable_reads", ""),
        row.get("n_targets", ""),
        row.get("target_start", ""),
        row.get("target_length", ""),
        row.get("k", ""),
    )


def production_case_key(row: dict[str, str]) -> tuple[str, str]:
    return (
        row.get("dataset", ""),
        row.get("records_per_sample", ""),
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
        if as_int(row.get("skipped_targets")) != 0:
            failures.append(f"public CRISPR GPU row has skipped targets for case {key}: {row.get('skipped_targets')}")
        if not row.get("device"):
            failures.append(f"public CRISPR GPU row missing device name for case {key}")


def production_crispr_row_gate(rows: list[dict[str, str]], failures: list[str]) -> None:
    if not rows:
        failures.append("crispr_cpu_metal.csv is empty")
        return
    cpu_rows = {
        production_case_key(row): row
        for row in rows
        if row.get("backend") == "cpu" and row.get("status") == "ok"
    }
    metal_ok = [
        row for row in rows
        if row.get("backend") == "gpu-metal-experimental" and row.get("status") == "ok"
    ]
    metal_validate_failures = [
        row for row in rows
        if row.get("backend") == "gpu-metal-experimental" and row.get("metal_validate") == "1"
        and row.get("status") == "validation_failed"
    ]
    sanson_ineligible = [
        row for row in rows
        if row.get("dataset") == "sanson_brunello"
        and row.get("backend") == "gpu-metal-experimental"
        and row.get("status") == "ineligible"
        and row.get("offset_mode") == "multi"
    ]
    if not cpu_rows:
        failures.append("missing production CRISPR CPU authority rows")
    if not metal_ok and not metal_validate_failures and not sanson_ineligible:
        failures.append("missing production CRISPR Metal rows, validation failures, or ineligible rows")
        return
    for row in metal_ok:
        key = production_case_key(row)
        if row.get("metal_validate") == "0" and row.get("count_match_cpu") != "1" and not metal_validate_failures:
            failures.append(f"unvalidated non-identical production Metal row lacks validation-failure evidence for case {key}")
        cpu = cpu_rows.get(key)
        if cpu is None:
            failures.append(f"missing production CRISPR CPU baseline for Metal case {key}")
            continue
        if row.get("metal_validate") == "1" and (
            row.get("metal_validation") != "passed" or row.get("count_match_cpu") != "1"
        ):
            failures.append(f"validated production Metal row is not CPU-identical for case {key}")
        if row.get("backend_effective") != "gpu-metal-experimental":
            failures.append(f"production Metal row missing effective Metal backend for case {key}")
        if not row.get("count_engine"):
            failures.append(f"production Metal row missing count engine for case {key}")
        if not row.get("device"):
            failures.append(f"production Metal row missing device for case {key}")
    if not any(row.get("dataset") == "mageck_yusa" for row in metal_validate_failures):
        failures.append("missing Mageck/Yusa validation_failed row for experimental production Metal")
    if not sanson_ineligible:
        failures.append("missing Sanson/Brunello multi-offset ineligible Metal row")


def report_gate(path: Path, failures: list[str]) -> None:
    if not path.exists():
        failures.append(f"missing GPU benchmark report: {path}")
        return
    text = path.read_text(encoding="utf-8")
    if "Experimental GPU Acceleration Benchmark" not in text:
        failures.append("GPU report missing expected title")
    if "not a production speed claim" not in text:
        failures.append("GPU report must keep the experimental scope boundary")


def production_crispr_report_gate(path: Path, failures: list[str]) -> None:
    if not path.exists():
        failures.append(f"missing production CRISPR CPU-vs-Metal report: {path}")
        return
    text = path.read_text(encoding="utf-8")
    if not text.strip():
        failures.append(f"empty production CRISPR CPU-vs-Metal report: {path}")
        return
    required_fragments = [
        "CPU remains the assignment authority",
        "Metal is opt-in via `--backend gpu-metal-experimental`",
        "Treat any Metal speedup as advisory until `metal_validation=passed` and guide-by-guide counts match the CPU shadow run",
        "Do not use Metal for Sanson/Brunello-style multi-offset counting",
        "not count-identical to CPU",
        "`--metal-validate` fails on Yusa today",
        "Production `auto` staying on CPU is consistent with these measurements",
        "Sanson/Brunello remains CPU-only",
    ]
    for required in required_fragments:
        if required not in text:
            failures.append(f"production CRISPR CPU-vs-Metal report must retain evidence boundary: {required}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default=str(RAW))
    parser.add_argument("--crispr-csv", default=str(CRISPR_RAW))
    parser.add_argument("--production-crispr-csv", default=str(PRODUCTION_CRISPR_RAW))
    parser.add_argument("--report", default=str(REPORT))
    parser.add_argument("--production-crispr-report", default=str(PRODUCTION_CRISPR_REPORT))
    args = parser.parse_args(argv)

    failures: list[str] = []
    row_gate(read_rows(Path(args.csv)), failures)
    real_row_gate(read_rows(Path(args.crispr_csv)), failures)
    production_crispr_row_gate(read_rows(Path(args.production_crispr_csv)), failures)
    report_gate(Path(args.report), failures)
    production_crispr_report_gate(Path(args.production_crispr_report), failures)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        print("GPU EVIDENCE: FAIL")
        return 1
    print("GPU EVIDENCE: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
