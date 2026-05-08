#!/usr/bin/env python3
"""Gate feature-barcode smoke and public assignment evidence lanes."""

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "feature_barcode.csv"


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def as_int(row: dict[str, str], key: str) -> int:
    try:
        return int(row.get(key) or 0)
    except ValueError:
        return 0


def row_gate(rows: list[dict[str, str]], failures: list[str]) -> None:
    passing = [
        row for row in rows
        if row.get("tool") == "dotmatch_count"
        and row.get("workflow") == "synthetic_feature_barcode_fixture"
        and row.get("status") == "smoke"
        and row.get("exit_code") == "0"
    ]
    if not passing:
        failures.append("missing successful synthetic feature-barcode smoke row")
        return
    if any(as_int(row, "validation_mismatches") != 0 for row in passing):
        failures.append("feature-barcode smoke row requires zero validation mismatches")
    validated = [row for row in passing if as_int(row, "validation_mismatches") == 0]
    if not validated or any(as_int(row, "assigned_unique") <= 0 for row in validated):
        failures.append("feature-barcode smoke row must assign at least one feature barcode")
    if not validated or any(as_int(row, "ambiguous_reads") <= 0 for row in validated):
        failures.append("feature-barcode smoke row must include an ambiguous-read diagnostic")
    if not validated or any(as_int(row, "unmatched_reads") <= 0 for row in validated):
        failures.append("feature-barcode smoke row must include an unmatched-read diagnostic")
    if not validated or any("dotmatch count" not in row.get("command", "") for row in validated):
        failures.append("feature-barcode smoke row must record the dotmatch count command")


def public_row_gate(rows: list[dict[str, str]], failures: list[str]) -> None:
    public_rows = [
        row for row in rows
        if row.get("workflow") == "public_10x_totalseq_b_feature_barcode"
        and row.get("status") == "supported"
        and row.get("exit_code") == "0"
    ]
    dotmatch_k0 = [
        row for row in public_rows
        if row.get("tool") == "dotmatch_count"
        and row.get("k") == "0"
        and row.get("metric") == "hamming"
    ]
    dotmatch_k1 = [
        row for row in public_rows
        if row.get("tool") == "dotmatch_count"
        and row.get("k") == "1"
        and row.get("metric") == "hamming"
    ]
    exact = [
        row for row in public_rows
        if row.get("tool") == "exact_slice_hash"
        and row.get("k") == "0"
        and row.get("metric") == "exact"
    ]
    if not dotmatch_k0:
        failures.append("missing public 10x feature-barcode DotMatch k=0 row")
    if not dotmatch_k1:
        failures.append("missing public 10x feature-barcode DotMatch k=1 row")
    if not exact:
        failures.append("missing public 10x feature-barcode exact-slice baseline row")
    if not dotmatch_k0 or not dotmatch_k1 or not exact:
        return

    k0 = dotmatch_k0[0]
    k1 = dotmatch_k1[0]
    baseline = exact[0]
    if any(as_int(row, "validation_mismatches") != 0 for row in [k0, k1, baseline]):
        failures.append("public feature-barcode rows require zero validation mismatches")
    if as_int(k0, "assigned_unique") <= 0 or as_int(baseline, "assigned_unique") <= 0:
        failures.append("public feature-barcode rows must assign at least one read")
    if as_int(k0, "assigned_unique") != as_int(baseline, "assigned_unique"):
        failures.append("public feature-barcode DotMatch k=0 row must agree with the exact-slice baseline")
    if as_int(k0, "assigned_exact") != as_int(baseline, "assigned_exact"):
        failures.append("public feature-barcode DotMatch k=0 exact count must agree with the exact-slice baseline")
    if as_int(k1, "assigned_unique") < as_int(k0, "assigned_unique"):
        failures.append("public feature-barcode k=1 row must not assign fewer reads than k=0")
    if k0.get("target_start") != "10" or k0.get("target_length") != "15":
        failures.append("public feature-barcode rows must use the 10x fixed window start 10 length 15")
    if "dotmatch count" not in k0.get("command", "") or "dotmatch count" not in k1.get("command", ""):
        failures.append("public feature-barcode DotMatch rows must record dotmatch count commands")
    if "bench_feature_barcode.py" not in baseline.get("command", ""):
        failures.append("public feature-barcode exact-slice baseline must record the benchmark command")
    if not k0.get("metadata"):
        failures.append("public feature-barcode rows must record metadata path")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default=str(RAW))
    args = parser.parse_args(argv)

    failures: list[str] = []
    rows = read_rows(Path(args.csv))
    row_gate(rows, failures)
    public_row_gate(rows, failures)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        print("FEATURE BARCODE GATE: FAIL")
        return 1
    print("FEATURE BARCODE GATE: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
