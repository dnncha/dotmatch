#!/usr/bin/env python3
"""Gate oligo/adapter fixed-window smoke and public assignment evidence."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "oligo_adapter.csv"
WORKFLOW = "synthetic_oligo_adapter_fixture"
PUBLIC_WORKFLOW = "public_fast_adapter_truseq_r1"


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
        and row.get("workflow") == WORKFLOW
        and row.get("status") == "smoke"
        and row.get("exit_code") == "0"
    ]
    if not passing:
        failures.append("missing successful synthetic oligo/adapter smoke row")
        return
    if any(as_int(row, "validation_mismatches") != 0 for row in passing):
        failures.append("oligo/adapter smoke row requires zero validation mismatches")
    validated = [row for row in passing if as_int(row, "validation_mismatches") == 0]
    if not validated or any(as_int(row, "assigned_unique") <= 0 for row in validated):
        failures.append("oligo/adapter smoke row must assign at least one oligo")
    if not validated or any(as_int(row, "corrected_reads") <= 0 for row in validated):
        failures.append("oligo/adapter smoke row must include a corrected one-substitution diagnostic")
    if not validated or any(as_int(row, "ambiguous_reads") <= 0 for row in validated):
        failures.append("oligo/adapter smoke row must include an ambiguous-read diagnostic")
    if not validated or any(as_int(row, "unmatched_reads") <= 0 for row in validated):
        failures.append("oligo/adapter smoke row must include an unmatched-read diagnostic")
    if not validated or any(row.get("target_start") != "0" or as_int(row, "target_length") <= 0 for row in validated):
        failures.append("oligo/adapter smoke row must record a fixed target window")
    if not validated or any(row.get("k") != "1" or row.get("metric") != "hamming" for row in validated):
        failures.append("oligo/adapter smoke row must use the hamming k=1 diagnostic lane")
    if not validated or any("dotmatch count" not in row.get("command", "") for row in validated):
        failures.append("oligo/adapter smoke row must record the dotmatch count command")


def public_row_gate(rows: list[dict[str, str]], failures: list[str]) -> None:
    public_rows = [
        row for row in rows
        if row.get("workflow") == PUBLIC_WORKFLOW
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
        failures.append("missing public adapter-prefix DotMatch k=0 row")
    if not dotmatch_k1:
        failures.append("missing public adapter-prefix DotMatch k=1 row")
    if not exact:
        failures.append("missing public adapter-prefix exact-slice baseline row")
    if not dotmatch_k0 or not dotmatch_k1 or not exact:
        return

    k0 = dotmatch_k0[0]
    k1 = dotmatch_k1[0]
    baseline = exact[0]
    if any(as_int(row, "validation_mismatches") != 0 for row in [k0, k1, baseline]):
        failures.append("public adapter-prefix rows require zero validation mismatches")
    if as_int(k0, "assigned_unique") <= 0 or as_int(baseline, "assigned_unique") <= 0:
        failures.append("public adapter-prefix rows must assign at least one read")
    if as_int(k0, "assigned_unique") != as_int(baseline, "assigned_unique"):
        failures.append("public adapter-prefix DotMatch k=0 row must agree with the exact-slice baseline")
    if as_int(k0, "assigned_exact") != as_int(baseline, "assigned_exact"):
        failures.append("public adapter-prefix DotMatch k=0 exact count must agree with the exact-slice baseline")
    if as_int(k1, "assigned_unique") < as_int(k0, "assigned_unique"):
        failures.append("public adapter-prefix k=1 row must not assign fewer reads than k=0")
    if any(row.get("target_start") != "229" or row.get("target_length") != "20" for row in [k0, k1, baseline]):
        failures.append("public adapter-prefix rows must use the fixed R1 window start 229 length 20")
    if "dotmatch count" not in k0.get("command", "") or "dotmatch count" not in k1.get("command", ""):
        failures.append("public adapter-prefix DotMatch rows must record dotmatch count commands")
    if "bench_oligo_adapter.py" not in baseline.get("command", ""):
        failures.append("public adapter-prefix exact-slice baseline must record the benchmark command")
    if not k0.get("metadata"):
        failures.append("public adapter-prefix rows must record metadata path")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default=str(RAW))
    parser.add_argument("--public", action="store_true")
    args = parser.parse_args(argv)

    failures: list[str] = []
    rows = read_rows(Path(args.csv))
    row_gate(rows, failures)
    if args.public:
        public_row_gate(rows, failures)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        print("OLIGO/ADAPTER GATE: FAIL")
        return 1
    print("OLIGO/ADAPTER GATE: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
