#!/usr/bin/env python3
"""Gate the amplicon/panel smoke evidence lane."""

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "amplicon_panel.csv"
PUBLIC_WORKFLOW = "public_nfcore_artic_v3_amplicon_primer"


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
        and row.get("workflow") == "synthetic_amplicon_panel_fixture"
        and row.get("status") == "smoke"
        and row.get("exit_code") == "0"
    ]
    if not passing:
        failures.append("missing successful synthetic amplicon/panel smoke row")
        return
    if any(as_int(row, "validation_mismatches") != 0 for row in passing):
        failures.append("amplicon/panel smoke row requires zero validation mismatches")
    validated = [row for row in passing if as_int(row, "validation_mismatches") == 0]
    if not any(as_int(row, "assigned_unique") > 0 for row in validated):
        failures.append("amplicon/panel smoke row must assign at least one target")
    if not any(as_int(row, "ambiguous_reads") > 0 for row in validated):
        failures.append("amplicon/panel smoke row must include an ambiguous-read diagnostic")
    if not any(as_int(row, "unmatched_reads") > 0 for row in validated):
        failures.append("amplicon/panel smoke row must include an unmatched-read diagnostic")
    if not any("dotmatch count" in row.get("command", "") for row in validated):
        failures.append("amplicon/panel smoke row must record the dotmatch count command")


def public_row_gate(rows: list[dict[str, str]], failures: list[str]) -> None:
    k0 = next(
        (
            row for row in rows
            if row.get("tool") == "dotmatch_count"
            and row.get("workflow") == PUBLIC_WORKFLOW
            and row.get("status") == "supported"
            and row.get("k") == "0"
            and row.get("exit_code") == "0"
        ),
        None,
    )
    k1 = next(
        (
            row for row in rows
            if row.get("tool") == "dotmatch_count"
            and row.get("workflow") == PUBLIC_WORKFLOW
            and row.get("status") == "supported"
            and row.get("k") == "1"
            and row.get("exit_code") == "0"
        ),
        None,
    )
    exact = next(
        (
            row for row in rows
            if row.get("tool") == "exact_prefix_hash"
            and row.get("workflow") == PUBLIC_WORKFLOW
            and row.get("status") == "supported"
            and row.get("k") == "0"
            and row.get("exit_code") == "0"
        ),
        None,
    )
    if k0 is None:
        failures.append("missing public nf-core ARTIC amplicon DotMatch k=0 row")
    if k1 is None:
        failures.append("missing public nf-core ARTIC amplicon DotMatch k=1 row")
    if exact is None:
        failures.append("missing public nf-core ARTIC amplicon exact-prefix baseline row")
    if k0 is None or k1 is None or exact is None:
        return
    public_rows = [k0, k1, exact]
    if any(as_int(row, "validation_mismatches") != 0 for row in public_rows):
        failures.append("public amplicon/panel rows require zero validation mismatches")
    if any(as_int(row, "assigned_unique") <= 0 for row in public_rows):
        failures.append("public amplicon/panel rows must assign at least one primer-start read")
    if as_int(k0, "assigned_unique") != as_int(exact, "assigned_unique"):
        failures.append("public amplicon/panel DotMatch k=0 row must agree with the exact-prefix baseline")
    if as_int(k0, "assigned_exact") != as_int(exact, "assigned_exact"):
        failures.append("public amplicon/panel DotMatch k=0 exact count must agree with the exact-prefix baseline")
    if as_int(k1, "assigned_unique") < as_int(k0, "assigned_unique"):
        failures.append("public amplicon/panel k=1 row must not assign fewer reads than k=0")
    if not k0.get("target_start") or not k0.get("target_length"):
        failures.append("public amplicon/panel rows must record the primer-start window")
    if any("dotmatch count" not in row.get("command", "") for row in [k0, k1]):
        failures.append("public amplicon/panel DotMatch rows must record dotmatch count commands")
    if "bench_amplicon_panel.py" not in exact.get("command", ""):
        failures.append("public amplicon/panel exact-prefix baseline must record the benchmark command")
    if any(not row.get("metadata") for row in public_rows):
        failures.append("public amplicon/panel rows must record metadata path")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default=str(RAW))
    parser.add_argument("--public", action="store_true", help="require public nf-core ARTIC amplicon rows")
    args = parser.parse_args(argv)

    failures: list[str] = []
    rows = read_rows(Path(args.csv))
    row_gate(rows, failures)
    if args.public:
        public_row_gate(rows, failures)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        print("AMPLICON/PANEL PUBLIC: FAIL" if args.public else "AMPLICON/PANEL SMOKE: FAIL")
        return 1
    print("AMPLICON/PANEL PUBLIC: PASS" if args.public else "AMPLICON/PANEL SMOKE: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
