#!/usr/bin/env python3
"""Gate the perturb-seq-style pair-assignment smoke evidence lane."""

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "perturb_seq.csv"
PUBLIC_WORKFLOW = "public_10x_crispr_guide_capture"


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
        if row.get("tool") == "dotmatch_pair_count"
        and row.get("workflow") == "synthetic_perturb_seq_fixture"
        and row.get("status") == "smoke"
        and row.get("exit_code") == "0"
    ]
    if not passing:
        failures.append("missing successful synthetic perturb-seq smoke row")
        return
    if any(as_int(row, "validation_mismatches") != 0 for row in passing):
        failures.append("perturb-seq smoke row requires zero validation mismatches")
    validated = [row for row in passing if as_int(row, "validation_mismatches") == 0]
    if not validated or any(as_int(row, "assigned_pairs") <= 0 for row in validated):
        failures.append("perturb-seq smoke row must assign at least one guide-feature pair")
    if not validated or any(as_int(row, "pair_ambiguous") <= 0 for row in validated):
        failures.append("perturb-seq smoke row must include an ambiguous pair diagnostic")
    if not validated or any(as_int(row, "left_unmatched") <= 0 for row in validated):
        failures.append("perturb-seq smoke row must include a left-unmatched diagnostic")
    if not validated or any(as_int(row, "right_unmatched") <= 0 for row in validated):
        failures.append("perturb-seq smoke row must include a right-unmatched diagnostic")
    if not validated or any(as_int(row, "invalid_reads") <= 0 for row in validated):
        failures.append("perturb-seq smoke row must include an invalid-read diagnostic")
    if not validated or any("dotmatch pair-count" not in row.get("command", "") for row in validated):
        failures.append("perturb-seq smoke row must record the dotmatch pair-count command")


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
            if row.get("tool") == "exact_slice_hash"
            and row.get("workflow") == PUBLIC_WORKFLOW
            and row.get("status") == "supported"
            and row.get("k") == "0"
            and row.get("exit_code") == "0"
        ),
        None,
    )
    if k0 is None:
        failures.append("missing public 10x CRISPR guide-capture DotMatch k=0 row")
    if k1 is None:
        failures.append("missing public 10x CRISPR guide-capture DotMatch k=1 row")
    if exact is None:
        failures.append("missing public 10x CRISPR guide-capture exact-slice baseline row")
    if k0 is None or k1 is None or exact is None:
        return
    public_rows = [k0, k1, exact]
    if any(as_int(row, "validation_mismatches") != 0 for row in public_rows):
        failures.append("public perturb-seq rows require zero validation mismatches")
    if any(as_int(row, "assigned_unique") <= 0 for row in public_rows):
        failures.append("public perturb-seq rows must assign at least one CRISPR guide read")
    if as_int(k0, "assigned_unique") != as_int(exact, "assigned_unique"):
        failures.append("public perturb-seq DotMatch k=0 row must agree with the exact-slice baseline")
    if as_int(k0, "assigned_exact") != as_int(exact, "assigned_exact"):
        failures.append("public perturb-seq DotMatch k=0 exact count must agree with the exact-slice baseline")
    if as_int(k1, "assigned_unique") < as_int(k0, "assigned_unique"):
        failures.append("public perturb-seq k=1 row must not assign fewer reads than k=0")
    if not k0.get("target_start") or not k0.get("target_length"):
        failures.append("public perturb-seq rows must record the observed fixed guide window")
    if any("dotmatch count" not in row.get("command", "") for row in [k0, k1]):
        failures.append("public perturb-seq DotMatch rows must record dotmatch count commands")
    if "bench_perturb_seq.py" not in exact.get("command", ""):
        failures.append("public perturb-seq exact-slice baseline must record the benchmark command")
    if any(not row.get("metadata") for row in public_rows):
        failures.append("public perturb-seq rows must record metadata path")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default=str(RAW))
    parser.add_argument("--public", action="store_true", help="require public 10x CRISPR guide-capture rows")
    args = parser.parse_args(argv)

    failures: list[str] = []
    rows = read_rows(Path(args.csv))
    row_gate(rows, failures)
    if args.public:
        public_row_gate(rows, failures)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        print("PERTURB-SEQ PUBLIC: FAIL" if args.public else "PERTURB-SEQ SMOKE: FAIL")
        return 1
    print("PERTURB-SEQ PUBLIC: PASS" if args.public else "PERTURB-SEQ SMOKE: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
