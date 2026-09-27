#!/usr/bin/env python3
"""Check a Yusa MAGeCK-format count table against MAGeCK's published tutorial rows."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFERENCE = ROOT / "benchmarks" / "reference" / "mageck_yusa_tutorial_rows.csv"


def read_unique_rows(path: Path, key: str) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        delimiter = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None or key not in reader.fieldnames:
            raise ValueError(f"{path} must contain a {key!r} column")
        rows: dict[str, dict[str, str]] = {}
        for row in reader:
            value = row.get(key, "")
            if not value:
                raise ValueError(f"{path} contains an empty {key!r}")
            if value in rows:
                raise ValueError(f"{path} contains duplicate {key!r}: {value}")
            rows[value] = row
    return rows


def compare_rows(reference_path: Path, counts_path: Path) -> list[dict[str, str]]:
    reference = read_unique_rows(reference_path, "guide_id")
    observed = read_unique_rows(counts_path, "sgRNA")
    details: list[dict[str, str]] = []
    for guide_id, expected in reference.items():
        actual = observed.get(guide_id)
        row = {
            "guide_id": guide_id,
            "expected_gene": expected.get("gene", ""),
            "observed_gene": "" if actual is None else actual.get("Gene", ""),
            "expected_plasmid": expected.get("plasmid", ""),
            "observed_plasmid": "" if actual is None else actual.get("plasmid", ""),
            "expected_ESC1": expected.get("ESC1", ""),
            "observed_ESC1": "" if actual is None else actual.get("ESC1", ""),
        }
        row["status"] = "match" if actual is not None and (
            row["expected_gene"], row["expected_plasmid"], row["expected_ESC1"]
        ) == (
            row["observed_gene"], row["observed_plasmid"], row["observed_ESC1"]
        ) else "mismatch"
        details.append(row)
    return details


def write_details(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("counts", type=Path, help="DotMatch MAGeCK-format count table")
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--details-out", type=Path)
    parser.add_argument("--summary-out", type=Path)
    args = parser.parse_args()

    try:
        details = compare_rows(args.reference, args.counts)
    except (OSError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc
    matched = sum(row["status"] == "match" for row in details)
    summary = {
        "status": "pass" if matched == len(details) and details else "fail",
        "reference_rows": len(details),
        "matching_rows": matched,
        "matching_sample_counts": matched * 2,
        "scope": "three-guide, two-sample external-reference spot check",
    }
    if args.details_out:
        write_details(args.details_out, details)
    if args.summary_out:
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        args.summary_out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))
    if summary["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
