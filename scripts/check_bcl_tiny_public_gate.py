#!/usr/bin/env python3
"""Gate the narrow public 10x tiny-BCL milestone evidence.

This verifier is intentionally narrower than `check_bcl_comparison_gate.py`.
It proves that the public classic-BCL demo row, DotMatch output counts, and
available bcl2fastq validation are present. It does not authorize broad BCL
Convert, CBCL, NovaSeq, or production demultiplexing comparison wording.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "bcl_demux.csv"
REPORT = ROOT / "docs" / "benchmarks" / "bcl_demux" / "README.md"
PUBLIC_WORKFLOW = "public_10x_tiny_bcl"


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


def _hash_present(row: dict[str, str], key: str) -> bool:
    value = row.get(key, "")
    return len(value) >= 32 and all(char in "0123456789abcdef" for char in value.lower())


def row_gate(rows: list[dict[str, str]], failures: list[str]) -> None:
    public_rows = [row for row in rows if row.get("workflow") == PUBLIC_WORKFLOW]
    dotmatch = [
        row for row in public_rows
        if row.get("tool") == "dotmatch_bcl_demux"
        and row.get("format") == "classic_bcl"
        and row.get("exit_code") == "0"
    ]
    validated_bcl2fastq = [
        row for row in public_rows
        if row.get("tool") == "bcl2fastq"
        and row.get("format") == "classic_bcl"
        and row.get("exit_code") == "0"
        and row.get("validation_exit_code") == "0"
        and row.get("validation_mismatches") == "0"
    ]
    if not dotmatch:
        failures.append("missing successful DotMatch public 10x tiny-BCL classic row")
        return
    row = dotmatch[0]
    if as_int(row, "clusters") <= 0 or as_int(row, "total_clusters") <= 0:
        failures.append("public tiny-BCL DotMatch row must record positive cluster counts")
    if as_int(row, "clusters") != as_int(row, "total_clusters"):
        failures.append("public tiny-BCL DotMatch clusters must match total_clusters")
    if as_int(row, "cycles") <= 0 or as_int(row, "samples") <= 0 or as_int(row, "tiles") <= 0:
        failures.append("public tiny-BCL DotMatch row must record cycles, samples, and tiles")
    if as_int(row, "assigned_reads") <= 0:
        failures.append("public tiny-BCL DotMatch row must assign reads")
    if as_int(row, "assigned_reads") + as_int(row, "undetermined_reads") + as_int(row, "filtered_clusters") != as_int(row, "total_clusters"):
        failures.append("public tiny-BCL DotMatch assigned, undetermined, and filtered counts must sum to total clusters")
    if not _hash_present(row, "output_sha256") or not _hash_present(row, "fastq_content_sha256"):
        failures.append("public tiny-BCL DotMatch row must record output and FASTQ-content hashes")
    command = row.get("command", "")
    if "dotmatch bcl-demux" not in command:
        failures.append("public tiny-BCL DotMatch row must record the dotmatch bcl-demux command")
    if "cellranger-tiny-bcl-1.2.0" not in command or "cellranger-tiny-bcl-samplesheet.normalized.csv" not in command:
        failures.append("public tiny-BCL DotMatch command must name the public run folder and normalized sample sheet")

    if not validated_bcl2fastq:
        failures.append("missing validated bcl2fastq public 10x tiny-BCL comparator row")
        return
    comparator = validated_bcl2fastq[0]
    if as_int(comparator, "assigned_reads") != as_int(row, "assigned_reads"):
        failures.append("public tiny-BCL bcl2fastq assigned read count must agree with DotMatch")
    if as_int(comparator, "undetermined_reads") != as_int(row, "undetermined_reads"):
        failures.append("public tiny-BCL bcl2fastq undetermined read count must agree with DotMatch")
    if comparator.get("validation_mode") != "count_totals":
        failures.append("public tiny-BCL bcl2fastq row must declare count_totals validation mode")
    if not _hash_present(comparator, "fastq_content_sha256"):
        failures.append("public tiny-BCL bcl2fastq row must record a FASTQ-content hash")


def report_gate(path: Path, failures: list[str]) -> None:
    if not path.exists():
        failures.append(f"missing BCL benchmark report: {path}")
        return
    text = path.read_text(encoding="utf-8")
    if PUBLIC_WORKFLOW not in text:
        failures.append("BCL benchmark report must name public_10x_tiny_bcl")
    if "not a comparison result by itself" not in text:
        failures.append("BCL benchmark report must state the public tiny-BCL row is not a comparison result by itself")
    if "make bcl-tiny-public-gate" not in text:
        failures.append("BCL benchmark report must name make bcl-tiny-public-gate for the narrow public milestone")
    if "make bcl-comparison-gate" not in text:
        failures.append("BCL benchmark report must keep broader BCL comparison wording behind make bcl-comparison-gate")
    if "successful DotMatch CBCL row" not in text:
        failures.append("BCL benchmark report must document that broader comparison requires a successful DotMatch CBCL row")
    if "distinct repeated" not in text:
        failures.append("BCL benchmark report must document that broader comparison requires distinct repeated timing")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default=str(RAW))
    parser.add_argument("--report", default=str(REPORT))
    args = parser.parse_args(argv)

    failures: list[str] = []
    row_gate(read_rows(Path(args.csv)), failures)
    report_gate(Path(args.report), failures)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        print("BCL TINY PUBLIC: FAIL")
        return 1
    print("BCL TINY PUBLIC: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
