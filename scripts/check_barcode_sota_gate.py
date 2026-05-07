#!/usr/bin/env python3
"""Fail unless inline barcode/demux evidence is real-data and comparator-backed."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw"
META = ROOT / "examples" / "barcode_demux" / "data" / "metadata.json"


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"missing required artifact: {path}")
    with path.open() as fh:
        return list(csv.DictReader(fh))


def as_int(value: str | None, default: int = 0) -> int:
    if not value:
        return default
    return int(float(value))


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def metadata_gate(path: Path, failures: list[str]) -> None:
    require(path.exists(), f"missing barcode metadata: {path}", failures)
    if not path.exists():
        return
    meta = json.loads(path.read_text())
    require(bool(meta.get("claim_grade_ready")), "barcode metadata is not claim-grade: real barcode sheet is missing", failures)
    require(as_int(str(meta.get("barcode_count", "0"))) > 0, "barcode metadata has zero parsed barcodes", failures)
    require(as_int(str(meta.get("barcode_length", "0"))) > 0,
            "barcode metadata must declare a fixed barcode length for the current fixed-position benchmark",
            failures)
    require(bool(meta.get("runs")), "barcode metadata has no ENA run metadata", failures)
    for run in meta.get("runs", []):
        require(bool(run.get("local_md5") or run.get("ena", {}).get("fastq_md5")),
                f"barcode run lacks checksum metadata: {run.get('accession', '')}", failures)


def row_gate(rows: list[dict[str, str]], min_repeats: int, require_cutadapt: bool,
             require_second_comparator: bool, allow_fixture: bool, failures: list[str]) -> None:
    ok = [r for r in rows if r.get("exit_code") == "0"]
    require(bool(ok), "barcode_demux.csv has no successful rows", failures)
    real = ok if allow_fixture else [
        r for r in ok if "fixture" not in r.get("workflow", "").lower() and "real" in r.get("workflow", "").lower()
    ]
    require(bool(real), "barcode SOTA rows must use a real FASTQ workflow, not fixture/smoke rows", failures)
    by_tool: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in real:
        by_tool[row.get("tool", "")].append(row)
    require(len(by_tool.get("dotmatch_demux", [])) >= min_repeats,
            f"dotmatch_demux needs >= {min_repeats} successful real-data repeats", failures)
    if require_cutadapt:
        require(len(by_tool.get("cutadapt_demux", [])) >= min_repeats,
                f"cutadapt_demux needs >= {min_repeats} successful real-data repeats", failures)
    if require_second_comparator:
        second = sum(len(by_tool.get(tool, [])) for tool in ["ultraplex_demux", "je_demux", "hash_splitter_exact"])
        require(second >= min_repeats,
                "barcode SOTA needs a second successful comparator row: Ultraplex, Je, or exact hash splitter",
                failures)
    for row in by_tool.get("dotmatch_demux", []):
        require(as_int(row.get("n_reads")) > 0, "DotMatch barcode row has zero reads", failures)
        require(as_int(row.get("n_barcodes")) > 0, "DotMatch barcode row has zero barcodes", failures)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", default=str(RAW / "barcode_demux.csv"))
    parser.add_argument("--metadata", default=str(META))
    parser.add_argument("--min-repeats", type=int, default=5)
    parser.add_argument("--require-cutadapt", action="store_true", default=True)
    parser.add_argument("--no-cutadapt", action="store_false", dest="require_cutadapt")
    parser.add_argument("--require-second-comparator", action="store_true", default=True)
    parser.add_argument("--no-second-comparator", action="store_false", dest="require_second_comparator")
    parser.add_argument("--skip-metadata", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.smoke:
        args.min_repeats = 1
        args.require_cutadapt = False
        args.require_second_comparator = False
        args.skip_metadata = True

    failures: list[str] = []
    if not args.skip_metadata:
        metadata_gate(Path(args.metadata), failures)
    row_gate(read_rows(Path(args.rows)), args.min_repeats, args.require_cutadapt,
             args.require_second_comparator, args.smoke, failures)
    if failures:
        print("BARCODE SOTA GATE: FAIL")
        for failure in failures:
            print(f"- {failure}")
        raise SystemExit(1)
    print("BARCODE SOTA GATE: PASS")


if __name__ == "__main__":
    main()
