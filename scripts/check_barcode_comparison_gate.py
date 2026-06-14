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
REPORT = ROOT / "docs" / "benchmarks" / "barcode_demux" / "README.md"


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"missing required artifact: {path}")
    with path.open() as fh:
        return list(csv.DictReader(fh))


def as_int(value: str | None, default: int = 0) -> int:
    if not value:
        return default
    return int(float(value))


def as_float(value: str | None, default: float = 0.0) -> float:
    if not value:
        return default
    return float(value)


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def mean_rate(rows: list[dict[str, str]]) -> float:
    rates = [as_float(row.get("reads_per_sec")) for row in rows if as_float(row.get("reads_per_sec")) > 0.0]
    return sum(rates) / len(rates) if rates else 0.0


def speedup(tool_a: str, rows_a: list[dict[str, str]], tool_b: str, rows_b: list[dict[str, str]],
            min_speedup: float, workflow: str, failures: list[str]) -> dict[str, str] | None:
    rate_a = mean_rate(rows_a)
    rate_b = mean_rate(rows_b)
    require(rate_a > 0.0 and rate_b > 0.0,
            f"{workflow}: {tool_a} and {tool_b} need positive reads/sec for speed comparison",
            failures)
    if rate_a > 0.0 and rate_b > 0.0:
        ratio = rate_a / rate_b
        require(ratio >= min_speedup,
                f"{workflow}: {tool_a}/{tool_b} speedup below {min_speedup:.2f}x: {ratio:.2f}x",
                failures)
        return {
            "workflow": workflow,
            "comparator": tool_b,
            "dotmatch_reads_per_sec": f"{rate_a:.1f}",
            "comparator_reads_per_sec": f"{rate_b:.1f}",
            "speedup": f"{ratio:.2f}x",
            "gate_floor": f"{min_speedup:.2f}x",
        }
    return None


def require_speedup(tool_a: str, rows_a: list[dict[str, str]], tool_b: str, rows_b: list[dict[str, str]],
                    min_speedup: float, workflow: str, failures: list[str]) -> None:
    speedup(tool_a, rows_a, tool_b, rows_b, min_speedup, workflow, failures)


def gated_speedups(rows: list[dict[str, str]], failures: list[str]) -> list[dict[str, str]]:
    ok = [
        row for row in rows
        if row.get("exit_code") == "0" and "fixture" not in row.get("workflow", "").lower() and
        "real" in row.get("workflow", "").lower()
    ]
    real_by_workflow: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in ok:
        real_by_workflow[row.get("workflow", "")].append(row)
    out: list[dict[str, str]] = []
    for workflow, workflow_rows in real_by_workflow.items():
        dotmatch_rows = [row for row in workflow_rows if row.get("tool") == "dotmatch_demux"]
        if not dotmatch_rows:
            continue
        cutadapt_rows = [row for row in workflow_rows if row.get("tool") == "cutadapt_demux"]
        if cutadapt_rows:
            entry = speedup("dotmatch_demux", dotmatch_rows, "cutadapt_demux", cutadapt_rows,
                            5.0, workflow, failures)
            if entry:
                out.append(entry)
        if any(row.get("k") == "0" for row in dotmatch_rows):
            hash_rows = [row for row in workflow_rows if row.get("tool") == "hash_splitter_exact"]
            if hash_rows:
                entry = speedup("dotmatch_demux", dotmatch_rows, "hash_splitter_exact", hash_rows,
                                3.0, workflow, failures)
                if entry:
                    out.append(entry)
        if any(row.get("k") == "1" and row.get("metric") == "hamming" for row in dotmatch_rows):
            hamming_rows = [row for row in workflow_rows if row.get("tool") == "hamming_radius_splitter"]
            if hamming_rows:
                entry = speedup("dotmatch_demux", dotmatch_rows, "hamming_radius_splitter", hamming_rows,
                                12.0, workflow, failures)
                if entry:
                    out.append(entry)
    return out


def report_gate(report: Path, speedups: list[dict[str, str]], failures: list[str]) -> None:
    require(report.exists(), f"missing barcode benchmark report: {report}", failures)
    if not report.exists():
        return
    text = report.read_text(encoding="utf-8")
    for required in [
        "## Gated Real-Data Speedups",
        "not public real-data evidence",
        "Levenshtein indel lane is fixture evidence",
    ]:
        require(required in text, f"barcode benchmark report must mention: {required}", failures)
    for row in speedups:
        expected = (
            f"| {row['workflow']} | {row['comparator']} | {row['dotmatch_reads_per_sec']} | "
            f"{row['comparator_reads_per_sec']} | {row['speedup']} | {row['gate_floor']} |"
        )
        require(expected in text, f"barcode benchmark report missing gated speedup row: {expected}", failures)


def metadata_gate(path: Path, failures: list[str]) -> None:
    require(path.exists(), f"missing barcode metadata: {path}", failures)
    if not path.exists():
        return
    meta = json.loads(path.read_text())
    require(bool(meta.get("evidence_ready")), "barcode metadata is not evidence-ready: real barcode sheet is missing", failures)
    require(as_int(str(meta.get("barcode_count", "0"))) > 0, "barcode metadata has zero parsed barcodes", failures)
    length_mode = str(meta.get("barcode_length_mode") or "").strip()
    has_fixed_length = as_int(str(meta.get("barcode_length", "0"))) > 0
    has_auto_lengths = length_mode == "auto" and bool(meta.get("barcode_lengths"))
    require(has_fixed_length or has_auto_lengths,
            "barcode metadata must declare a fixed barcode length or barcode length mode auto with parsed lengths",
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
    require(bool(real), "barcode comparison rows must use a real FASTQ workflow, not fixture/smoke rows", failures)
    by_tool: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in real:
        by_tool[row.get("tool", "")].append(row)
    require(len(by_tool.get("dotmatch_demux", [])) >= min_repeats,
            f"dotmatch_demux needs >= {min_repeats} successful real-data repeats", failures)
    if require_cutadapt:
        require(len(by_tool.get("cutadapt_demux", [])) >= min_repeats,
                f"cutadapt_demux needs >= {min_repeats} successful real-data repeats", failures)
    if require_second_comparator:
        second = sum(len(by_tool.get(tool, [])) for tool in ["ultraplex_demux", "je_demux", "hamming_radius_splitter"])
        second += sum(1 for row in by_tool.get("hash_splitter_exact", []) if as_int(row.get("k")) == 0)
        require(second >= min_repeats,
                "barcode comparison needs a second successful comparator row: Ultraplex, Je, hamming-radius splitter, or exact hash splitter for k=0",
                failures)
    for row in by_tool.get("dotmatch_demux", []):
        require(as_int(row.get("n_reads")) > 0, "DotMatch barcode row has zero reads", failures)
        require(as_int(row.get("n_barcodes")) > 0, "DotMatch barcode row has zero barcodes", failures)
        require(as_int(row.get("assigned_reads")) > 0, "DotMatch barcode row assigned zero reads", failures)
        if row.get("barcode_length") not in {"", "auto"} and row.get("metric") == "hamming" and row.get("k") in {"0", "1"}:
            require(str(row.get("assignment_engine", "")).startswith("hamming_"),
                    "fixed-length DotMatch Hamming barcode rows must record a direct hamming assignment_engine",
                    failures)
        if row.get("barcode_length") not in {"", "auto"} and row.get("metric") == "levenshtein" and row.get("k") == "1":
            require(row.get("assignment_engine") == "levenshtein_k1_lookup_direct",
                    "fixed-length DotMatch Levenshtein barcode rows must record assignment_engine=levenshtein_k1_lookup_direct",
                    failures)
    for tool in ["cutadapt_demux", "hash_splitter_exact", "hamming_radius_splitter",
                 "levenshtein_radius_splitter", "ultraplex_demux", "je_demux"]:
        for row in by_tool.get(tool, []):
            require(as_int(row.get("assigned_reads")) > 0, f"{tool} barcode row assigned zero reads", failures)
    for row in by_tool.get("hamming_radius_splitter", []):
        matching = [
            dotmatch for dotmatch in by_tool.get("dotmatch_demux", [])
            if dotmatch.get("workflow") == row.get("workflow") and dotmatch.get("k") == row.get("k")
        ]
        require(bool(matching), "hamming-radius splitter row has no matching DotMatch row", failures)
        for dotmatch in matching:
            require(as_int(dotmatch.get("assigned_reads")) == as_int(row.get("assigned_reads")),
                    "DotMatch barcode assigned_reads must match hamming-radius splitter", failures)
            require(as_int(dotmatch.get("corrected_reads")) == as_int(row.get("corrected_reads")),
                    "DotMatch barcode corrected_reads must match hamming-radius splitter", failures)
            require(as_int(dotmatch.get("ambiguous_reads")) == as_int(row.get("ambiguous_reads")),
                    "DotMatch barcode ambiguous_reads must match hamming-radius splitter", failures)
    for row in by_tool.get("levenshtein_radius_splitter", []):
        matching = [
            dotmatch for dotmatch in by_tool.get("dotmatch_demux", [])
            if dotmatch.get("workflow") == row.get("workflow") and dotmatch.get("k") == row.get("k") and
            dotmatch.get("metric") == "levenshtein"
        ]
        require(bool(matching), "Levenshtein splitter row has no matching DotMatch row", failures)
        for dotmatch in matching:
            require(as_int(dotmatch.get("assigned_reads")) == as_int(row.get("assigned_reads")),
                    "DotMatch barcode assigned_reads must match Levenshtein splitter", failures)
            require(as_int(dotmatch.get("corrected_reads")) == as_int(row.get("corrected_reads")),
                    "DotMatch barcode corrected_reads must match Levenshtein splitter", failures)
            require(as_int(dotmatch.get("ambiguous_reads")) == as_int(row.get("ambiguous_reads")),
                    "DotMatch barcode ambiguous_reads must match Levenshtein splitter", failures)
    gated_speedups(real, failures)


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
    parser.add_argument("--report", default=str(REPORT))
    parser.add_argument("--skip-report", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.smoke:
        args.min_repeats = 1
        args.require_cutadapt = False
        args.require_second_comparator = False
        args.skip_metadata = True
        args.skip_report = True

    failures: list[str] = []
    if not args.skip_metadata:
        metadata_gate(Path(args.metadata), failures)
    rows = read_rows(Path(args.rows))
    row_gate(rows, args.min_repeats, args.require_cutadapt,
             args.require_second_comparator, args.smoke, failures)
    if not args.skip_report:
        report_gate(Path(args.report), gated_speedups(rows, failures), failures)
    if failures:
        print("BARCODE comparison GATE: FAIL")
        for failure in failures:
            print(f"- {failure}")
        raise SystemExit(1)
    print("BARCODE comparison GATE: PASS")


if __name__ == "__main__":
    main()
