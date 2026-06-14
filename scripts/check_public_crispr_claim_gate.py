#!/usr/bin/env python3
"""Fail unless public CRISPR benchmark artifacts support public evidence statements."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw"
REPORT = ROOT / "docs" / "benchmarks" / "public_crispr" / "README.md"


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"missing required artifact: {path}")
    with path.open() as fh:
        return list(csv.DictReader(fh))


def as_int(value: str | None, default: int = 0) -> int:
    if value is None or value == "":
        return default
    if value == "full":
        return 10**18
    return int(float(value))


def as_float(value: str | None, default: float = 0.0) -> float:
    if value is None or value == "":
        return default
    return float(value)


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def repeated_stats(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row.get("exit_code") != "0":
            continue
        key = (row.get("tool", ""), row.get("semantics", ""), row.get("requested_records_per_sample", row.get("n_reads", "")))
        groups[key].append(row)
    out: list[dict[str, str]] = []
    for (tool, semantics, requested), group in sorted(groups.items()):
        values = [as_float(row.get("reads_per_sec")) for row in group]
        if not values:
            continue
        out.append({
            "tool": tool,
            "semantics": semantics,
            "records_per_sample": requested,
            "mean_reads_per_sec": f"{statistics.mean(values):.1f}",
        })
    return out


def speedup_rows(stats: list[dict[str, str]], dotmatch_tool: str, baseline_tools: set[str], label: str) -> list[dict[str, str]]:
    by_size = {row["records_per_sample"]: row for row in stats if row["tool"] == dotmatch_tool}
    out: list[dict[str, str]] = []
    for row in stats:
        size = row["records_per_sample"]
        if size not in by_size or row["tool"] not in baseline_tools:
            continue
        dotmatch_rps = as_float(by_size[size].get("mean_reads_per_sec"))
        baseline_rps = as_float(row.get("mean_reads_per_sec"))
        if dotmatch_rps <= 0.0 or baseline_rps <= 0.0:
            continue
        out.append({
            "baseline": row["tool"],
            "records_per_sample": size,
            f"{label}_reads_per_sec": f"{dotmatch_rps:.1f}",
            "baseline_reads_per_sec": f"{baseline_rps:.1f}",
            "speedup": f"{dotmatch_rps / baseline_rps:.2f}x",
        })
    return out


def report_gate(report: Path, repeated_rows: list[dict[str, str]], require_guide_counter: bool,
                failures: list[str]) -> None:
    require(report.exists(), f"missing public CRISPR report: {report}", failures)
    if not report.exists():
        return
    text = report.read_text(encoding="utf-8")
    for required in [
        "## DotMatch Hamming Speedup",
        "## DotMatch Exact Count Speedup",
        "## Edlib Oracle Validation",
        "bounded Edlib validation",
        "guide-counter is fast",
        "DotMatch assigns at most one target per read",
    ]:
        require(required in text, f"public CRISPR report must mention: {required}", failures)
    stats = repeated_stats(repeated_rows)
    hamming_speedups = speedup_rows(
        stats,
        "dotmatch_hamming_k1",
        {"guide_counter_one_mismatch"} if require_guide_counter else set(),
        "dotmatch_hamming",
    )
    exact_speedups = speedup_rows(
        stats,
        "dotmatch_exact_k0",
        {"mageck_count_exact", "guide_counter_exact"} if require_guide_counter else {"mageck_count_exact"},
        "dotmatch_exact",
    )
    for row in hamming_speedups:
        expected = (
            f"| {row['baseline']} | {row['records_per_sample']} | {row['dotmatch_hamming_reads_per_sec']} | "
            f"{row['baseline_reads_per_sec']} | {row['speedup']} |"
        )
        require(expected in text, f"public CRISPR report missing Hamming speedup row: {expected}", failures)
    for row in exact_speedups:
        expected = (
            f"| {row['baseline']} | {row['records_per_sample']} | {row['dotmatch_exact_reads_per_sec']} | "
            f"{row['baseline_reads_per_sec']} | {row['speedup']} |"
        )
        require(expected in text, f"public CRISPR report missing exact speedup row: {expected}", failures)


def repeated_gate(rows: list[dict[str, str]], min_records: int, min_repeats: int,
                  require_guide_counter: bool, failures: list[str]) -> None:
    ok_rows = [r for r in rows if r.get("exit_code") == "0"]
    require(bool(ok_rows), "public_crispr_repeated.csv has no successful rows", failures)

    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in ok_rows:
        requested = row.get("requested_records_per_sample") or row.get("n_reads") or "0"
        if as_int(requested) >= min_records:
            groups[(row.get("tool", ""), requested)].append(row)

    required_tools = ["dotmatch_exact_k0", "dotmatch_hamming_k1", "dotmatch_levenshtein_k1", "mageck_count_exact"]
    if require_guide_counter:
        required_tools.append("guide_counter_one_mismatch")
        required_tools.append("guide_counter_exact")

    for tool in required_tools:
        best = max((len(group) for (t, _), group in groups.items() if t == tool), default=0)
        require(best >= min_repeats,
                f"{tool} needs >= {min_repeats} successful repeats at >= {min_records} records/sample; found {best}",
                failures)

    for (tool, requested), group in groups.items():
        if not tool.startswith("dotmatch"):
            continue
        verified = [as_float(r.get("verified_per_read"), -1.0) for r in group if r.get("verified_per_read")]
        if tool == "dotmatch_levenshtein_k1" and verified:
            require(max(verified) <= 5.0,
                    f"{tool} verified_per_read exceeds 5.0 at {requested}: max {max(verified):.3f}",
                    failures)


def validation_gate(rows: list[dict[str, str]], min_checked: int, failures: list[str]) -> None:
    require(bool(rows), "public_crispr_edlib_validation.csv is empty", failures)
    for row in rows:
        mismatches = as_int(row.get("mismatches"))
        checked = as_int(row.get("checked_reads"))
        require(mismatches == 0,
                f"Edlib validation mismatch for {row.get('sample', '')}: {mismatches}",
                failures)
        require(checked >= min_checked,
                f"Edlib validation checked_reads below {min_checked} for {row.get('sample', '')}: {checked}",
                failures)
        if "oracle_strategy" in row:
            require(row.get("oracle_strategy") == "bounded_edlib_candidates",
                    f"Edlib validation must use bounded_edlib_candidates for {row.get('sample', '')}: {row.get('oracle_strategy', '')}",
                    failures)
        if "bounded_windows" in row:
            require(as_int(row.get("bounded_windows")) > 0,
                    f"Edlib validation bounded_windows missing for {row.get('sample', '')}",
                    failures)
        if "fallback_windows" in row:
            fallback_windows = as_int(row.get("fallback_windows"))
            max_fallback = max(1, int(checked * 0.05))
            require(fallback_windows <= max_fallback,
                    f"Edlib validation fallback_windows exceeds 5% for {row.get('sample', '')}: {fallback_windows}/{checked}",
                    failures)
        if "edlib_alignments" in row:
            require(as_int(row.get("edlib_alignments")) > 0,
                    f"Edlib validation edlib_alignments missing for {row.get('sample', '')}",
                    failures)


def count_agreement_gate(rows: list[dict[str, str]], require_guide_counter: bool, failures: list[str]) -> None:
    by_name = {row.get("comparison", ""): row for row in rows}
    exact = by_name.get("dotmatch_exact_vs_mageck_exact")
    require(exact is not None and exact.get("status") == "ok", "missing exact DotMatch-vs-MAGeCK count agreement", failures)
    if exact is not None and exact.get("status") == "ok":
        require(as_int(exact.get("total_delta")) == 0, "DotMatch exact count total differs from MAGeCK exact", failures)
        require(as_int(exact.get("differing_guides")) == 0, "DotMatch exact has guide-level differences vs MAGeCK exact", failures)

    if require_guide_counter:
        hamming = by_name.get("dotmatch_hamming_vs_guide_counter")
        require(hamming is not None and hamming.get("status") == "ok", "missing DotMatch-vs-guide-counter count agreement", failures)
        if hamming is not None and hamming.get("status") == "ok":
            require(as_float(hamming.get("pearson")) >= 0.90,
                    f"DotMatch-vs-guide-counter Pearson below 0.90: {hamming.get('pearson')}",
                    failures)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeated", default=str(RAW / "public_crispr_repeated.csv"))
    parser.add_argument("--validation", default=str(RAW / "public_crispr_edlib_validation.csv"))
    parser.add_argument("--count-agreement", default=str(RAW / "count_agreement_summary.csv"))
    parser.add_argument("--report", default=str(REPORT))
    parser.add_argument("--skip-report", action="store_true")
    parser.add_argument("--min-records", type=int, default=10000)
    parser.add_argument("--min-repeats", type=int, default=5)
    parser.add_argument("--min-edlib-checked", type=int, default=1000)
    parser.add_argument("--require-guide-counter", action="store_true", default=True)
    parser.add_argument("--no-guide-counter", action="store_false", dest="require_guide_counter")
    parser.add_argument("--smoke", action="store_true", help="development gate: lower thresholds, no competitor requirement")
    args = parser.parse_args()

    if args.smoke:
        args.min_records = 25
        args.min_repeats = 1
        args.min_edlib_checked = 25
        args.require_guide_counter = False
        args.skip_report = True

    failures: list[str] = []
    repeated_rows = read_rows(Path(args.repeated))
    repeated_gate(repeated_rows, args.min_records, args.min_repeats,
                  args.require_guide_counter, failures)
    validation_gate(read_rows(Path(args.validation)), args.min_edlib_checked, failures)
    count_agreement_gate(read_rows(Path(args.count_agreement)), args.require_guide_counter, failures)
    if not args.skip_report:
        report_gate(Path(args.report), repeated_rows, args.require_guide_counter, failures)

    if failures:
        print("PUBLIC CRISPR EVIDENCE GATE: FAIL")
        for failure in failures:
            print(f"- {failure}")
        raise SystemExit(1)
    print("PUBLIC CRISPR EVIDENCE GATE: PASS")


if __name__ == "__main__":
    main()
