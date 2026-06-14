#!/usr/bin/env python3
"""Generate native Edlib comparison graphs for README benchmarks."""

from __future__ import annotations

import platform
import subprocess
import os
import csv
import statistics
from io import StringIO
from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "docs" / "benchmarks" / "native"
RAW_DIR = ROOT / "benchmarks" / "raw"
FIG_DIR = ROOT / "benchmarks" / "figures"
REPORT_READS = int(os.environ.get("DOTMATCH_NATIVE_REPORT_READS", os.environ.get("QDALN_NATIVE_REPORT_READS", "1000")))
REPORT_REPEATS = int(os.environ.get("DOTMATCH_NATIVE_REPEATS", "3"))


def run_command(cmd: list[str]) -> str:
    result = subprocess.run(cmd, cwd=ROOT, check=True, text=True, capture_output=True)
    return result.stdout


def _fnum(row: dict[str, object], key: str) -> float:
    try:
        return float(row.get(key, "") or 0.0)
    except (TypeError, ValueError):
        return 0.0


def _snum(row: dict[str, object], key: str) -> str:
    value = row.get(key, "")
    return "" if value is None else str(value)


def _quantile(values: list[float], q: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    pos = (len(ordered) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(ordered) - 1)
    frac = pos - lo
    return ordered[lo] * (1.0 - frac) + ordered[hi] * frac


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, object]], columns: list[str], floatfmt: str = ".2f") -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        vals = []
        for col in columns:
            value = row.get(col, "")
            if col in {"n_reads", "n_targets", "len", "k"}:
                vals.append(str(int(value)))
            elif col in {"err", "indel_rate"}:
                vals.append(format(float(value), ".3f"))
            elif isinstance(value, float):
                vals.append(format(value, floatfmt))
            else:
                vals.append(str(value))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def load_native() -> list[dict[str, object]]:
    run_command(["make", "build/bench_edlib_native"])
    rows: list[dict[str, object]] = []
    raw_outputs = []
    for repeat in range(REPORT_REPEATS):
        output = run_command(["./build/bench_edlib_native", str(REPORT_READS)])
        raw_outputs.append(output)
        for row in csv.DictReader(StringIO(output)):
            row["repeat"] = str(repeat)
            rows.append(row)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    _write_csv(OUT_DIR / "native_edlib_assignment.csv", rows)
    _write_csv(RAW_DIR / "native_edlib_assignment.csv", rows)
    (RAW_DIR / "native_edlib_assignment_raw_runs.csv").write_text("\n".join(raw_outputs), encoding="utf-8")
    return rows


def speedup_frame(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    keys = ["workload", "error_mode", "n_reads", "n_targets", "len", "k", "err", "indel_rate", "repeat"]
    selected = [row for row in rows if row.get("tool") == "dotmatch_indexed" and row.get("k") != "0"]
    exact_by_key: dict[tuple[str, ...], dict[str, object]] = {}
    for row in rows:
        if row.get("tool") not in {"dotmatch_exact_direct", "dotmatch_exact_batch"}:
            continue
        key = tuple(_snum(row, field) for field in keys)
        current = exact_by_key.get(key)
        if current is None or _fnum(row, "reads_per_sec") > _fnum(current, "reads_per_sec"):
            exact_by_key[key] = row
    selected.extend(exact_by_key.values())
    edlib = {tuple(_snum(row, field) for field in keys): row for row in rows if row.get("tool") == "edlib_native_scan"}
    out: list[dict[str, object]] = []
    for row in selected:
        key = tuple(_snum(row, field) for field in keys)
        edlib_row = edlib.get(key)
        if edlib_row is None or _fnum(edlib_row, "reads_per_sec") <= 0.0:
            continue
        merged = {field: row.get(field, "") for field in keys}
        merged.update({
            "reads_per_sec_dotmatch": _fnum(row, "reads_per_sec"),
            "reads_per_sec_edlib": _fnum(edlib_row, "reads_per_sec"),
            "candidates_per_read": _fnum(row, "candidates_per_read"),
            "verified_per_read": _fnum(row, "verified_per_read"),
            "peak_rss_kb": _fnum(row, "peak_rss_kb"),
            "mismatches": _fnum(row, "mismatches"),
            "dotmatch_tool": row.get("tool", ""),
            "speedup_vs_edlib_native": _fnum(row, "reads_per_sec") / _fnum(edlib_row, "reads_per_sec"),
        })
        out.append(merged)
    return out


def aggregate_stats(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    keys = ["tool", "workload", "error_mode", "n_targets", "len", "k", "err", "indel_rate"]
    groups: dict[tuple[str, ...], list[dict[str, object]]] = {}
    for row in rows:
        groups.setdefault(tuple(_snum(row, field) for field in keys), []).append(row)
    out: list[dict[str, object]] = []
    for key, group in sorted(groups.items()):
        rps = [_fnum(row, "reads_per_sec") for row in group]
        seconds = [_fnum(row, "seconds") for row in group]
        verified = [_fnum(row, "verified_per_read") for row in group]
        peak = [_fnum(row, "peak_rss_kb") for row in group]
        mean = statistics.mean(rps) if rps else 0.0
        std = statistics.stdev(rps) if len(rps) > 1 else 0.0
        item: dict[str, object] = dict(zip(keys, key))
        item.update({
            "repeats": len({_snum(row, "repeat") for row in group}),
            "reads_per_sec_mean": mean,
            "reads_per_sec_p50": statistics.median(rps) if rps else 0.0,
            "reads_per_sec_p95": _quantile(rps, 0.95),
            "reads_per_sec_std": std,
            "seconds_mean": statistics.mean(seconds) if seconds else 0.0,
            "verified_per_read_median": statistics.median(verified) if verified else 0.0,
            "peak_rss_kb_max": max(peak) if peak else 0.0,
            "mismatches_sum": sum(_fnum(row, "mismatches") for row in group),
            "reads_per_sec_cv": std / mean if mean else 0.0,
        })
        out.append(item)
    return out


def gated_evidence_summary(speedups: list[dict[str, object]]) -> list[dict[str, object]]:
    definitions = [
        {
            "claim": "k=1 substitution indexed rows",
            "k": "1",
            "error_mode": "one_substitution",
            "min_speedup_required": 10.0,
            "max_verified_required": 1.05,
        },
        {
            "claim": "k=2 substitution indexed rows",
            "k": "2",
            "error_mode": "one_substitution",
            "min_speedup_required": 10.0,
            "max_verified_required": 1.05,
        },
        {
            "claim": "Levenshtein k=2 insertion/deletion rows",
            "k": "2",
            "error_mode": "one_insertion|one_deletion",
            "min_speedup_required": 8.0,
            "max_verified_required": 25.0,
        },
    ]
    out: list[dict[str, object]] = []
    for definition in definitions:
        modes = set(str(definition["error_mode"]).split("|"))
        rows = [
            row for row in speedups
            if _snum(row, "dotmatch_tool") == "dotmatch_indexed"
            and _snum(row, "k") == definition["k"]
            and _snum(row, "error_mode") in modes
            and int(float(_snum(row, "n_targets") or "0")) >= 4096
        ]
        if not rows:
            continue
        speedups_values = [_fnum(row, "speedup_vs_edlib_native") for row in rows]
        verified_values = [_fnum(row, "verified_per_read") for row in rows]
        out.append({
            "claim": definition["claim"],
            "large_library_rows": len(rows),
            "min_speedup_vs_edlib": min(speedups_values),
            "median_speedup_vs_edlib": statistics.median(speedups_values),
            "max_verified_per_read": max(verified_values),
            "min_speedup_required": definition["min_speedup_required"],
            "max_verified_required": definition["max_verified_required"],
        })
    return out


def _group_values(rows: list[dict[str, object]], group_key: str, value_key: str) -> dict[str, list[float]]:
    grouped: dict[str, list[float]] = {}
    for row in rows:
        grouped.setdefault(_snum(row, group_key), []).append(_fnum(row, value_key))
    return grouped


def plot_speedup(speedups: list[dict[str, object]]) -> None:
    if plt is None:
        return
    subset = [
        row for row in speedups
        if _snum(row, "len") in {"16", "32"}
        and _snum(row, "err") in {"0.000", "0.010", "0", "0.01"}
        and _snum(row, "error_mode") in {"exact", "one_substitution"}
    ]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    series: dict[tuple[str, str, str], list[float]] = {}
    for row in subset:
        key = (_snum(row, "len"), _snum(row, "k"), _snum(row, "n_targets"))
        series.setdefault(key, []).append(_fnum(row, "speedup_vs_edlib_native"))
    by_line: dict[tuple[str, str], list[tuple[int, float]]] = {}
    for (length, k, n_targets), values in series.items():
        by_line.setdefault((length, k), []).append((int(n_targets), statistics.median(values)))
    for (length, k), points in sorted(by_line.items()):
        points.sort()
        ax.plot([p[0] for p in points], [p[1] for p in points], marker="o", label=f"len={length} k={k}")
    ax.axhline(10.0, color="#8b0000", linestyle="--", linewidth=1, label="10x target")
    ax.set_xscale("log")
    ax.set_title("Indexed assignment speedup vs native Edlib scan")
    ax.set_xlabel("number of targets")
    ax.set_ylabel("reads/sec speedup")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, ncols=2)
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / "native_speedup_vs_edlib.svg")
    fig.savefig(FIG_DIR / "native_speedup_vs_edlib.svg")
    fig.savefig(FIG_DIR / "native_speedup_vs_edlib.pdf")
    plt.close(fig)


def plot_candidates(speedups: list[dict[str, object]]) -> None:
    if plt is None:
        return
    subset = [
        row for row in speedups
        if _snum(row, "len") == "32"
        and _snum(row, "err") in {"0.010", "0.01"}
        and _snum(row, "error_mode") == "one_substitution"
    ]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    series: dict[tuple[str, str], list[float]] = {}
    for row in subset:
        series.setdefault((_snum(row, "k"), _snum(row, "n_targets")), []).append(_fnum(row, "verified_per_read"))
    by_line: dict[str, list[tuple[int, float]]] = {}
    for (k, n_targets), values in series.items():
        by_line.setdefault(k, []).append((int(n_targets), statistics.median(values)))
    for k, points in sorted(by_line.items()):
        points.sort()
        ax.plot([p[0] for p in points], [p[1] for p in points], marker="o", label=f"k={k}")
    ax.set_xscale("log")
    ax.set_title("DotMatch indexed candidates verified per read")
    ax.set_xlabel("number of targets")
    ax.set_ylabel("verified candidates/read")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / "native_candidates_per_read.svg")
    fig.savefig(FIG_DIR / "native_candidates_per_read.svg")
    fig.savefig(FIG_DIR / "native_candidates_per_read.pdf")
    plt.close(fig)


def plot_throughput(rows: list[dict[str, object]]) -> None:
    if plt is None:
        return
    subset = [
        row for row in rows
        if _snum(row, "len") == "32"
        and _snum(row, "k") == "1"
        and _snum(row, "err") in {"0.010", "0.01"}
        and _snum(row, "error_mode") == "one_substitution"
    ]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    series: dict[tuple[str, str], list[float]] = {}
    for row in subset:
        series.setdefault((_snum(row, "tool"), _snum(row, "n_targets")), []).append(_fnum(row, "reads_per_sec"))
    by_line: dict[str, list[tuple[int, float]]] = {}
    for (tool, n_targets), values in series.items():
        by_line.setdefault(tool, []).append((int(n_targets), statistics.median(values)))
    for tool, points in sorted(by_line.items()):
        points.sort()
        ax.plot([p[0] for p in points], [p[1] for p in points], marker="o", label=tool)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("Assignment throughput, len=32 k=1 err=1%")
    ax.set_xlabel("number of targets")
    ax.set_ylabel("reads/sec")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / "native_assignment_throughput.svg")
    fig.savefig(FIG_DIR / "native_assignment_throughput.svg")
    fig.savefig(FIG_DIR / "native_assignment_throughput.pdf")
    plt.close(fig)


def write_report(rows: list[dict[str, object]], speedups: list[dict[str, object]]) -> None:
    agg = aggregate_stats(rows)
    _write_csv(RAW_DIR / "native_edlib_assignment_summary.csv", agg)

    best_groups: dict[tuple[str, str, str, str, str], list[dict[str, object]]] = {}
    for row in speedups:
        key = tuple(_snum(row, field) for field in ["n_targets", "len", "k", "err", "error_mode"])
        best_groups.setdefault(key, []).append(row)
    best: list[dict[str, object]] = []
    for (n_targets, length, k, err, error_mode), group in best_groups.items():
        best.append({
            "dotmatch_tool": group[0].get("dotmatch_tool", ""),
            "n_targets": n_targets,
            "len": length,
            "k": k,
            "error_mode": error_mode,
            "err": err,
            "reads_per_sec_dotmatch": statistics.median([_fnum(row, "reads_per_sec_dotmatch") for row in group]),
            "reads_per_sec_edlib": statistics.median([_fnum(row, "reads_per_sec_edlib") for row in group]),
            "verified_per_read": statistics.median([_fnum(row, "verified_per_read") for row in group]),
            "peak_rss_kb": max(_fnum(row, "peak_rss_kb") for row in group),
            "mismatches": sum(_fnum(row, "mismatches") for row in group),
            "speedup_vs_edlib_native": statistics.median([_fnum(row, "speedup_vs_edlib_native") for row in group]),
        })
    best = sorted(best, key=lambda row: _fnum(row, "speedup_vs_edlib_native"), reverse=True)[:12]

    summary_groups: dict[tuple[str, str, str, str], list[float]] = {}
    for row in speedups:
        key = tuple(_snum(row, field) for field in ["len", "k", "n_targets", "error_mode"])
        summary_groups.setdefault(key, []).append(_fnum(row, "speedup_vs_edlib_native"))
    summary = [
        {"len": length, "k": k, "n_targets": n_targets, "error_mode": error_mode,
         "speedup_vs_edlib_native": statistics.median(values)}
        for (length, k, n_targets, error_mode), values in summary_groups.items()
    ]
    summary = sorted(summary, key=lambda row: _fnum(row, "speedup_vs_edlib_native"), reverse=True)[:12]
    gated = gated_evidence_summary(speedups)
    zero_mismatch = int(sum(_fnum(row, "mismatches") for row in rows))
    lines = [
        "# Native Edlib Benchmark Report",
        "",
        f"- Platform: `{platform.platform()}`",
        f"- Python: `{platform.python_version()}`",
        f"- Reads per benchmark case: `{REPORT_READS}`",
        f"- Repetitions per benchmark case: `{REPORT_REPEATS}`",
        "- Comparator: native Edlib C/C++ API, `EDLIB_MODE_NW`, `EDLIB_TASK_DISTANCE`, fixed threshold `k`.",
        "- Additional baselines: exact hash lookup for `k=0`; BK-tree and neighbor lookup approximate baselines for `k=1`.",
        "- Gate: `make native-exact-gate` requires zero mismatches, large-library exact rows to beat `exact_hash_lookup`, large-library indexed `k=1` rows to beat exhaustive Edlib by >10x, beat the best BK-tree/neighbor baseline, and verify no more than 1.05 candidates/read, plus large-library `k=2` substitution rows to beat exhaustive Edlib by >10x with no more than 1.05 verified candidates/read and Levenshtein `k=2` insertion/deletion rows to beat exhaustive Edlib by >8x while verifying no more than 25 candidates/read.",
        f"- Assignment mismatches recorded across all rows: `{zero_mismatch}`.",
        "- Every benchmark run aborts on assignment disagreement between DotMatch and native Edlib scan.",
        "",
        "![Native speedup vs Edlib](native_speedup_vs_edlib.svg)",
        "",
        "![Native candidates per read](native_candidates_per_read.svg)",
        "",
        "![Native assignment throughput](native_assignment_throughput.svg)",
        "",
        "## Gated Native Scaling Claims",
        "",
        markdown_table(gated, ["claim", "large_library_rows", "min_speedup_vs_edlib", "median_speedup_vs_edlib", "max_verified_per_read", "min_speedup_required", "max_verified_required"], ".2f"),
        "",
        "## Highest Observed Microbenchmark Speedups",
        "",
        markdown_table(best, ["dotmatch_tool", "n_targets", "len", "k", "error_mode", "err", "reads_per_sec_dotmatch", "reads_per_sec_edlib", "verified_per_read", "peak_rss_kb", "speedup_vs_edlib_native"], ".2f"),
        "",
        "## Median Speedup Summary",
        "",
        markdown_table(summary, ["len", "k", "n_targets", "error_mode", "speedup_vs_edlib_native"], ".2f"),
        "",
        "## Repeated-Run Statistics",
        "",
        markdown_table(
            [row for row in agg if row.get("tool") in {"dotmatch_exact_direct", "dotmatch_exact_batch", "dotmatch_indexed", "exact_hash_lookup", "edlib_native_scan"}][:24],
            ["tool", "error_mode", "n_targets", "len", "k", "err", "reads_per_sec_mean", "reads_per_sec_p50", "reads_per_sec_p95", "reads_per_sec_cv", "peak_rss_kb_max", "mismatches_sum"],
            ".2f",
        ),
        "",
        "## Evidence Boundary",
        "",
        "These are native Edlib scan microbenchmarks for exact short-DNA assignment workloads, plus simple exact-hash and BK-tree/neighbor baselines. The largest rows are useful for understanding algorithmic scaling against exhaustive scan, but they are not end-to-end workflow speed claims. Exact `k=0` lookup should be judged against hash-table baselines: broad exact-hash superiority is not claimed unless the exact gate proves it, while large-library exact rows (`n_targets >= 4096`) may be described only when `make native-exact-gate` records a >1.0 ratio against `exact_hash_lookup`. For `k=1`, large-library indexed rows may be described as non-exhaustive only when the same gate records zero Edlib disagreements, >10x speedup over exhaustive Edlib scan, >1.0 speedup over the best BK-tree/neighbor baseline, and no more than 1.05 verified candidates/read. Fixed-length `k=2` substitution rows may be described as non-exhaustive only when the same gate records zero Edlib disagreements, >10x speedup over exhaustive Edlib scan, and no more than 1.05 verified candidates/read. Levenshtein `k=2` insertion/deletion rows may be described as non-exhaustive only when the same gate records zero Edlib disagreements, >8x speedup over exhaustive Edlib scan, and no more than 25 verified candidates/read. This remains scoped to packed A/C/G/T fixed-window assignment up to 32 bases with fallback preserving semantics for unsupported cases.",
        "",
    ]
    (OUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    df = load_native()
    speedups = speedup_frame(df)
    plot_speedup(speedups)
    plot_candidates(speedups)
    plot_throughput(df)
    write_report(df, speedups)
    print(f"wrote native benchmark report to {OUT_DIR}")


if __name__ == "__main__":
    main()
