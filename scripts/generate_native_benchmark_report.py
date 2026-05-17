#!/usr/bin/env python3
"""Generate native Edlib comparison graphs for README benchmarks."""

from __future__ import annotations

import platform
import subprocess
import os
from io import StringIO
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "docs" / "benchmarks" / "native"
RAW_DIR = ROOT / "benchmarks" / "raw"
FIG_DIR = ROOT / "benchmarks" / "figures"
REPORT_READS = int(os.environ.get("DOTMATCH_NATIVE_REPORT_READS", os.environ.get("QDALN_NATIVE_REPORT_READS", "1000")))
REPORT_REPEATS = int(os.environ.get("DOTMATCH_NATIVE_REPEATS", "3"))


def run_command(cmd: list[str]) -> str:
    result = subprocess.run(cmd, cwd=ROOT, check=True, text=True, capture_output=True)
    return result.stdout


def markdown_table(df: pd.DataFrame, floatfmt: str = ".2f") -> str:
    columns = list(df.columns)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for _, row in df.iterrows():
        vals = []
        for col in columns:
            value = row[col]
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


def load_native() -> pd.DataFrame:
    run_command(["make", "build/bench_edlib_native"])
    frames = []
    raw_outputs = []
    for repeat in range(REPORT_REPEATS):
        output = run_command(["./build/bench_edlib_native", str(REPORT_READS)])
        raw_outputs.append(output)
        df = pd.read_csv(StringIO(output))
        df["repeat"] = repeat
        frames.append(df)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    all_df = pd.concat(frames, ignore_index=True)
    all_df.to_csv(OUT_DIR / "native_edlib_assignment.csv", index=False)
    all_df.to_csv(RAW_DIR / "native_edlib_assignment.csv", index=False)
    (RAW_DIR / "native_edlib_assignment_raw_runs.csv").write_text("\n".join(raw_outputs), encoding="utf-8")
    df = all_df
    return df


def speedup_frame(df: pd.DataFrame) -> pd.DataFrame:
    keys = ["workload", "error_mode", "n_reads", "n_targets", "len", "k", "err", "indel_rate", "repeat"]
    q = df[df["tool"] == "dotmatch_indexed"][keys + ["reads_per_sec", "candidates_per_read", "verified_per_read", "peak_rss_kb", "mismatches"]]
    e = df[df["tool"] == "edlib_native_scan"][keys + ["reads_per_sec"]]
    merged = q.merge(e, on=keys, suffixes=("_dotmatch", "_edlib"))
    merged["speedup_vs_edlib_native"] = merged["reads_per_sec_dotmatch"] / merged["reads_per_sec_edlib"]
    return merged


def aggregate_stats(df: pd.DataFrame) -> pd.DataFrame:
    keys = ["tool", "workload", "error_mode", "n_targets", "len", "k", "err", "indel_rate"]
    grouped = df.groupby(keys)
    out = grouped.agg(
        repeats=("repeat", "nunique"),
        reads_per_sec_mean=("reads_per_sec", "mean"),
        reads_per_sec_p50=("reads_per_sec", "median"),
        reads_per_sec_p95=("reads_per_sec", lambda s: s.quantile(0.95)),
        reads_per_sec_std=("reads_per_sec", "std"),
        seconds_mean=("seconds", "mean"),
        verified_per_read_median=("verified_per_read", "median"),
        peak_rss_kb_max=("peak_rss_kb", "max"),
        mismatches_sum=("mismatches", "sum"),
    ).reset_index()
    out["reads_per_sec_cv"] = out["reads_per_sec_std"].fillna(0.0) / out["reads_per_sec_mean"]
    return out


def plot_speedup(speedups: pd.DataFrame) -> None:
    subset = speedups[(speedups["len"].isin([16, 32])) & (speedups["err"].isin([0.0, 0.01])) & (speedups["error_mode"].isin(["exact", "one_substitution"]))]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for (length, k), group in subset.groupby(["len", "k"]):
        summary = group.groupby("n_targets")["speedup_vs_edlib_native"].median().reset_index()
        ax.plot(summary["n_targets"], summary["speedup_vs_edlib_native"], marker="o", label=f"len={length} k={k}")
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


def plot_candidates(speedups: pd.DataFrame) -> None:
    subset = speedups[(speedups["len"] == 32) & (speedups["err"] == 0.01) & (speedups["error_mode"] == "one_substitution")]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for k, group in subset.groupby("k"):
        summary = group.groupby("n_targets")["verified_per_read"].median().reset_index()
        ax.plot(summary["n_targets"], summary["verified_per_read"], marker="o", label=f"k={k}")
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


def plot_throughput(df: pd.DataFrame) -> None:
    subset = df[(df["len"] == 32) & (df["k"] == 1) & (df["err"] == 0.01) & (df["error_mode"] == "one_substitution")]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for tool, group in subset.groupby("tool"):
        summary = group.groupby("n_targets")["reads_per_sec"].median().reset_index()
        ax.plot(summary["n_targets"], summary["reads_per_sec"], marker="o", label=tool)
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


def write_report(df: pd.DataFrame, speedups: pd.DataFrame) -> None:
    agg = aggregate_stats(df)
    agg.to_csv(RAW_DIR / "native_edlib_assignment_summary.csv", index=False)
    best = speedups.groupby(["n_targets", "len", "k", "err", "error_mode"]).agg(
        reads_per_sec_dotmatch=("reads_per_sec_dotmatch", "median"),
        reads_per_sec_edlib=("reads_per_sec_edlib", "median"),
        verified_per_read=("verified_per_read", "median"),
        peak_rss_kb=("peak_rss_kb", "max"),
        mismatches=("mismatches", "sum"),
        speedup_vs_edlib_native=("speedup_vs_edlib_native", "median"),
    ).reset_index().sort_values("speedup_vs_edlib_native", ascending=False).head(12)
    summary = speedups.groupby(["len", "k", "n_targets", "error_mode"])["speedup_vs_edlib_native"].median().reset_index()
    summary = summary.sort_values("speedup_vs_edlib_native", ascending=False).head(12)
    zero_mismatch = int(df["mismatches"].sum())
    lines = [
        "# Native Edlib Benchmark Report",
        "",
        f"- Platform: `{platform.platform()}`",
        f"- Python: `{platform.python_version()}`",
        f"- Reads per benchmark case: `{REPORT_READS}`",
        f"- Repetitions per benchmark case: `{REPORT_REPEATS}`",
        "- Comparator: native Edlib C/C++ API, `EDLIB_MODE_NW`, `EDLIB_TASK_DISTANCE`, fixed threshold `k`.",
        "- Additional baselines: exact hash lookup for `k=0`; BK-tree and neighbor lookup approximate baselines for `k=1`.",
        f"- Assignment mismatches recorded across all rows: `{zero_mismatch}`.",
        "- Every benchmark run aborts on assignment disagreement between DotMatch and native Edlib scan.",
        "",
        "![Native speedup vs Edlib](native_speedup_vs_edlib.svg)",
        "",
        "![Native candidates per read](native_candidates_per_read.svg)",
        "",
        "![Native assignment throughput](native_assignment_throughput.svg)",
        "",
        "## Highest Observed Microbenchmark Speedups",
        "",
        markdown_table(best[["n_targets", "len", "k", "error_mode", "err", "reads_per_sec_dotmatch", "reads_per_sec_edlib", "verified_per_read", "peak_rss_kb", "speedup_vs_edlib_native"]], ".2f"),
        "",
        "## Median Speedup Summary",
        "",
        markdown_table(summary, ".2f"),
        "",
        "## Repeated-Run Statistics",
        "",
        markdown_table(agg[agg["tool"].isin(["dotmatch_indexed", "edlib_native_scan"])].head(16)[["tool", "error_mode", "n_targets", "len", "k", "reads_per_sec_mean", "reads_per_sec_p50", "reads_per_sec_p95", "reads_per_sec_cv", "peak_rss_kb_max", "mismatches_sum"]], ".2f"),
        "",
        "## Evidence Boundary",
        "",
        "These are native Edlib scan microbenchmarks for exact short-DNA assignment workloads, plus simple exact-hash and BK-tree baselines. The largest rows are useful for understanding algorithmic scaling against exhaustive scan, but they are not end-to-end workflow speed claims. Exact `k=0` lookup should be judged against hash-table baselines. For `k=1`, the indexed path is reported only when it has zero correctness disagreements against the exhaustive comparator. Levenshtein `k=2` uses the exhaustive assignment path.",
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
