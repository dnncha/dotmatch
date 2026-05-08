#!/usr/bin/env python3
"""Generate benchmark tables and README-ready graphs.

This script runs the current reproducible benchmarks, parses their CSV-shaped
output, and writes a small report under docs/benchmarks/.
"""

from __future__ import annotations

import csv
import platform
import subprocess
from io import StringIO
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "docs" / "benchmarks"
REPORT_BATCH_READS = 200


def run_command(cmd: list[str]) -> str:
    result = subprocess.run(cmd, cwd=ROOT, check=True, text=True, capture_output=True)
    return result.stdout


def parse_edlib_output(text: str) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    reader = csv.DictReader(
        line for line in StringIO(text) if line.strip() and not line.startswith("#")
    )
    for row in reader:
        rows.append(row)
    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("no Edlib benchmark rows parsed")
    df["len"] = df["len"].astype(int)
    df["k"] = pd.to_numeric(df["k"], errors="coerce")
    df["err"] = df["err"].astype(float)
    df["calls"] = df["calls"].astype(int)
    df["seconds"] = df["seconds"].astype(float)
    df["calls_per_sec"] = df["calls_per_sec"].astype(float)
    df["ns_per_call"] = df["ns_per_call"].astype(float)
    return df


def parse_batch_output(text: str) -> pd.DataFrame:
    df = pd.read_csv(StringIO(text))
    for col in ["n_reads", "n_targets", "len", "k"]:
        df[col] = df[col].astype(int)
    for col in ["err", "seconds", "reads_per_sec", "assignments_per_sec"]:
        df[col] = df[col].astype(float)
    return df


def speedup_frame(df: pd.DataFrame, kind: str) -> pd.DataFrame:
    part = df[df["kind"] == kind]
    key = ["kind", "len", "k", "err"] if kind == "threshold" else ["kind", "len", "err"]
    qda = part[part["tool"] == "dotmatch"][key + ["ns_per_call"]].rename(
        columns={"ns_per_call": "dotmatch_ns"}
    )
    edlib = part[part["tool"] == "edlib"][key + ["ns_per_call"]].rename(
        columns={"ns_per_call": "edlib_ns"}
    )
    merged = qda.merge(edlib, on=key)
    merged["speedup_vs_edlib"] = merged["edlib_ns"] / merged["dotmatch_ns"]
    return merged


def plot_exact_speedup(exact: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for err, group in exact.groupby("err"):
        ax.plot(group["len"], group["speedup_vs_edlib"], marker="o", label=f"{err:.0%} errors")
    ax.axhline(1.0, color="#444444", linewidth=1, linestyle="--")
    ax.set_title("Exact global edit distance vs Edlib Python binding")
    ax.set_xlabel("sequence length (bp)")
    ax.set_ylabel("speedup vs Edlib")
    ax.set_xticks(sorted(exact["len"].unique()))
    ax.legend(frameon=False, ncols=2)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "exact_speedup_vs_edlib.svg")
    plt.close(fig)


def plot_threshold_speedup(threshold: pd.DataFrame) -> None:
    pivot = threshold.pivot_table(
        index="len", columns="k", values="speedup_vs_edlib", aggfunc="mean"
    )
    fig, ax = plt.subplots(figsize=(7, 4.5))
    image = ax.imshow(pivot.values, cmap="viridis", aspect="auto")
    ax.set_title("Mean threshold speedup vs Edlib Python binding")
    ax.set_xlabel("threshold k")
    ax.set_ylabel("sequence length (bp)")
    ax.set_xticks(range(len(pivot.columns)), [str(int(k)) for k in pivot.columns])
    ax.set_yticks(range(len(pivot.index)), [str(int(v)) for v in pivot.index])
    for y in range(pivot.shape[0]):
        for x in range(pivot.shape[1]):
            ax.text(x, y, f"{pivot.values[y, x]:.2f}x", ha="center", va="center", color="white")
    fig.colorbar(image, ax=ax, label="speedup")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "threshold_speedup_heatmap.svg")
    plt.close(fig)


def plot_batch_speed(batch: pd.DataFrame) -> None:
    subset = batch[
        (batch["n_reads"] == batch["n_reads"].max())
        & (batch["len"].isin([16, 32]))
        & (batch["k"].isin([1, 2]))
        & (batch["err"].isin([0.0, 0.01]))
    ].copy()
    subset["case"] = (
        "len="
        + subset["len"].astype(str)
        + " k="
        + subset["k"].astype(str)
        + " err="
        + (subset["err"] * 100).round(1).astype(str)
        + "%"
        + " targets="
        + subset["n_targets"].astype(str)
    )
    subset = subset.sort_values(["n_targets", "len", "k", "err", "tool"])

    fig, ax = plt.subplots(figsize=(10, 5))
    tools = ["dotmatch_indexed", "dotmatch_scan", "dotmatch_naive"]
    width = 0.38
    cases = list(dict.fromkeys(subset["case"]))
    x_positions = range(len(cases))
    for offset, tool in enumerate(tools):
        values = []
        for case in cases:
            row = subset[(subset["case"] == case) & (subset["tool"] == tool)]
            values.append(float(row["reads_per_sec"].iloc[0]) if not row.empty else 0.0)
        shifted = [x + (offset - 0.5) * width for x in x_positions]
        ax.bar(shifted, values, width=width, label=tool)

    ax.set_title("Synthetic barcode assignment throughput")
    ax.set_ylabel("reads per second")
    ax.set_xticks(list(x_positions), cases, rotation=45, ha="right")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "batch_assignment_throughput.svg")
    plt.close(fig)


def markdown_table(df: pd.DataFrame, floatfmt: str = ".2f") -> str:
    columns = list(df.columns)
    rows = []
    for _, row in df.iterrows():
        cells = []
        for col in columns:
            value = row[col]
            if pd.isna(value):
                cells.append("")
            elif col in {"len", "k", "n_targets", "n_reads", "calls"}:
                cells.append(str(int(value)))
            elif isinstance(value, float):
                cells.append(format(value, floatfmt))
            else:
                cells.append(str(value))
        rows.append(cells)

    out = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        out.append("| " + " | ".join(row) + " |")
    return "\n".join(out)


def write_report(edlib_df: pd.DataFrame, batch_df: pd.DataFrame) -> None:
    exact = speedup_frame(edlib_df, "exact")
    threshold = speedup_frame(edlib_df, "threshold")

    exact_summary = exact.groupby("len")["speedup_vs_edlib"].mean().reset_index()
    threshold_summary = threshold.groupby(["len", "k"])["speedup_vs_edlib"].mean().reset_index()

    best_threshold = threshold.sort_values("speedup_vs_edlib", ascending=False).head(8)
    batch_summary = (
        batch_df.groupby(["tool", "len", "k", "n_targets"])["reads_per_sec"]
        .median()
        .reset_index()
        .sort_values("reads_per_sec", ascending=False)
        .head(12)
    )

    lines = [
        "# Benchmark Report",
        "",
        f"- Platform: `{platform.platform()}`",
        f"- Python: `{platform.python_version()}`",
        "- External exact edit-distance baseline: Edlib Python binding with `mode=\"NW\"`, `task=\"distance\"`.",
        "- These graphs compare short global edit-distance workloads. They are not broad aligner benchmarks.",
        "",
        "## Graphs",
        "",
        "![Exact speedup vs Edlib](exact_speedup_vs_edlib.svg)",
        "",
        "![Threshold speedup heatmap](threshold_speedup_heatmap.svg)",
        "",
        "![Batch assignment throughput](batch_assignment_throughput.svg)",
        "",
        "## Exact Distance Mean Speedup",
        "",
        markdown_table(exact_summary, ".2f"),
        "",
        "## Threshold Best Cases",
        "",
        markdown_table(best_threshold[["len", "k", "err", "dotmatch_ns", "edlib_ns", "speedup_vs_edlib"]], ".2f"),
        "",
        "## Batch Assignment Top Throughput Rows",
        "",
        markdown_table(batch_summary, ".1f"),
        "",
        "## Evidence Boundary",
        "",
        "These results cover short-DNA global edit-distance and threshold matching against Edlib's Python binding.",
        "Comparative performance wording should use the native Edlib assignment report for checked native-library claims. SeqAn/Parasail comparisons are not claimed until docs/native-comparator-scope.md requirements are met.",
        "",
    ]
    (OUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    run_command(["make", "shared"])
    run_command(["make", "build/bench_batch"])

    edlib_output = run_command(["python3", "scripts/bench_vs_edlib.py"])
    batch_output = run_command(["./build/bench_batch", str(REPORT_BATCH_READS)])

    (OUT_DIR / "edlib_python.csv").write_text(
        "\n".join(line for line in edlib_output.splitlines() if line and not line.startswith("#")) + "\n",
        encoding="utf-8",
    )
    (OUT_DIR / "batch_assignment.csv").write_text(batch_output, encoding="utf-8")

    edlib_df = parse_edlib_output(edlib_output)
    batch_df = parse_batch_output(batch_output)

    exact = speedup_frame(edlib_df, "exact")
    threshold = speedup_frame(edlib_df, "threshold")
    plot_exact_speedup(exact)
    plot_threshold_speedup(threshold)
    plot_batch_speed(batch_df)
    write_report(edlib_df, batch_df)
    print(f"wrote benchmark report to {OUT_DIR}")


if __name__ == "__main__":
    main()
