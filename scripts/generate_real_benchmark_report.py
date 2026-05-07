#!/usr/bin/env python3
"""Generate real public CRISPR benchmark data and a README-ready report."""

from __future__ import annotations

import os
import shutil
import subprocess
from io import StringIO
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "benchmarks" / "real" / "data"
RAW = ROOT / "benchmarks" / "raw"
FIG = ROOT / "benchmarks" / "figures"
DOC = ROOT / "docs" / "benchmarks" / "real"
READ_LIMIT = int(os.environ.get("DOTMATCH_REAL_READS", "25"))
FETCH_RECORDS = int(os.environ.get("DOTMATCH_REAL_FETCH_RECORDS", str(max(READ_LIMIT, 25))))


def run(cmd: list[str], capture: bool = True, check: bool = True) -> str:
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=capture, check=check)
    return result.stdout if capture else ""


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
            if isinstance(value, float):
                vals.append(format(value, floatfmt))
            else:
                vals.append(str(value))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def plot_real(df: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.bar(df["tool"], df["reads_per_sec"], color=["#175a7a", "#8b3a3a"])
    ax.set_yscale("log")
    ax.set_ylabel("reads/sec")
    ax.set_title("Real public CRISPR reads: k=1 assignment")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "real_crispr_throughput.svg")
    fig.savefig(FIG / "real_crispr_throughput.pdf")
    plt.close(fig)


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    DOC.mkdir(parents=True, exist_ok=True)
    run(["python3", "scripts/fetch_mageck_demo.py", "--subsample", str(FETCH_RECORDS), "--out", str(DATA)])
    run(["make", "build/bench_real_edlib"], capture=False)
    read1 = DATA / "ERR376998.fastq.gz"
    read2 = DATA / "ERR376999.fastq.gz"
    targets = DATA / "yusa_library.csv"
    output = run(
        [
            "./build/bench_real_edlib",
            "--targets",
            str(targets),
            "--reads",
            str(read1),
            "--reads",
            str(read2),
            "--target-start",
            "23",
            "--target-length",
            "19",
            "--k",
            "1",
            "--limit",
            str(READ_LIMIT),
        ]
    )
    (RAW / "real_crispr_edlib.csv").write_text(output, encoding="utf-8")
    df = pd.read_csv(StringIO(output))
    plot_real(df)

    comparator_rows = []
    for tool in ["mageck", "cutadapt", "bowtie2"]:
        comparator_rows.append(
            {
                "tool": tool,
                "status": "available" if shutil.which(tool) else "not installed in this environment",
            }
        )
    comparators = pd.DataFrame(comparator_rows)
    comparators.to_csv(RAW / "real_workflow_comparator_availability.csv", index=False)

    speedup = None
    if {"dotmatch_indexed", "edlib_native_scan"}.issubset(set(df["tool"])):
        q = float(df[df["tool"] == "dotmatch_indexed"]["reads_per_sec"].iloc[0])
        e = float(df[df["tool"] == "edlib_native_scan"]["reads_per_sec"].iloc[0])
        speedup = q / e if e > 0 else None

    lines = [
        "# Real CRISPR Benchmark",
        "",
        "This report uses real public MAGeCK/Yusa CRISPR screen inputs, not synthetic reads.",
        "",
        f"- Guide library: `yusa_library.csv`",
        f"- FASTQ samples: `ERR376998.fastq.gz`, `ERR376999.fastq.gz`",
        "- Extraction: `--target-start 23 --target-length 19`",
        "- Assignment threshold: exact Levenshtein `k=1`",
        f"- Real reads benchmarked against Edlib scan: `{READ_LIMIT}` extracted reads",
        f"- Real FASTQ records downloaded per sample: `{FETCH_RECORDS}`",
        "- Correctness comparator: native Edlib exhaustive scan over the guide library",
        "",
        "![Real CRISPR throughput](../../../benchmarks/figures/real_crispr_throughput.svg)",
        "",
        "## Native Real-Data Results",
        "",
        markdown_table(df, ".2f"),
        "",
        "## Speedup",
        "",
        f"DotMatch indexed speedup vs native Edlib scan on this real-data subset: `{speedup:.2f}x`." if speedup is not None else "Speedup unavailable.",
        "",
        "## Workflow Comparator Availability",
        "",
        markdown_table(comparators),
        "",
        "MAGeCK, Cutadapt, and Bowtie2 are workflow comparators and are intentionally optional. Run `python3 scripts/run_public_crispr_benchmark.py --run-mageck --run-cutadapt --run-bowtie2` in an environment with those tools installed to populate those rows.",
        "",
        "## Claim Discipline",
        "",
        "This benchmark supports the real-workflow claim only for known-target CRISPR guide assignment under the listed extraction rules. Native Edlib remains the exact semantic oracle; MAGeCK/Cutadapt/Bowtie2 comparisons should be described as workflow comparisons, not identical semantic oracles.",
        "",
    ]
    (DOC / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print(DOC / "README.md")


if __name__ == "__main__":
    main()
