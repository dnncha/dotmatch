#!/usr/bin/env python3
"""Generate the experimental GPU acceleration benchmark report."""

from __future__ import annotations

import csv
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "gpu_acceleration.csv"
OUT_DIR = ROOT / "docs" / "benchmarks" / "gpu"
FIG_DIR = ROOT / "benchmarks" / "figures"


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def fnum(value: str | None) -> float:
    try:
        return float(value or 0.0)
    except ValueError:
        return 0.0


def case_key(row: dict[str, str]) -> tuple[str, str, str, str, str]:
    return (
        row.get("n_reads", ""),
        row.get("n_targets", ""),
        row.get("len", ""),
        row.get("k", ""),
        row.get("error_rate", ""),
    )


def speed_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    cpu = {
        case_key(row): row
        for row in rows
        if row.get("tool") == "dotmatch_cpu_index" and row.get("status") == "ok"
    }
    out: list[dict[str, str]] = []
    for row in rows:
        if row.get("tool") != "dotmatch_gpu_metal" or row.get("status") != "ok":
            continue
        base = cpu.get(case_key(row))
        if base is None:
            continue
        gpu_rps = fnum(row.get("reads_per_sec"))
        cpu_rps = fnum(base.get("reads_per_sec"))
        gpu_total_rps = fnum(row.get("total_reads_per_sec"))
        cpu_total_rps = fnum(base.get("total_reads_per_sec"))
        out.append({
            "n_reads": row.get("n_reads", ""),
            "n_targets": row.get("n_targets", ""),
            "len": row.get("len", ""),
            "k": row.get("k", ""),
            "gpu_reads_per_sec": f"{gpu_rps:.1f}",
            "cpu_reads_per_sec": f"{cpu_rps:.1f}",
            "gpu_vs_cpu_kernel": f"{(gpu_rps / cpu_rps):.2f}x" if cpu_rps else "",
            "gpu_vs_cpu_total": f"{(gpu_total_rps / cpu_total_rps):.2f}x" if cpu_total_rps else "",
            "mismatches": row.get("mismatches", ""),
        })
    return out


def markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    if not rows:
        return "_No rows available._\n"
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines) + "\n"


def svg_speedup(rows: list[dict[str, str]], path: Path) -> None:
    selected = speed_rows(rows)
    if not selected:
        return
    labels = [f"{r['n_reads']} reads / {r['n_targets']} targets" for r in selected]
    values = []
    for row in selected:
        text = row["gpu_vs_cpu_total"].replace("x", "")
        values.append(fnum(text))
    width = 900
    left = 320
    row_h = 34
    height = 76 + row_h * len(selected)
    max_v = max(max(values), 1.0)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,sans-serif;font-size:13px}.title{font-size:18px;font-weight:700}.axis{fill:#444}.bar{fill:#476fb3}</style>',
        '<text class="title" x="20" y="28">Experimental Metal GPU total throughput ratio</text>',
        '<text class="axis" x="20" y="50">GPU total reads/sec divided by CPU indexed total reads/sec; higher is better</text>',
    ]
    for i, (label, value) in enumerate(zip(labels, values)):
        y = 78 + i * row_h
        w = max(1, int((width - left - 110) * value / max_v))
        parts.append(f'<text x="20" y="{y + 15}">{label}</text>')
        parts.append(f'<rect class="bar" x="{left}" y="{y}" width="{w}" height="21" rx="2"/>')
        parts.append(f'<text x="{left + w + 8}" y="{y + 15}">{value:.2f}x</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    rows = read_rows(RAW)
    if not rows:
        raise SystemExit(f"missing GPU benchmark CSV: {RAW}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    speedup_svg = FIG_DIR / "gpu_metal_speedup.svg"
    svg_speedup(rows, speedup_svg)

    unavailable = [row for row in rows if row.get("status") == "unavailable"]
    lines = [
        "# Experimental GPU Acceleration Benchmark",
        "",
        "This report is a skunk-works GPU evidence lane. It is intentionally not a production speed claim: the Metal path brute-forces packed Hamming `k=1` distances and is compared against DotMatch's existing CPU indexed Hamming assignment with identical output checks.",
        "",
        "The decision rule is simple: GPU rows must have zero mismatches before any speed result is considered, and CPU-indexed throughput remains the production baseline unless the GPU path wins end-to-end on real workloads.",
        "",
    ]
    if unavailable:
        lines.extend([
            "## Availability",
            "",
            markdown_table(unavailable, ["tool", "backend", "status", "device", "notes"]),
            "",
        ])
    if speedup_svg.exists():
        lines.extend([
            "## Figure",
            "",
            f"![GPU speedup]({os.path.relpath(speedup_svg, OUT_DIR)})",
            "",
        ])
    lines.extend([
        "## CPU vs Metal Rows",
        "",
        markdown_table(speed_rows(rows), [
            "n_reads", "n_targets", "len", "k", "gpu_reads_per_sec", "cpu_reads_per_sec",
            "gpu_vs_cpu_kernel", "gpu_vs_cpu_total", "mismatches",
        ]),
        "",
        "## Raw Rows",
        "",
        markdown_table(rows, [
            "tool", "backend", "status", "workload", "n_reads", "n_targets", "len", "k",
            "prep_seconds", "seconds", "total_seconds", "reads_per_sec",
            "total_reads_per_sec", "pairs_per_sec", "checksum", "mismatches", "device", "notes",
        ]),
        "",
        "## Scope",
        "",
        "This lane tests whether GPU compute is worth productizing. It currently covers fixed-length packed A/C/G/T Hamming `k=1` assignment only. It does not cover Levenshtein indels, BCL conversion, GPU-resident FASTQ parsing, CUDA deployment, or production scheduling.",
        "",
        "A future production GPU path needs a real workload gate that includes sequence extraction, packing, transfer or shared-memory preparation, dispatch, readback, and downstream count/QC generation.",
        "",
    ])
    (OUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print(OUT_DIR / "README.md")


if __name__ == "__main__":
    main()
