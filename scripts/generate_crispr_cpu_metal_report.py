#!/usr/bin/env python3
"""Generate the production crispr-count CPU vs Metal benchmark report."""

from __future__ import annotations

import csv
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "crispr_cpu_metal.csv"
OUT = ROOT / "docs" / "benchmarks" / "gpu" / "production_crispr_cpu_metal.md"
FIG = ROOT / "benchmarks" / "figures" / "crispr_cpu_metal_speedup.svg"


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


def markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    if not rows:
        return "_No rows available._\n"
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines) + "\n"


def speedup_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    cpu = {
        (row["dataset"], row["records_per_sample"]): row
        for row in rows
        if row.get("backend") == "cpu" and row.get("status") == "ok"
    }
    out: list[dict[str, str]] = []
    for row in rows:
        if row.get("backend") != "gpu-metal-experimental" or row.get("status") != "ok":
            continue
        base = cpu.get((row["dataset"], row["records_per_sample"]))
        if base is None:
            continue
        metal_rps = fnum(row.get("reads_per_sec"))
        cpu_rps = fnum(base.get("reads_per_sec"))
        metal_wall = fnum(row.get("wall_seconds"))
        cpu_wall = fnum(base.get("wall_seconds"))
        out.append({
            "dataset": row.get("dataset", ""),
            "records_per_sample": row.get("records_per_sample", ""),
            "metal_validate": row.get("metal_validate", ""),
            "cpu_reads_per_sec": f"{cpu_rps:.1f}",
            "metal_reads_per_sec": f"{metal_rps:.1f}",
            "metal_vs_cpu_wall": f"{(cpu_wall / metal_wall):.2f}x" if metal_wall else "",
            "metal_vs_cpu_throughput": f"{(metal_rps / cpu_rps):.2f}x" if cpu_rps else "",
            "count_match_cpu": row.get("count_match_cpu", ""),
            "metal_validation": row.get("metal_validation", ""),
            "count_engine": row.get("count_engine", ""),
        })
    return out


def svg_speedup(rows: list[dict[str, str]]) -> None:
    selected = [row for row in speedup_rows(rows) if row.get("metal_validate") == "0" and row.get("dataset") == "mageck_yusa"]
    if not selected:
        return
    labels = [f"Yusa {row['records_per_sample']} reads/sample" for row in selected]
    values = [fnum(row["metal_vs_cpu_wall"].replace("x", "")) for row in selected]
    width = 900
    left = 300
    row_h = 34
    height = 76 + row_h * len(selected)
    max_v = max(max(values, default=1.0), 1.0)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,sans-serif;font-size:13px}.title{font-size:18px;font-weight:700}.axis{fill:#444}.bar{fill:#8c4b2f}</style>',
        '<text class="title" x="20" y="28">Production crispr-count CPU wall time / Metal wall time</text>',
        '<text class="axis" x="20" y="50">Fixed-offset MAGeCK/Yusa Hamming k=1; higher means Metal is faster</text>',
    ]
    for i, (label, value) in enumerate(zip(labels, values)):
        y = 78 + i * row_h
        w = max(1, int((width - left - 110) * value / max_v))
        parts.append(f'<text x="20" y="{y + 15}">{label}</text>')
        parts.append(f'<rect class="bar" x="{left}" y="{y}" width="{w}" height="21" rx="2"/>')
        parts.append(f'<text x="{left + w + 8}" y="{y + 15}">{value:.2f}x</text>')
    parts.append("</svg>")
    FIG.parent.mkdir(parents=True, exist_ok=True)
    FIG.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    rows = read_rows(RAW)
    if not rows:
        raise SystemExit(f"missing benchmark CSV: {RAW}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    svg_speedup(rows)
    speedups = speedup_rows(rows)
    ineligible = [
        row for row in rows
        if row.get("status") in {"ineligible", "unavailable", "validation_failed", "error"}
        and (row.get("backend") == "gpu-metal-experimental" or row.get("status") != "error")
    ]
    lines = [
        "# Production CRISPR CPU vs Metal",
        "",
        "This report benchmarks the shipping `dotmatch crispr-count` CLI path, not the standalone `bench_gpu_crispr_metal` microbench. CPU remains the assignment authority. Metal is opt-in via `--backend gpu-metal-experimental` and should be paired with `--metal-validate` before production use.",
        "",
        "## Decision Rule",
        "",
        "- Use CPU when the workload is not Metal-eligible, when `auto` is selected, or when validation fails.",
        "- Treat any Metal speedup as advisory until `metal_validation=passed` and guide-by-guide counts match the CPU shadow run.",
        "- Do not use Metal for Sanson/Brunello-style multi-offset counting; the production path rejects or blocks that lane today.",
        "",
    ]
    if FIG.exists():
        lines.extend([
            "## Yusa Speedup Figure",
            "",
            f"![Production CRISPR CPU vs Metal]({os.path.relpath(FIG, OUT.parent)})",
            "",
        ])
    lines.extend([
        "## CPU vs Metal Speedup",
        "",
        markdown_table(speedups, [
            "dataset", "records_per_sample", "metal_validate", "cpu_reads_per_sec",
            "metal_reads_per_sec", "metal_vs_cpu_wall", "metal_vs_cpu_throughput",
            "count_match_cpu", "metal_validation", "count_engine",
        ]),
        "",
        "## Raw Rows",
        "",
        markdown_table(rows, [
            "dataset", "backend", "metal_validate", "status", "records_per_sample",
            "total_reads", "n_targets", "offset_mode", "auto_offset", "wall_seconds",
            "reads_per_sec", "backend_effective", "count_engine", "metal_validation",
            "count_match_cpu", "exit_code", "device", "notes",
        ]),
        "",
    ])
    if ineligible:
        lines.extend([
            "## Ineligible Or Failed Rows",
            "",
            markdown_table(ineligible, [
                "dataset", "backend", "metal_validate", "status", "records_per_sample",
                "offset_mode", "auto_offset", "notes",
            ]),
            "",
        ])
    yusa_metal = [
        row for row in rows
        if row.get("dataset") == "mageck_yusa" and row.get("backend") == "gpu-metal-experimental" and row.get("metal_validate") == "0"
    ]
    yusa_validate = [row for row in rows if row.get("status") == "validation_failed"]
    lines.extend([
        "## Interpretation",
        "",
    ])
    if yusa_metal and any(row.get("count_match_cpu") == "0" for row in yusa_metal):
        lines.append(
            "- **Yusa Metal without `--metal-validate` can be faster but is not count-identical to CPU** on the recorded rows (`count_match_cpu=0`). Do not use those counts for downstream analysis."
        )
    if yusa_validate:
        lines.append(
            "- **`--metal-validate` fails on Yusa today** because the experimental Metal path disagrees with the CPU authority checksum. The CLI exits non-zero and deletes trust in Metal counts until the mismatch is fixed."
        )
    cpu_100k = next((row for row in rows if row.get("dataset") == "mageck_yusa" and row.get("backend") == "cpu" and row.get("records_per_sample") == "100000"), None)
    metal_100k = next((row for row in yusa_metal if row.get("records_per_sample") == "100000"), None)
    if cpu_100k and metal_100k and fnum(cpu_100k.get("wall_seconds")) < fnum(metal_100k.get("wall_seconds")):
        lines.append(
            "- **At 100k reads/sample on this host, CPU is faster than Metal even before validation overhead.** Production `auto` staying on CPU is consistent with these measurements."
        )
    lines.append(
        "- **Sanson/Brunello remains CPU-only** because `--offset-mode multi` is outside the production Metal eligibility contract."
    )
    lines.extend(["", "## Workload Notes",
        "",
        "### MAGeCK/Yusa (`mageck_yusa`)",
        "",
        "- Fixed guide window: `--guide-start 23 --guide-length 19`.",
        "- Hamming `k=1` with `--ambiguity-policy best`.",
        "- No auto-offset, so the workload stays in the Metal-eligible single-offset lane.",
        "- `auto` intentionally stays on CPU in the shipping CLI; use explicit `--backend gpu-metal-experimental` to test Metal.",
        "",
        "### Sanson/Brunello (`sanson_brunello`)",
        "",
        "- Uses `--offset-mode multi` and auto-offset detection, matching the public guide-counter lane.",
        "- That multi-offset path is CPU-only today; Metal rows are recorded as `ineligible`.",
        "",
        "## Commands",
        "",
        "```bash",
        "make bench-crispr-cpu-metal",
        "make crispr-cpu-metal-report",
        "```",
        "",
        "Optional full FASTQs:",
        "",
        "```bash",
        "DOTMATCH_CPU_METAL_FULL=1 make bench-crispr-cpu-metal",
        "```",
        "",
    ])
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()