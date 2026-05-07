#!/usr/bin/env python3
"""Generate an honest atlas CRISPR progress report from the latest raw rows."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "atlas"
FIG = ROOT / "benchmarks" / "figures"
OUT = ROOT / "docs" / "benchmarks" / "crispr_sota" / "atlas_progress.md"


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh))


def fnum(row: dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, "") or 0.0)
    except ValueError:
        return 0.0


def label(row: dict[str, str]) -> str:
    dataset = row.get("dataset_id", "")
    records = row.get("requested_records_per_sample") or ("full" if row.get("n_reads") == "246950411" else "")
    return f"{dataset} {records} {row.get('tool', '')}"


def svg_bar(rows: list[dict[str, str]], key: str, title: str, path: Path, lower_is_better: bool = False) -> None:
    selected = [r for r in rows if fnum(r, key) > 0]
    if not selected:
        return
    width = 1180
    left = 470
    row_h = 30
    height = 72 + row_h * len(selected)
    max_v = max(fnum(r, key) for r in selected) or 1.0
    color = "#7a4f9f" if lower_is_better else "#2f7d68"
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,sans-serif;font-size:12px}.title{font-size:18px;font-weight:700}.sub{fill:#444}.bar{fill:' + color + '}</style>',
        f'<text class="title" x="20" y="28">{title}</text>',
        f'<text class="sub" x="20" y="50">{key}; {"lower is better" if lower_is_better else "higher is better"}</text>',
    ]
    for i, row in enumerate(selected):
        y = 76 + i * row_h
        value = fnum(row, key)
        w = max(1, int((width - left - 140) * value / max_v))
        parts.append(f'<text x="20" y="{y + 14}">{label(row)}</text>')
        parts.append(f'<rect class="bar" x="{left}" y="{y}" width="{w}" height="19" rx="2"/>')
        parts.append(f'<text x="{left + w + 8}" y="{y + 14}">{value:.2f}</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def table(rows: list[dict[str, str]]) -> str:
    cols = ["dataset_id", "run_level", "requested_records_per_sample", "tool", "n_reads", "seconds",
            "reads_per_sec", "peak_rss_kb", "assigned_reads", "corrected_reads", "ambiguous_reads",
            "verified_per_read"]
    lines = ["|" + "|".join(cols) + "|", "|" + "|".join(["---"] * len(cols)) + "|"]
    for row in rows:
        lines.append("|" + "|".join(row.get(c, "") for c in cols) + "|")
    return "\n".join(lines)


def row_by(rows: list[dict[str, str]], dataset: str, tool: str, run_level: str | None = None) -> dict[str, str] | None:
    for row in rows:
        if row.get("dataset_id") != dataset or row.get("tool") != tool:
            continue
        if run_level is not None and row.get("run_level") != run_level:
            continue
        return row
    return None


def speedup(a: dict[str, str] | None, b: dict[str, str] | None) -> str:
    if not a or not b or fnum(a, "seconds") <= 0:
        return "n/a"
    return f"{fnum(b, 'seconds') / fnum(a, 'seconds'):.2f}x"


def main() -> None:
    latest_100k = read_rows(RAW / "crispr_sota_after_reserved_slab_100k.csv")
    full_all = read_rows(RAW / "crispr_sota_full_hamming_after_rolling.csv")
    sanson_full_dotmatch = read_rows(RAW / "sanson_full_compact_dotmatch_only.csv")

    rows: list[dict[str, str]] = []
    rows.extend(latest_100k)
    full_sanson_competitors = [
        r for r in full_all
        if r.get("dataset_id") == "sanson_brunello"
        and r.get("run_level") == "full"
        and r.get("tool") in {"mageck_count_exact", "guide_counter_one_mismatch"}
    ]
    for row in sanson_full_dotmatch:
        row = dict(row)
        row.setdefault("run_level", "full")
        row.setdefault("requested_records_per_sample", "full")
        rows.append(row)
    rows.extend(full_sanson_competitors)

    FIG.mkdir(parents=True, exist_ok=True)
    runtime_fig = FIG / "crispr_atlas_progress_runtime.svg"
    throughput_fig = FIG / "crispr_atlas_progress_throughput.svg"
    memory_fig = FIG / "crispr_atlas_progress_memory.svg"
    svg_bar(rows, "seconds", "Atlas CRISPR real-data runtime", runtime_fig, lower_is_better=True)
    svg_bar(rows, "reads_per_sec", "Atlas CRISPR real-data throughput", throughput_fig)
    svg_bar(rows, "peak_rss_kb", "Atlas CRISPR peak RSS", memory_fig, lower_is_better=True)

    yusa_dm = row_by(latest_100k, "mageck_yusa", "dotmatch_hamming_k1")
    yusa_gc = row_by(latest_100k, "mageck_yusa", "guide_counter_one_mismatch")
    sanson_dm_100k = row_by(latest_100k, "sanson_brunello", "dotmatch_hamming_k1")
    sanson_gc_100k = row_by(latest_100k, "sanson_brunello", "guide_counter_one_mismatch")
    sanson_dm_full = row_by(rows, "sanson_brunello", "dotmatch_hamming_k1", "full")
    sanson_gc_full = row_by(rows, "sanson_brunello", "guide_counter_one_mismatch", "full")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "\n".join([
            "# Atlas CRISPR Progress",
            "",
            "This report is intentionally not a SOTA claim. It summarizes the latest atlas real-data rows after the direct Hamming path, fused offset detection, block gzip reader, exact-only k=0 index, compact Hamming hash table, adaptive seed/precompute strategy, FASTQ length-return reader, and reserved sampling buffer work.",
            "",
            "## Current Read",
            "",
            f"- Yusa Hamming k=1 speedup vs guide-counter on the 100k/sample row: **{speedup(yusa_dm, yusa_gc)}**.",
            f"- Sanson/Brunello Hamming k=1 speedup vs guide-counter on the 100k/sample row: **{speedup(sanson_dm_100k, sanson_gc_100k)}**.",
            f"- Sanson/Brunello Hamming k=1 speedup vs guide-counter on the full 246.95M-read row: **{speedup(sanson_dm_full, sanson_gc_full)}**.",
            "- The full Sanson row is the serious stress test. It is useful evidence, but it is not 10x.",
            "- DotMatch uses less memory than guide-counter in these rows, but the current 2-vCPU atlas host limits how much memory-efficient parallelism can be demonstrated.",
            "",
            "## Figures",
            "",
            f"![Runtime]({runtime_fig.resolve()})",
            "",
            f"![Throughput]({throughput_fig.resolve()})",
            "",
            f"![Peak memory]({memory_fig.resolve()})",
            "",
            "## Rows",
            "",
            table(rows),
            "",
            "## Gate Status",
            "",
            "**CRISPR SOTA/10x gate: fail.** The data show real utility and some speed/memory wins, but not a defensible 10x claim against guide-counter on the full Sanson/Brunello workload.",
        ]) + "\n",
        encoding="utf-8",
    )
    print(OUT)


if __name__ == "__main__":
    main()
