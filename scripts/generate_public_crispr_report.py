#!/usr/bin/env python3
"""Generate the public CRISPR workflow comparator report."""

from __future__ import annotations

import csv
import math
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "public_crispr_workflow.csv"
RAW_REPEATED = ROOT / "benchmarks" / "raw" / "public_crispr_repeated.csv"
COUNT_AGREEMENT = ROOT / "benchmarks" / "raw" / "count_agreement_summary.csv"
EDLIB_VALIDATION = ROOT / "benchmarks" / "raw" / "public_crispr_edlib_validation.csv"
SAMPLE_SCALING = ROOT / "benchmarks" / "raw" / "public_crispr_sample_scaling.csv"
OUT_DIR = ROOT / "docs" / "benchmarks" / "public_crispr"
FIG_DIR = ROOT / "benchmarks" / "figures"


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh))


def svg_bar(rows: list[dict[str, str]], path: Path) -> None:
    selected = [r for r in rows if r.get("exit_code") == "0"]
    if not selected:
        return
    labels = [r["tool"] for r in selected]
    values = [float(r.get("reads_per_sec") or 0.0) for r in selected]
    width = 900
    row_h = 34
    left = 260
    height = 70 + row_h * len(selected)
    max_v = max(values) if values else 1.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,sans-serif;font-size:13px}.title{font-size:18px;font-weight:700}.axis{fill:#444}.bar{fill:#2f7d68}</style>',
        '<text class="title" x="20" y="28">Public CRISPR workflow throughput</text>',
        '<text class="axis" x="20" y="50">Small real MAGeCK/Yusa FASTQ subset; higher is better</text>',
    ]
    for i, (label, value) in enumerate(zip(labels, values)):
        y = 75 + i * row_h
        bar_w = 1 if max_v == 0 else int((width - left - 80) * value / max_v)
        parts.append(f'<text x="20" y="{y + 15}">{label}</text>')
        parts.append(f'<rect class="bar" x="{left}" y="{y}" width="{bar_w}" height="20" rx="2"/>')
        parts.append(f'<text x="{left + bar_w + 8}" y="{y + 15}">{value:.1f} reads/s</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def p95(values: list[float]) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    idx = max(0, min(len(ordered) - 1, math.ceil(0.95 * len(ordered)) - 1))
    return ordered[idx]


def repeated_stats(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for row in rows:
        if row.get("exit_code") != "0":
            continue
        key = (row.get("tool", ""), row.get("semantics", ""), row.get("requested_records_per_sample", row.get("n_reads", "")))
        groups.setdefault(key, []).append(row)
    out: list[dict[str, str]] = []
    for (tool, semantics, requested), group in sorted(groups.items()):
        values = [float(r.get("reads_per_sec") or 0.0) for r in group]
        seconds_values = [float(r.get("seconds") or 0.0) for r in group]
        rss_values = [float(r.get("peak_rss_kb") or 0.0) for r in group]
        verified = [float(r.get("verified_per_read") or 0.0) for r in group if r.get("verified_per_read")]
        mean = statistics.mean(values)
        stdev = statistics.stdev(values) if len(values) > 1 else 0.0
        out.append({
            "tool": tool,
            "semantics": semantics,
            "records_per_sample": requested,
            "repeats": str(len(values)),
            "mean_reads_per_sec": f"{mean:.1f}",
            "p50_reads_per_sec": f"{statistics.median(values):.1f}",
            "p95_reads_per_sec": f"{p95(values):.1f}",
            "mean_seconds": f"{statistics.mean(seconds_values):.4f}",
            "p50_seconds": f"{statistics.median(seconds_values):.4f}",
            "cv": f"{(stdev / mean):.4f}" if mean else "0.0000",
            "max_peak_rss_mb": f"{(max(rss_values) / 1024.0):.1f}" if rss_values else "",
            "mean_verified_per_read": f"{statistics.mean(verified):.3f}" if verified else "",
        })
    return out


def svg_repeated(rows: list[dict[str, str]], path: Path) -> None:
    stats = repeated_stats(rows)
    selected = [r for r in stats if r["tool"].startswith("dotmatch") or r["tool"].startswith("guide_counter") or r["tool"].startswith("mageck")]
    if not selected:
        return
    labels = [f"{r['tool']} ({r['records_per_sample']}/sample)" for r in selected]
    values = [float(r["mean_reads_per_sec"]) for r in selected]
    width = 1100
    row_h = 32
    left = 430
    height = 70 + row_h * len(selected)
    max_v = max(values) if values else 1.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,sans-serif;font-size:13px}.title{font-size:18px;font-weight:700}.axis{fill:#444}.bar{fill:#476fb3}</style>',
        '<text class="title" x="20" y="28">Repeated public CRISPR throughput</text>',
        '<text class="axis" x="20" y="50">Mean reads/s across repeated MAGeCK/Yusa real FASTQ subsamples</text>',
    ]
    for i, (label, value) in enumerate(zip(labels, values)):
        y = 75 + i * row_h
        bar_w = 1 if max_v == 0 else int((width - left - 90) * value / max_v)
        parts.append(f'<text x="20" y="{y + 15}">{label}</text>')
        parts.append(f'<rect class="bar" x="{left}" y="{y}" width="{bar_w}" height="19" rx="2"/>')
        parts.append(f'<text x="{left + bar_w + 8}" y="{y + 14}">{value:.1f}</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def svg_count_agreement(rows: list[dict[str, str]], path: Path) -> None:
    ok = [r for r in rows if r.get("status") == "ok"]
    if not ok:
        return
    width = 900
    row_h = 34
    left = 360
    height = 70 + row_h * len(ok)
    values = [float(r.get("pearson") or 0.0) for r in ok]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,sans-serif;font-size:13px}.title{font-size:18px;font-weight:700}.axis{fill:#444}.bar{fill:#7f5aa2}</style>',
        '<text class="title" x="20" y="28">Count-table agreement</text>',
        '<text class="axis" x="20" y="50">Pearson correlation by guide count table; 1.0 is exact linear agreement</text>',
    ]
    for i, (row, value) in enumerate(zip(ok, values)):
        y = 75 + i * row_h
        bar_w = int((width - left - 80) * max(0.0, min(1.0, value)))
        parts.append(f'<text x="20" y="{y + 15}">{row.get("comparison", "")}</text>')
        parts.append(f'<rect class="bar" x="{left}" y="{y}" width="{bar_w}" height="20" rx="2"/>')
        parts.append(f'<text x="{left + bar_w + 8}" y="{y + 15}">{value:.6f}</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def svg_metric_bar(rows: list[dict[str, str]], path: Path, title: str, subtitle: str,
                   value_key: str, label_suffix: str = "", lower_is_better: bool = False) -> None:
    selected = [r for r in rows if r.get("exit_code") == "0"]
    values = [float(r.get(value_key) or 0.0) for r in selected]
    if not selected or not values:
        return
    labels = [r.get("tool", "") for r in selected]
    width = 1050
    row_h = 34
    left = 310
    height = 74 + row_h * len(selected)
    max_v = max(values) if values else 1.0
    color = "#a45f3f" if lower_is_better else "#2f7d68"
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<style>text{{font-family:Arial,sans-serif;font-size:13px}}.title{{font-size:18px;font-weight:700}}.axis{{fill:#444}}.bar{{fill:{color}}}</style>',
        f'<text class="title" x="20" y="28">{title}</text>',
        f'<text class="axis" x="20" y="50">{subtitle}</text>',
    ]
    for i, (label, value) in enumerate(zip(labels, values)):
        y = 78 + i * row_h
        bar_w = 1 if max_v == 0 else int((width - left - 95) * value / max_v)
        parts.append(f'<text x="20" y="{y + 15}">{label}</text>')
        parts.append(f'<rect class="bar" x="{left}" y="{y}" width="{bar_w}" height="20" rx="2"/>')
        parts.append(f'<text x="{left + bar_w + 8}" y="{y + 15}">{value:.2f}{label_suffix}</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def svg_repeated_metric(stats: list[dict[str, str]], path: Path, title: str, subtitle: str,
                        value_key: str, label_suffix: str = "", color: str = "#476fb3") -> None:
    selected = [r for r in stats if r["tool"].startswith("dotmatch") or r["tool"].startswith("guide_counter") or r["tool"].startswith("mageck")]
    values = [float(r.get(value_key) or 0.0) for r in selected]
    if not selected or not values:
        return
    labels = [f"{r['tool']} ({r['records_per_sample']}/sample)" for r in selected]
    width = 1120
    row_h = 32
    left = 430
    height = 72 + row_h * len(selected)
    max_v = max(values) if values else 1.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<style>text{{font-family:Arial,sans-serif;font-size:13px}}.title{{font-size:18px;font-weight:700}}.axis{{fill:#444}}.bar{{fill:{color}}}</style>',
        f'<text class="title" x="20" y="28">{title}</text>',
        f'<text class="axis" x="20" y="50">{subtitle}</text>',
    ]
    for i, (label, value) in enumerate(zip(labels, values)):
        y = 76 + i * row_h
        bar_w = 1 if max_v == 0 else int((width - left - 90) * value / max_v)
        parts.append(f'<text x="20" y="{y + 15}">{label}</text>')
        parts.append(f'<rect class="bar" x="{left}" y="{y}" width="{bar_w}" height="19" rx="2"/>')
        parts.append(f'<text x="{left + bar_w + 8}" y="{y + 14}">{value:.2f}{label_suffix}</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def svg_assignment_impact(rows: list[dict[str, str]], path: Path) -> None:
    selected = [r for r in rows if r.get("exit_code") == "0" and (r.get("assigned_reads") or r.get("rejected_reads"))]
    if not selected:
        return
    width = 1120
    row_h = 38
    left = 310
    height = 78 + row_h * len(selected)
    max_total = max(float(r.get("n_reads") or 0.0) for r in selected) or 1.0
    colors = {
        "exact": "#2f7d68",
        "corrected": "#476fb3",
        "ambiguous": "#d19a3c",
        "rejected": "#8f4f4f",
    }
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,sans-serif;font-size:13px}.title{font-size:18px;font-weight:700}.axis{fill:#444}</style>',
        '<text class="title" x="20" y="28">Assignment impact on real CRISPR reads</text>',
        '<text class="axis" x="20" y="50">Exact assigned, one-edit rescued, ambiguous, and rejected reads</text>',
    ]
    for i, row in enumerate(selected):
        y = 82 + i * row_h
        exact = float(row.get("exact_reads") or 0.0)
        corrected = float(row.get("corrected_reads") or 0.0)
        ambiguous = float(row.get("ambiguous_reads") or 0.0)
        rejected = float(row.get("rejected_reads") or 0.0)
        values = [("exact", exact), ("corrected", corrected), ("ambiguous", ambiguous), ("rejected", rejected)]
        x = left
        parts.append(f'<text x="20" y="{y + 15}">{row.get("tool", "")}</text>')
        for key, value in values:
            w = int((width - left - 245) * value / max_total)
            if w > 0:
                parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="21" rx="2" fill="{colors[key]}"/>')
            x += w
        parts.append(f'<text x="{left + width - left - 230}" y="{y + 15}">exact {exact:.0f}, rescued {corrected:.0f}, ambig {ambiguous:.0f}, rejected {rejected:.0f}</text>')
    legend_y = 62
    legend_x = left
    for key in ["exact", "corrected", "ambiguous", "rejected"]:
        parts.append(f'<rect x="{legend_x}" y="{legend_y - 11}" width="12" height="12" fill="{colors[key]}"/>')
        parts.append(f'<text x="{legend_x + 16}" y="{legend_y}">{key}</text>')
        legend_x += 105
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def speedup_rows(stats: list[dict[str, str]]) -> list[dict[str, str]]:
    by_size = {r["records_per_sample"]: r for r in stats if r["tool"] == "dotmatch_hamming_k1"}
    out: list[dict[str, str]] = []
    for row in stats:
        size = row["records_per_sample"]
        if size not in by_size or row["tool"].startswith("dotmatch"):
            continue
        dot = float(by_size[size].get("mean_reads_per_sec") or 0.0)
        other = float(row.get("mean_reads_per_sec") or 0.0)
        if dot <= 0 or other <= 0:
            continue
        out.append({
            "baseline": row["tool"],
            "records_per_sample": size,
            "dotmatch_hamming_reads_per_sec": f"{dot:.1f}",
            "baseline_reads_per_sec": f"{other:.1f}",
            "speedup": f"{dot / other:.2f}x",
        })
    return out


def svg_scaling(rows: list[dict[str, str]], path: Path, metric: str, title: str, subtitle: str,
                suffix: str = "", color_a: str = "#2f7d68", color_b: str = "#8f4f4f") -> None:
    selected = [r for r in rows if r.get("exit_code") == "0"]
    if not selected:
        return
    sample_counts = sorted({int(r["n_samples"]) for r in selected})
    tools = sorted({r["tool"] for r in selected})
    width = 980
    height = 90 + 58 * len(sample_counts)
    left = 120
    group_w = width - left - 80
    max_v = max(float(r.get(metric) or 0.0) for r in selected) or 1.0
    colors = {tools[0]: color_a}
    if len(tools) > 1:
        colors[tools[1]] = color_b
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,sans-serif;font-size:13px}.title{font-size:18px;font-weight:700}.axis{fill:#444}</style>',
        f'<text class="title" x="20" y="28">{title}</text>',
        f'<text class="axis" x="20" y="50">{subtitle}</text>',
    ]
    for si, n_samples in enumerate(sample_counts):
        y = 82 + si * 58
        parts.append(f'<text x="20" y="{y + 23}">{n_samples} samples</text>')
        for ti, tool in enumerate(tools):
            row = next((r for r in selected if int(r["n_samples"]) == n_samples and r["tool"] == tool), None)
            if row is None:
                continue
            value = float(row.get(metric) or 0.0)
            bar_w = int(group_w * value / max_v)
            yy = y + ti * 24
            parts.append(f'<rect x="{left}" y="{yy}" width="{bar_w}" height="18" rx="2" fill="{colors.get(tool, "#476fb3")}"/>')
            parts.append(f'<text x="{left + bar_w + 8}" y="{yy + 14}">{tool}: {value:.1f}{suffix}</text>')
        if si == 0:
            parts.append(f'<text x="{left}" y="70">higher is better unless noted</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def md_table(rows: list[dict[str, str]], cols: list[str]) -> list[str]:
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for r in rows:
        lines.append("| " + " | ".join((r.get(c, "") or "").replace("|", "\\|") for c in cols) + " |")
    return lines


def main() -> None:
    rows = read_rows(RAW)
    repeated_rows = read_rows(RAW_REPEATED)
    count_rows = read_rows(COUNT_AGREEMENT)
    validation_rows = read_rows(EDLIB_VALIDATION)
    scaling_rows = read_rows(SAMPLE_SCALING)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = FIG_DIR / "public_crispr_throughput.svg"
    repeated_fig = FIG_DIR / "public_crispr_repeated_throughput.svg"
    agreement_fig = FIG_DIR / "public_crispr_count_agreement.svg"
    runtime_fig = FIG_DIR / "public_crispr_runtime_seconds.svg"
    memory_fig = FIG_DIR / "public_crispr_peak_memory.svg"
    impact_fig = FIG_DIR / "public_crispr_assignment_impact.svg"
    candidates_fig = FIG_DIR / "public_crispr_verified_candidates.svg"
    repeated_runtime_fig = FIG_DIR / "public_crispr_repeated_runtime_seconds.svg"
    repeated_memory_fig = FIG_DIR / "public_crispr_repeated_peak_memory.svg"
    repeated_candidates_fig = FIG_DIR / "public_crispr_repeated_verified_candidates.svg"
    scaling_throughput_fig = FIG_DIR / "public_crispr_sample_scaling_throughput.svg"
    scaling_memory_fig = FIG_DIR / "public_crispr_sample_scaling_memory.svg"
    svg_bar(rows, fig)
    svg_repeated(repeated_rows, repeated_fig)
    svg_count_agreement(count_rows, agreement_fig)
    svg_metric_bar(rows, runtime_fig, "Public CRISPR end-to-end runtime",
                   "Small real MAGeCK/Yusa FASTQ subset; lower is better", "seconds", " s", True)
    svg_metric_bar(rows, memory_fig, "Public CRISPR peak RSS",
                   "Measured per command with /usr/bin/time", "peak_rss_kb", " KB", True)
    svg_metric_bar(rows, candidates_fig, "Verified candidates per read",
                   "Mechanism figure: lookup-like candidate verification vs scan/workflow baselines",
                   "verified_per_read", " verified/read", True)
    svg_assignment_impact(rows, impact_fig)
    repeated_summary_for_figs = repeated_stats(repeated_rows)
    svg_repeated_metric(repeated_summary_for_figs, repeated_runtime_fig, "Repeated public CRISPR runtime",
                        "Mean end-to-end seconds across repeated real FASTQ subsamples; lower is better",
                        "mean_seconds", " s", "#a45f3f")
    svg_repeated_metric(repeated_summary_for_figs, repeated_memory_fig, "Repeated public CRISPR peak RSS",
                        "Maximum observed peak resident memory across repeated runs",
                        "max_peak_rss_mb", " MB", "#7f5aa2")
    svg_repeated_metric(repeated_summary_for_figs, repeated_candidates_fig, "Repeated verified candidates per read",
                        "DotMatch mechanism: exact verification work per input read",
                        "mean_verified_per_read", " verified/read", "#5d7b3f")
    svg_scaling(scaling_rows, scaling_throughput_fig, "reads_per_sec", "Multi-sample CRISPR scaling",
                "Real MAGeCK/Yusa FASTQs replicated as sample batches; DotMatch uses sample threads", " reads/s")
    svg_scaling(scaling_rows, scaling_memory_fig, "peak_rss_kb", "Multi-sample peak RSS",
                "Memory remains the lever that allows threaded sample batches", " KB", "#7f5aa2", "#a45f3f")

    lines: list[str] = [
        "# Public CRISPR Workflow Comparator",
        "",
        "This report tracks the MAGeCK/Yusa public CRISPR benchmark. The single-run table below is a smoke/latest wiring check only; repeated rows and SOTA-gated rows are the only rows intended to support publication-style performance claims.",
        "",
        "## Smoke/Latest Wiring Table",
        "",
        "**Not claim-grade.** These rows may be stale or reduced and must not be used for headline claims. Use the repeated-run statistics below, and use `docs/benchmarks/crispr_sota/README.md` once `make crispr-sota-gate` passes for two real CRISPR datasets.",
        "",
        *md_table(
            rows,
            [
                "tool",
                "version",
                "semantics",
                "n_reads",
                "n_targets",
                "seconds",
                "reads_per_sec",
                "peak_rss_kb",
                "assigned_reads",
                "corrected_reads",
                "ambiguous_reads",
                "rejected_reads",
                "overcount_reads",
                "verified_per_read",
                "exit_code",
            ],
        ),
        "",
        "![Public CRISPR throughput](../../../benchmarks/figures/public_crispr_throughput.svg)",
        "",
        "![Public CRISPR runtime](../../../benchmarks/figures/public_crispr_runtime_seconds.svg)",
        "",
        "![Public CRISPR memory](../../../benchmarks/figures/public_crispr_peak_memory.svg)",
        "",
        "![Public CRISPR assignment impact](../../../benchmarks/figures/public_crispr_assignment_impact.svg)",
        "",
        "![Public CRISPR verified candidates](../../../benchmarks/figures/public_crispr_verified_candidates.svg)",
        "",
        "## Repeated-Run Statistics",
        "",
    ]
    repeated_summary = repeated_summary_for_figs
    if repeated_summary:
        max_requested = 0
        for row in repeated_summary:
            try:
                max_requested = max(max_requested, int(row["records_per_sample"]))
            except ValueError:
                max_requested = 10000
        if max_requested < 10000:
            lines.extend([
                "**Current repeated CSV is a reduced harness check, not a claim-grade 10k/100k run.**",
                "",
            ])
        lines.extend(md_table(
            repeated_summary,
            [
                "tool",
                "semantics",
                "records_per_sample",
                "repeats",
                "mean_reads_per_sec",
                "p50_reads_per_sec",
                "p95_reads_per_sec",
                "mean_seconds",
                "p50_seconds",
                "cv",
                "max_peak_rss_mb",
                "mean_verified_per_read",
            ],
        ))
        lines.extend([
            "",
            "![Repeated public CRISPR throughput](../../../benchmarks/figures/public_crispr_repeated_throughput.svg)",
            "",
            "![Repeated public CRISPR runtime](../../../benchmarks/figures/public_crispr_repeated_runtime_seconds.svg)",
            "",
            "![Repeated public CRISPR peak memory](../../../benchmarks/figures/public_crispr_repeated_peak_memory.svg)",
            "",
            "![Repeated public CRISPR verified candidates](../../../benchmarks/figures/public_crispr_repeated_verified_candidates.svg)",
            "",
        ])
        speedups = speedup_rows(repeated_summary)
        if speedups:
            lines.extend([
                "## DotMatch Hamming Speedup",
                "",
                "This table keeps the fair CRISPR speed lane separate: DotMatch Hamming `k=1` versus tools with one-mismatch/no-indel or exact-count semantics.",
                "",
                *md_table(speedups, [
                    "baseline",
                    "records_per_sample",
                    "dotmatch_hamming_reads_per_sec",
                    "baseline_reads_per_sec",
                    "speedup",
                ]),
                "",
            ])
    else:
        lines.extend([
            "`benchmarks/raw/public_crispr_repeated.csv` has not been generated yet. Run `make bench-public-crispr-repeated`.",
            "",
        ])
    lines.extend([
        "## Count Agreement",
        "",
    ])
    if count_rows:
        lines.extend(md_table(
            count_rows,
            [
                "comparison",
                "status",
                "n_guides",
                "total_left",
                "total_right",
                "total_delta",
                "differing_guides",
                "max_abs_delta",
                "pearson",
                "spearman",
            ],
        ))
        lines.extend([
            "",
            "![Public CRISPR count agreement](../../../benchmarks/figures/public_crispr_count_agreement.svg)",
            "",
        ])
    else:
        lines.extend([
            "`benchmarks/raw/count_agreement_summary.csv` has not been generated yet. Run `make count-agreement` after a benchmark that includes guide-counter and MAGeCK.",
            "",
        ])
    lines.extend([
        "## Multi-Sample Scaling",
        "",
    ])
    if scaling_rows:
        lines.extend(md_table(
            scaling_rows,
            [
                "tool",
                "n_samples",
                "records_per_sample",
                "total_reads",
                "threads",
                "seconds",
                "reads_per_sec",
                "peak_rss_kb",
                "assigned_reads",
                "overcount_reads",
                "exit_code",
            ],
        ))
        lines.extend([
            "",
            "![Public CRISPR sample scaling throughput](../../../benchmarks/figures/public_crispr_sample_scaling_throughput.svg)",
            "",
            "![Public CRISPR sample scaling memory](../../../benchmarks/figures/public_crispr_sample_scaling_memory.svg)",
            "",
        ])
    else:
        lines.extend([
            "`benchmarks/raw/public_crispr_sample_scaling.csv` has not been generated yet. Run `make bench-public-crispr-scaling`.",
            "",
        ])
    lines.extend([
        "## Edlib Oracle Validation",
        "",
    ])
    if validation_rows:
        lines.extend(md_table(
            validation_rows,
            [
                "dataset",
                "sample",
                "oracle",
                "checked_reads",
                "mismatches",
                "indel_window",
                "stratum_exact",
                "stratum_corrected",
                "stratum_ambiguous",
                "stratum_unmatched",
                "stratum_contains_n",
            ],
        ))
        lines.append("")
    else:
        lines.extend([
            "`benchmarks/raw/public_crispr_edlib_validation.csv` has not been generated yet. Run `make validate-public-crispr-edlib`.",
            "",
        ])
    lines.extend([
        "## Interpretation",
        "",
        "- `dotmatch_hamming_k1` is the fair lane for guide-counter-style one-mismatch/no-indel guide counting.",
        "- `dotmatch_levenshtein_k1` is DotMatch's stronger lane: substitutions plus one-base insertions/deletions with explicit ambiguity reporting.",
        "- `dotmatch_exact_k0` is the fair exact-count lane for MAGeCK's direct FASTQ counting mode.",
        "- MAGeCK is run as exact FASTQ counting with `--trim-5 23`, matching the public Yusa demo workflow.",
        "- guide-counter is fast, but on the 10k Yusa run its own stats report more mapped reads than input reads, consistent with its multi-offset counting loop; DotMatch assigns at most one target per read and reports ambiguity instead.",
        "- The multi-sample scaling table is the current best memory-leverage evidence: DotMatch can process sample batches with threads while remaining tens of MB, whereas guide-counter uses roughly half a GB and its count total grows beyond input reads.",
        "- Cutadapt and Bowtie2 rows are workflow comparators on extracted guide windows; they are not exact assignment oracles.",
        "- Native Edlib scan remains the exact semantic oracle for assignment correctness.",
        "- README and manuscript speed claims should cite only repeated rows with zero validation mismatches and explicit semantics.",
        "",
        "## Raw Commands",
        "",
        *md_table(rows, ["tool", "command"]),
        "",
    ])
    (OUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print(OUT_DIR / "README.md")


if __name__ == "__main__":
    main()
