#!/usr/bin/env python3
"""Generate CRISPR comparison evidence report from checked raw CSVs."""

from __future__ import annotations

import csv
import math
import os
import statistics
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
from check_crispr_comparison_gate import FULL_FASTQ_SAMPLE_READS


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw"
OUT_DIR = ROOT / "docs" / "benchmarks" / "crispr_comparison"
FIG_DIR = ROOT / "benchmarks" / "figures"


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh))


def fnum(value: str | None) -> float:
    if not value:
        return 0.0
    try:
        return float(value)
    except ValueError:
        return 0.0


def aggregate_full_sample_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    complete_keys: set[tuple[str, str, str]] = set()
    groups: dict[tuple[str, str, str], dict[str, dict[str, str]]] = {}
    for row in rows:
        if row.get("exit_code") != "0":
            continue
        if row.get("requested_records_per_sample") != "full" or row.get("run_level") != "full_sample":
            continue
        dataset = row.get("dataset_id", "")
        expected = FULL_FASTQ_SAMPLE_READS.get(dataset, {})
        sample_id = row.get("sample_id", "")
        if sample_id not in expected:
            continue
        if fnum(row.get("n_reads")) < expected[sample_id]:
            continue
        key = (dataset, row.get("tool", ""), row.get("repeat", ""))
        current = groups.setdefault(key, {}).get(sample_id)
        if current is None or fnum(row.get("reads_per_sec")) > fnum(current.get("reads_per_sec")):
            groups[key][sample_id] = row

    aggregate_rows: list[dict[str, str]] = []
    for key, by_sample in groups.items():
        dataset, tool, repeat = key
        expected = FULL_FASTQ_SAMPLE_READS.get(dataset, {})
        if set(by_sample) != set(expected):
            continue
        complete_keys.add(key)
        sample_rows = list(by_sample.values())
        total_reads = sum(fnum(row.get("n_reads")) for row in sample_rows)
        total_seconds = sum(fnum(row.get("seconds")) for row in sample_rows)
        weighted_verified = [
            (fnum(row.get("verified_per_read")), fnum(row.get("n_reads")))
            for row in sample_rows
            if row.get("verified_per_read")
        ]
        row = dict(sample_rows[0])
        row["tool"] = tool
        row["dataset_id"] = dataset
        row["repeat"] = repeat
        row["run_level"] = "full_sample_aggregate"
        row["sample_id"] = ""
        row["n_reads"] = f"{total_reads:.0f}"
        row["seconds"] = f"{total_seconds:.6f}"
        row["reads_per_sec"] = f"{(total_reads / total_seconds):.1f}" if total_seconds > 0 else "0.0"
        row["peak_rss_kb"] = f"{max(fnum(r.get('peak_rss_kb')) for r in sample_rows):.0f}"
        if weighted_verified:
            numerator = sum(value * reads for value, reads in weighted_verified)
            denominator = sum(reads for _, reads in weighted_verified)
            row["verified_per_read"] = f"{(numerator / denominator):.4f}" if denominator else ""
        aggregate_rows.append(row)

    out: list[dict[str, str]] = []
    for row in rows:
        key = (row.get("dataset_id", ""), row.get("tool", ""), row.get("repeat", ""))
        if row.get("requested_records_per_sample") == "full" and row.get("run_level") == "full_sample" and key in complete_keys:
            continue
        out.append(row)
    out.extend(aggregate_rows)
    return out


def p95(values: list[float]) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    return ordered[max(0, min(len(ordered) - 1, math.ceil(0.95 * len(ordered)) - 1))]


def repeated_stats(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for row in aggregate_full_sample_rows(rows):
        if row.get("exit_code") != "0":
            continue
        key = (row.get("dataset_id", ""), row.get("tool", ""), row.get("requested_records_per_sample", ""))
        groups.setdefault(key, []).append(row)
    out: list[dict[str, str]] = []
    for (dataset, tool, requested), group in sorted(groups.items()):
        reads_s = [fnum(r.get("reads_per_sec")) for r in group]
        seconds = [fnum(r.get("seconds")) for r in group]
        rss = [fnum(r.get("peak_rss_kb")) / 1024.0 for r in group]
        verified = [fnum(r.get("verified_per_read")) for r in group if r.get("verified_per_read")]
        mean = statistics.mean(reads_s) if reads_s else 0.0
        stdev = statistics.stdev(reads_s) if len(reads_s) > 1 else 0.0
        out.append({
            "dataset": dataset,
            "tool": tool,
            "records_per_sample": requested,
            "repeats": str(len(group)),
            "mean_reads_per_sec": f"{mean:.1f}",
            "p50_reads_per_sec": f"{statistics.median(reads_s):.1f}" if reads_s else "0.0",
            "p95_reads_per_sec": f"{p95(reads_s):.1f}",
            "mean_seconds": f"{statistics.mean(seconds):.4f}" if seconds else "0.0000",
            "cv": f"{(stdev / mean):.4f}" if mean else "0.0000",
            "max_peak_rss_mb": f"{max(rss):.1f}" if rss else "",
            "mean_verified_per_read": f"{statistics.mean(verified):.3f}" if verified else "",
        })
    return out


def markdown_table(rows: list[dict[str, str]], cols: list[str]) -> str:
    if not rows:
        return "_No rows available._\n"
    lines = ["|" + "|".join(cols) + "|", "|" + "|".join(["---"] * len(cols)) + "|"]
    for row in rows:
        lines.append("|" + "|".join(str(row.get(c, "")) for c in cols) + "|")
    return "\n".join(lines) + "\n"


def full_hamming_speed_rows(stats: list[dict[str, str]]) -> list[dict[str, str]]:
    by_key = {(r.get("dataset", ""), r.get("tool", ""), r.get("records_per_sample", "")): r for r in stats}
    datasets = sorted({r.get("dataset", "") for r in stats if r.get("records_per_sample") == "full"})
    out: list[dict[str, str]] = []
    for dataset in datasets:
        dotmatch = by_key.get((dataset, "dotmatch_hamming_k1", "full"))
        guide_counter = by_key.get((dataset, "guide_counter_one_mismatch", "full"))
        dm_rps = fnum(dotmatch.get("mean_reads_per_sec") if dotmatch else "")
        gc_rps = fnum(guide_counter.get("mean_reads_per_sec") if guide_counter else "")
        speedup = dm_rps / gc_rps if dm_rps > 0.0 and gc_rps > 0.0 else 0.0
        status = "pass" if speedup >= 1.0 else "blocked"
        if not dotmatch or not guide_counter:
            status = "missing"
        out.append({
            "dataset": dataset,
            "dotmatch_hamming_reads_per_sec": f"{dm_rps:.1f}" if dotmatch else "",
            "guide_counter_reads_per_sec": f"{gc_rps:.1f}" if guide_counter else "",
            "speedup": f"{speedup:.2f}" if dotmatch and guide_counter else "",
            "status": status,
        })
    return out


def svg_bars(stats: list[dict[str, str]], path: Path) -> None:
    selected = [r for r in stats if r["records_per_sample"] != "full"]
    if not selected:
        return
    labels = [f"{r['dataset']} {r['tool']} {r['records_per_sample']}" for r in selected]
    values = [fnum(r["mean_reads_per_sec"]) for r in selected]
    width = 1220
    row_h = 28
    left = 560
    height = 70 + row_h * len(selected)
    max_v = max(values) or 1.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,sans-serif;font-size:12px}.title{font-size:18px;font-weight:700}.axis{fill:#444}.bar{fill:#2f7d68}</style>',
        '<text class="title" x="20" y="28">CRISPR comparison repeated real-data throughput</text>',
        '<text class="axis" x="20" y="50">Mean reads/s; rows are separated by dataset, tool, and records/sample</text>',
    ]
    for i, (label, value) in enumerate(zip(labels, values)):
        y = 75 + i * row_h
        w = max(1, int((width - left - 120) * value / max_v))
        parts.append(f'<text x="20" y="{y + 14}">{label}</text>')
        parts.append(f'<rect class="bar" x="{left}" y="{y}" width="{w}" height="18" rx="2"/>')
        parts.append(f'<text x="{left + w + 8}" y="{y + 14}">{value:.1f}</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def markdown_link_path(target: Path, base: Path) -> str:
    return Path(os.path.relpath(target, start=base)).as_posix()


def main() -> None:
    repeated = read_rows(RAW / "crispr_comparison_repeated.csv")
    validation = read_rows(RAW / "crispr_comparison_edlib_validation.csv")
    agreement = read_rows(RAW / "crispr_comparison_count_agreement_summary.csv")
    stats = repeated_stats(repeated)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    svg_bars(stats, FIG_DIR / "crispr_comparison_throughput.svg")

    full_rows = [r for r in stats if r["records_per_sample"] == "full"]
    subsample_rows = [r for r in stats if r["records_per_sample"] != "full"]
    content = [
        "# CRISPR Comparison Evidence",
        "",
        "This report is generated from raw CSV files. It is intentionally stricter than the public quick-check report: comparison rows require both MAGeCK/Yusa and Sanson/Brunello real-data rows, competitor rows, count agreement, and Edlib validation.",
        "",
        "## Evidence Boundary",
        "",
        "- Hamming `k=1` rows are the fair guide-counter lane: one mismatch, no indels.",
        "- Levenshtein `k=1` rows are the DotMatch differentiator lane: substitutions plus single-base insertions/deletions, with Edlib validation.",
        "- Full FASTQ rows are reported separately from repeated subsamples.",
        "- Broad comparisons require the `make crispr-comparison-gate` check to pass.",
        "",
        "## Throughput Figure",
        "",
        "![CRISPR comparison throughput](" +
        markdown_link_path(FIG_DIR / "crispr_comparison_throughput.svg", OUT_DIR) + ")",
        "",
        "## Repeated Subsample Rows",
        "",
        markdown_table(subsample_rows, [
            "dataset", "tool", "records_per_sample", "repeats", "mean_reads_per_sec",
            "p50_reads_per_sec", "p95_reads_per_sec", "cv", "max_peak_rss_mb",
            "mean_verified_per_read",
        ]),
        "",
        "## Full FASTQ Rows",
        "",
        markdown_table(full_rows, [
            "dataset", "tool", "records_per_sample", "repeats", "mean_reads_per_sec",
            "mean_seconds", "max_peak_rss_mb", "mean_verified_per_read",
        ]),
        "",
        "## Full Hamming Speed Check",
        "",
        markdown_table(full_hamming_speed_rows(stats), [
            "dataset", "dotmatch_hamming_reads_per_sec", "guide_counter_reads_per_sec",
            "speedup", "status",
        ]),
        "",
        "## Edlib Validation",
        "",
        markdown_table(validation, [
            "dataset", "sample", "checked_reads", "mismatches", "oracle_strategy",
            "edlib_alignments", "bounded_windows", "fallback_windows", "selected_target_start",
            "stratum_exact", "stratum_corrected", "stratum_ambiguous", "stratum_unmatched",
            "stratum_contains_n",
        ]),
        "",
        "## Count Agreement",
        "",
        markdown_table(agreement, [
            "dataset", "comparison", "status", "n_guides", "total_delta",
            "differing_guides", "max_abs_delta", "pearson", "spearman",
        ]),
        "",
        "## Raw Inputs",
        "",
        "- `benchmarks/raw/crispr_comparison_repeated.csv`",
        "- `benchmarks/raw/crispr_comparison_edlib_validation.csv`",
        "- `benchmarks/raw/crispr_comparison_count_agreement_summary.csv`",
    ]
    (OUT_DIR / "README.md").write_text("\n".join(content) + "\n", encoding="utf-8")
    print(OUT_DIR / "README.md")


if __name__ == "__main__":
    main()
