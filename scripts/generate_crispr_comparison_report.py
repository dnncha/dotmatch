#!/usr/bin/env python3
"""Generate CRISPR comparison evidence report from gate-grade raw CSVs."""

from __future__ import annotations

import csv
import html
import json
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
HAMMING_K23_COMPARATOR_CSV = RAW / "crispr_comparison_hamming_k23_comparators.csv"
OUT_DIR = ROOT / "docs" / "benchmarks" / "crispr_comparison"
FIG_DIR = ROOT / "benchmarks" / "figures"
OPTIMIZER_ARTIFACTS = {
    "sanson_brunello": RAW / "crispr_sanson_brunello_backend_optimization_atlas_latest_dotmatch.json",
}


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


def first_value(row: dict[str, str], keys: list[str]) -> str:
    for key in keys:
        value = row.get(key, "")
        if value not in (None, ""):
            return str(value)
    return ""


def row_k(row: dict[str, str]) -> str:
    for key in ("k", "hamming_k", "max_mismatches", "mismatches"):
        value = first_value(row, [key]).strip()
        if value:
            return value.removeprefix("k")
    for key in ("comparison", "dotmatch_tool", "tool", "semantics"):
        text = str(row.get(key, "")).lower()
        for k in ("2", "3"):
            if f"k{k}" in text or f"k={k}" in text or f"hamming_{k}" in text:
                return k
    return ""


def has_dotmatch_hamming(row: dict[str, str]) -> bool:
    text = " ".join(str(row.get(key, "")) for key in ("comparison", "dotmatch_tool", "tool", "left_tool", "right_tool")).lower()
    return "dotmatch" in text and "hamming" in text


def has_bowtie1(row: dict[str, str]) -> bool:
    text = " ".join(str(row.get(key, "")) for key in ("comparison", "bowtie1_tool", "bowtie_tool", "tool", "left_tool", "right_tool", "comparator")).lower()
    return "bowtie1" in text or "bowtie_1" in text or "bowtie 1" in text


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


def optimizer_rows() -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for dataset, path in sorted(OPTIMIZER_ARTIFACTS.items()):
        if not path.exists():
            continue
        with path.open() as fh:
            data = json.load(fh)
        out.append({
            "dataset": dataset,
            "optimizer": str(data.get("optimizer", "")),
            "authority": str(data.get("authority", "")),
            "selected_backend": str(data.get("selected_backend", "")),
            "candidate_backend": str(data.get("candidate_backend", "")),
            "recommendation": str(data.get("recommendation", "")),
            "expected_speedup_band": str(data.get("expected_speedup_band", "")),
            "estimated_total_speedup": str(data.get("estimated_total_speedup", "")),
        })
    return out


def full_hamming_ratio_rows(stats: list[dict[str, str]]) -> list[dict[str, str]]:
    by_key = {(r.get("dataset", ""), r.get("tool", ""), r.get("records_per_sample", "")): r for r in stats}
    datasets = sorted({r.get("dataset", "") for r in stats if r.get("records_per_sample") == "full"})
    out: list[dict[str, str]] = []
    for dataset in datasets:
        dotmatch = by_key.get((dataset, "dotmatch_hamming_k1", "full"))
        guide_counter = by_key.get((dataset, "guide_counter_one_mismatch", "full"))
        dm_rps = fnum(dotmatch.get("mean_reads_per_sec") if dotmatch else "")
        gc_rps = fnum(guide_counter.get("mean_reads_per_sec") if guide_counter else "")
        speedup = dm_rps / gc_rps if dm_rps > 0.0 and gc_rps > 0.0 else 0.0
        status = "reported"
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


def hamming_k23_comparator_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for row in rows:
        k = row_k(row)
        if k not in {"2", "3"}:
            continue
        if not (has_dotmatch_hamming(row) and has_bowtie1(row)):
            continue
        dotmatch_rps = fnum(first_value(row, ["dotmatch_reads_per_sec", "dotmatch_mean_reads_per_sec", "left_reads_per_sec"]))
        bowtie_rps = fnum(first_value(row, ["bowtie1_reads_per_sec", "bowtie_reads_per_sec", "bowtie1_mean_reads_per_sec", "right_reads_per_sec"]))
        speedup_value = first_value(row, ["speedup", "dotmatch_vs_bowtie1_speedup"])
        speedup = fnum(speedup_value) if speedup_value else (dotmatch_rps / bowtie_rps if dotmatch_rps > 0 and bowtie_rps > 0 else 0.0)
        out.append({
            "dataset": first_value(row, ["dataset", "dataset_id", "workflow"]),
            "k": k,
            "records_per_sample": first_value(row, ["records_per_sample", "requested_records_per_sample", "n_reads"]),
            "dotmatch_tool": first_value(row, ["dotmatch_tool", "left_tool"]) or f"dotmatch_hamming_k{k}",
            "bowtie1_tool": first_value(row, ["bowtie1_tool", "bowtie_tool", "right_tool"]) or "bowtie1",
            "dotmatch_reads_per_sec": f"{dotmatch_rps:.1f}" if dotmatch_rps else "",
            "bowtie1_reads_per_sec": f"{bowtie_rps:.1f}" if bowtie_rps else "",
            "speedup": f"{speedup:.2f}" if speedup else "",
            "status": first_value(row, ["status"]) or ("reported" if dotmatch_rps and bowtie_rps else ""),
            "semantics": first_value(row, ["semantics"]) or f"Hamming k={k}, no indels",
            "artifact": first_value(row, ["artifact", "source_artifact"]),
        })
    return sorted(out, key=lambda r: (r["dataset"], r["k"], r["records_per_sample"]))


def guide_counter_style_rows(stats: list[dict[str, str]], agreement: list[dict[str, str]]) -> list[dict[str, str]]:
    by_key = {(r.get("dataset", ""), r.get("tool", ""), r.get("records_per_sample", "")): r for r in stats}
    agreement_by_dataset: dict[str, dict[str, str]] = {}
    for row in agreement:
        if row.get("comparison", "").endswith("dotmatch_hamming_vs_guide_counter"):
            agreement_by_dataset[row.get("dataset", "")] = row

    keys = sorted(
        (dataset, records)
        for dataset, tool, records in by_key
        if tool == "dotmatch_hamming_k1" and (dataset, "guide_counter_one_mismatch", records) in by_key
    )
    out: list[dict[str, str]] = []
    for dataset, records in keys:
        dotmatch = by_key[(dataset, "dotmatch_hamming_k1", records)]
        guide_counter = by_key[(dataset, "guide_counter_one_mismatch", records)]
        dm_rps = fnum(dotmatch.get("mean_reads_per_sec"))
        gc_rps = fnum(guide_counter.get("mean_reads_per_sec"))
        speedup = dm_rps / gc_rps if dm_rps > 0.0 and gc_rps > 0.0 else 0.0
        agreement_row = agreement_by_dataset.get(dataset, {})
        out.append({
            "dataset": dataset,
            "records_per_sample": records,
            "dotmatch_hamming_reads_per_sec": f"{dm_rps:.1f}",
            "guide_counter_reads_per_sec": f"{gc_rps:.1f}",
            "speedup": f"{speedup:.2f}" if speedup else "",
            "count_agreement_status": agreement_row.get("status", ""),
            "count_total_delta": agreement_row.get("total_delta", ""),
            "semantics": "one mismatch, no indels",
        })
    return out


def svg_bars(stats: list[dict[str, str]], path: Path) -> None:
    selected = [r for r in stats if fnum(r.get("mean_reads_per_sec")) > 0.0]
    if not selected:
        return
    labels = [
        f"{r['dataset']} {r['tool']} {r['records_per_sample']}"
        + (" FASTQs" if r["records_per_sample"] == "full" else "")
        for r in selected
    ]
    values = [fnum(r["mean_reads_per_sec"]) for r in selected]
    width = 1220
    row_h = 28
    left = 560
    height = 70 + row_h * len(selected)
    max_v = max(values) or 1.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,sans-serif;font-size:12px}.title{font-size:18px;font-weight:700}.axis{fill:#444}.bar{fill:#2f7d68}.bar-full{fill:#8b5e34}</style>',
        '<text class="title" x="20" y="28">CRISPR guide-counting throughput comparison</text>',
        '<text class="axis" x="20" y="50">Mean reads/s; includes repeated subsamples and available full FASTQ paper-data rows</text>',
    ]
    for i, (row, label, value) in enumerate(zip(selected, labels, values)):
        y = 75 + i * row_h
        w = max(1, int((width - left - 120) * value / max_v))
        klass = "bar-full" if row["records_per_sample"] == "full" else "bar"
        parts.append(f'<text x="20" y="{y + 14}">{html.escape(label)}</text>')
        parts.append(f'<rect class="{klass}" x="{left}" y="{y}" width="{w}" height="18" rx="2"/>')
        parts.append(f'<text x="{left + w + 8}" y="{y + 14}">{value:.1f}</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def svg_hamming_k23_comparators(rows: list[dict[str, str]], path: Path) -> None:
    selected = [r for r in rows if fnum(r.get("dotmatch_reads_per_sec")) > 0.0 and fnum(r.get("bowtie1_reads_per_sec")) > 0.0]
    if not selected:
        return
    values: list[tuple[str, str, float]] = []
    for row in selected:
        label = f"{row['dataset']} k{row['k']} {row['records_per_sample']}"
        values.append((label, "DotMatch", fnum(row["dotmatch_reads_per_sec"])))
        values.append((label, "Bowtie 1", fnum(row["bowtie1_reads_per_sec"])))
    width = 1220
    row_h = 26
    left = 420
    height = 74 + row_h * len(values)
    max_v = max(value for _, _, value in values) or 1.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,sans-serif;font-size:12px}.title{font-size:18px;font-weight:700}.axis{fill:#444}.dotmatch{fill:#2f7d68}.bowtie{fill:#5b6f95}</style>',
        '<text class="title" x="20" y="28">Hamming k2/k3 fixed-window comparator throughput</text>',
        '<text class="axis" x="20" y="50">Reads/s on extracted guide-window semantics; Bowtie 1 uses -v K --best --strata --norc -a</text>',
    ]
    for i, (label, tool, value) in enumerate(values):
        y = 75 + i * row_h
        w = max(1, int((width - left - 120) * value / max_v))
        klass = "dotmatch" if tool == "DotMatch" else "bowtie"
        parts.append(f'<text x="20" y="{y + 14}">{html.escape(label)} {tool}</text>')
        parts.append(f'<rect class="{klass}" x="{left}" y="{y}" width="{w}" height="17" rx="2"/>')
        parts.append(f'<text x="{left + w + 8}" y="{y + 13}">{value:.1f}</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def markdown_link_path(target: Path, base: Path) -> str:
    return Path(os.path.relpath(target, start=base)).as_posix()


def main() -> None:
    repeated = read_rows(RAW / "crispr_comparison_repeated.csv")
    validation = read_rows(RAW / "crispr_comparison_edlib_validation.csv")
    agreement = read_rows(RAW / "crispr_comparison_count_agreement_summary.csv")
    hamming_k23_comparators = read_rows(HAMMING_K23_COMPARATOR_CSV)
    stats = repeated_stats(repeated)
    hamming_k23_rows = hamming_k23_comparator_rows(hamming_k23_comparators)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    svg_bars(stats, FIG_DIR / "crispr_comparison_throughput.svg")
    svg_hamming_k23_comparators(hamming_k23_rows, FIG_DIR / "crispr_hamming_k23_comparison.svg")

    full_rows = [r for r in stats if r["records_per_sample"] == "full"]
    subsample_rows = [r for r in stats if r["records_per_sample"] != "full"]
    content = [
        "# CRISPR Comparison Evidence",
        "",
        "This report is generated from raw CSV artifacts. It is intentionally stricter than the public smoke report: comparison rows require both MAGeCK/Yusa and Sanson/Brunello real-data rows, competitor rows, count agreement, and Edlib validation.",
        "",
        "## Evidence Boundary",
        "",
        "- Hamming `k=1` rows are the fair guide-counter lane: one mismatch, no indels.",
        "- Hamming `k=2` and `k=3` external comparator rows are reported only from the separate DotMatch-vs-Bowtie 1 artifact when present.",
        "- Levenshtein `k=1` rows are the DotMatch differentiator lane: substitutions plus single-base insertions/deletions, with Edlib validation.",
        "- Full FASTQ rows are reported separately from repeated subsamples.",
        "- guide-counter speed ratios are reported when present; they are not universal replacement gates.",
        "- Broad comparisons require `make crispr-comparison-gate` to pass.",
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
        "## Guide-Counter-Style Public Paper-Data Lane",
        "",
        "DotMatch `dotmatch_hamming_k1` versus `guide_counter_one_mismatch` on the public paper-data inputs. This lane uses best-distance Hamming assignment with guide-counter's offset threshold, is limited to one mismatch and no indels, and keeps Levenshtein rows as a separate DotMatch capability lane.",
        "",
        markdown_table(guide_counter_style_rows(stats, agreement), [
            "dataset", "records_per_sample", "dotmatch_hamming_reads_per_sec",
            "guide_counter_reads_per_sec", "speedup", "count_agreement_status",
            "count_total_delta", "semantics",
        ]),
        "",
        "## Full Hamming k1 Guide-Counter Ratio",
        "",
        markdown_table(full_hamming_ratio_rows(stats), [
            "dataset", "dotmatch_hamming_reads_per_sec", "guide_counter_reads_per_sec",
            "speedup", "status",
        ]),
        "",
        "## Hamming k2/k3 External Comparator Rows",
        "",
        "DotMatch Hamming `k=2` and `k=3` evidence is kept separate from guide-counter claims. Rows in this table must come from `benchmarks/raw/crispr_comparison_hamming_k23_comparators.csv` and compare DotMatch directly with Bowtie 1.",
        "",
        "![Hamming k2/k3 comparator throughput](" +
        markdown_link_path(FIG_DIR / "crispr_hamming_k23_comparison.svg", OUT_DIR) + ")",
        "",
        markdown_table(hamming_k23_rows, [
            "dataset", "k", "records_per_sample", "dotmatch_tool", "bowtie1_tool",
            "dotmatch_reads_per_sec", "bowtie1_reads_per_sec", "speedup", "status",
            "semantics", "artifact",
        ]),
        "",
        "## Backend Optimizer",
        "",
        "The optimizer is advisory and CPU-authoritative: it records the fastest eligible candidate backend, but speed claims still require CPU checksum agreement.",
        "",
        markdown_table(optimizer_rows(), [
            "dataset", "optimizer", "authority", "selected_backend", "candidate_backend",
            "recommendation", "expected_speedup_band", "estimated_total_speedup",
        ]),
        "",
        "## Edlib Oracle Validation",
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
        "- `benchmarks/raw/crispr_comparison_full_sanson_atlas_latest_dotmatch.csv`",
        "- `benchmarks/raw/crispr_comparison_edlib_validation.csv`",
        "- `benchmarks/raw/crispr_comparison_count_agreement_summary.csv`",
        "- `benchmarks/raw/crispr_comparison_hamming_k23_comparators.csv`",
        "- `benchmarks/raw/crispr_sanson_brunello_backend_optimization_atlas_latest_dotmatch.json`",
    ]
    (OUT_DIR / "README.md").write_text("\n".join(content) + "\n", encoding="utf-8")
    print(OUT_DIR / "README.md")


if __name__ == "__main__":
    main()
