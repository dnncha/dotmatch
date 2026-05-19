#!/usr/bin/env python3
"""Generate BCL-demultiplexing benchmark report."""

from __future__ import annotations

import csv
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "bcl_demux.csv"
OUT_DIR = ROOT / "docs" / "benchmarks" / "bcl_demux"
FIG_DIR = ROOT / "benchmarks" / "figures"


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh))


def svg_bar(rows: list[dict[str, str]], key: str, title: str, subtitle: str, path: Path,
            suffix: str = "", color: str = "#2f7d68") -> None:
    selected = [r for r in rows if r.get("exit_code") == "0"]
    if not selected:
        return
    values = [float(r.get(key) or 0.0) for r in selected]
    labels = [r["tool"] for r in selected]
    width = 980
    left = 300
    row_h = 36
    height = 78 + row_h * len(selected)
    max_v = max(values) if values else 1.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<style>text{{font-family:Arial,sans-serif;font-size:13px}}.title{{font-size:18px;font-weight:700}}.axis{{fill:#444}}.bar{{fill:{color}}}</style>',
        f'<text class="title" x="20" y="28">{title}</text>',
        f'<text class="axis" x="20" y="50">{subtitle}</text>',
    ]
    for i, (label, value) in enumerate(zip(labels, values)):
        y = 82 + i * row_h
        w = 1 if max_v == 0 else int((width - left - 120) * value / max_v)
        parts.append(f'<text x="20" y="{y + 15}">{label}</text>')
        parts.append(f'<rect class="bar" x="{left}" y="{y}" width="{w}" height="21" rx="2"/>')
        parts.append(f'<text x="{left + w + 8}" y="{y + 15}">{value:.2f}{suffix}</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    rows = read_rows(RAW)
    if not rows:
        raise SystemExit(f"missing benchmark CSV: {RAW}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    throughput = FIG_DIR / "bcl_demux_throughput.svg"
    memory = FIG_DIR / "bcl_demux_peak_memory.svg"
    assigned = FIG_DIR / "bcl_demux_assigned_reads.svg"
    svg_bar(rows, "clusters_per_sec", "BCL demultiplexing throughput",
            "Classic BCL clusters processed per second; higher is better", throughput, " clusters/s")
    svg_bar(rows, "peak_rss_kb", "BCL demultiplexing peak memory",
            "Peak resident memory; lower is better", memory, " KB", "#a45f3f")
    svg_bar(rows, "assigned_reads", "BCL demultiplexing assigned reads",
            "Reads assigned to sample FASTQ outputs", assigned, " reads", "#476fb3")
    lines = [
        "# BCL Demultiplexing Benchmark",
        "",
        "This is the raw Illumina run-folder barcode demultiplexing evidence track.",
        "",
        "Current status: DotMatch supports a first lane-1 classic per-cycle BCL milestone. CBCL/NovaSeq-style input and broader multi-lane BCL operation are still future work. Broader comparisons need real run folders and zero-mismatch validation against BCL Convert or bcl2fastq where available.",
        "",
        "The public `public_10x_tiny_bcl` row uses the 10x Genomics Cell Ranger `tiny-bcl` mkfastq demo run folder. The bundled fetch script downloads the 10x sample sheet, the public run-folder archive, and the official Chromium i7 index CSV used to normalize the legacy `SI-P03-C9` index-set alias into concrete index sequences.",
        "",
        "Rows with exit code `127` are environment records for missing competitors, not runtime comparisons. They are kept in the raw table so the report is explicit about what was unavailable on this machine.",
        "",
    ]
    successful_dotmatch = [r for r in rows if r.get("tool") == "dotmatch_bcl_demux" and r.get("exit_code") == "0"]
    validated_competitors = [
        r for r in rows
        if r.get("tool") in {"bcl-convert", "bcl2fastq", "cuda-demux"}
        and r.get("exit_code") == "0"
        and r.get("validation_mismatches") == "0"
        and r.get("validation_exit_code") == "0"
    ]
    if successful_dotmatch and validated_competitors:
        best_dot = max(successful_dotmatch, key=lambda r: float(r.get("clusters_per_sec") or 0.0))
        best_comp = max(validated_competitors, key=lambda r: float(r.get("clusters_per_sec") or 0.0))
        dot_speed = float(best_dot.get("clusters_per_sec") or 0.0)
        comp_speed = float(best_comp.get("clusters_per_sec") or 0.0)
        speedup = dot_speed / comp_speed if comp_speed > 0 else 0.0
        lines.extend([
            "## Current Best Validated Comparison",
            "",
            f"Within this CSV, the fastest DotMatch row has {speedup:.2f}x higher throughput than the fastest validated installed comparator row ({best_comp.get('tool', '')}) on the same host and workflow.",
            "",
            "This is not a comparison result by itself: the broader comparison check still requires distinct repeated real runs, a successful DotMatch CBCL row, validated competitor rows, and stricter output validation where read names/paths are comparable.",
            "",
        ])
    lines.extend([
        "## Figures",
        "",
    ])
    for label, path in [("Throughput", throughput), ("Peak memory", memory), ("Assigned reads", assigned)]:
        if path.exists():
            lines.append(f"![{label}]({os.path.relpath(path, OUT_DIR)})")
            lines.append("")
    lines.extend([
        "## Raw Rows",
        "",
        "| tool | workflow | format | clusters | cycles | samples | threads | gzip | seconds | clusters/sec | peak RSS KB | output MB | assigned | undetermined | filtered | validation | mode | content hash | exit | logs |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: | --- |",
    ])
    for row in rows:
        row.setdefault("validation_mismatches", "")
        row.setdefault("requested_threads", "")
        output_mb = ""
        try:
            if row.get("output_bytes"):
                output_mb = f"{int(row['output_bytes']) / (1024 * 1024):.2f}"
        except ValueError:
            output_mb = ""
        content_hash = row.get("fastq_content_sha256", "")
        if content_hash:
            content_hash = content_hash[:12]
        logs = ""
        if row.get("stderr_log"):
            logs = os.path.relpath(row["stderr_log"], OUT_DIR)
        row.setdefault("validation_mode", "")
        row.setdefault("gzip_level", "")
        lines.append("| {tool} | {workflow} | {format} | {clusters} | {cycles} | {samples} | {requested_threads} | {gzip_level} | {seconds} | {clusters_per_sec} | {peak_rss_kb} | {output_mb} | {assigned_reads} | {undetermined_reads} | {filtered_clusters} | {validation_mismatches} | {validation_mode} | {content_hash} | {exit_code} | {logs} |".format(output_mb=output_mb, content_hash=content_hash, logs=logs, **row))
    lines.extend([
        "",
        "## Scope",
        "",
        "Run `make bcl-tiny-public-gate` to verify the narrow public 10x tiny-BCL classic per-cycle milestone. This check covers the committed DotMatch row, output hashes, count totals, and available bcl2fastq count-total validation.",
        "",
        "Broader raw-BCL evaluation needs real classic-BCL and CBCL run-folder rows, a successful DotMatch CBCL row, competitor rows for BCL Convert/bcl2fastq/CUDA-Demux where installable, distinct repeated timing, and `dotmatch bcl-validate` zero-mismatch evidence.",
        "",
        "Use `make bcl-comparison-gate` to check those rows when they are available.",
        "",
    ])
    (OUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print(OUT_DIR / "README.md")


if __name__ == "__main__":
    main()
