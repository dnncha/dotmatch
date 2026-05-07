#!/usr/bin/env python3
"""Generate barcode-demultiplexing benchmark figures and report."""

from __future__ import annotations

import csv
import statistics
from pathlib import Path
import os


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "barcode_demux.csv"
OUT_DIR = ROOT / "docs" / "benchmarks" / "barcode_demux"
FIG_DIR = ROOT / "benchmarks" / "figures"


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh))


def svg_bar(rows: list[dict[str, str]], value_key: str, title: str, subtitle: str,
            path: Path, suffix: str = "", color: str = "#2f7d68") -> None:
    selected = [r for r in rows if r.get("exit_code") == "0" and r.get(value_key)]
    if not selected:
        return
    labels = [r["tool"] for r in selected]
    values = [float(r.get(value_key) or 0.0) for r in selected]
    width = 980
    left = 280
    row_h = 36
    height = 76 + row_h * len(selected)
    max_v = max(values) if values else 1.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<style>text{{font-family:Arial,sans-serif;font-size:13px}}.title{{font-size:18px;font-weight:700}}.axis{{fill:#444}}.bar{{fill:{color}}}</style>',
        f'<text class="title" x="20" y="28">{title}</text>',
        f'<text class="axis" x="20" y="50">{subtitle}</text>',
    ]
    for i, (label, value) in enumerate(zip(labels, values)):
        y = 80 + i * row_h
        bar_w = 1 if max_v == 0 else int((width - left - 110) * value / max_v)
        parts.append(f'<text x="20" y="{y + 15}">{label}</text>')
        parts.append(f'<rect class="bar" x="{left}" y="{y}" width="{bar_w}" height="21" rx="2"/>')
        parts.append(f'<text x="{left + bar_w + 8}" y="{y + 15}">{value:.2f}{suffix}</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def aggregate_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for row in rows:
        key = (row.get("tool", ""), row.get("workflow", ""), row.get("semantics", ""))
        groups.setdefault(key, []).append(row)
    out: list[dict[str, str]] = []
    for (tool, workflow, semantics), group in sorted(groups.items()):
        ok = [r for r in group if r.get("exit_code") == "0"]
        base = ok[0] if ok else group[0]
        values = [float(r.get("reads_per_sec") or 0.0) for r in ok]
        seconds = [float(r.get("seconds") or 0.0) for r in ok]
        rss = [float(r["peak_rss_kb"]) for r in ok if r.get("peak_rss_kb")]
        mean = statistics.mean(values) if values else 0.0
        stdev = statistics.stdev(values) if len(values) > 1 else 0.0
        merged = dict(base)
        merged.update({
            "tool": tool,
            "workflow": workflow,
            "semantics": semantics,
            "repeat": str(len(group)),
            "seconds": f"{statistics.mean(seconds):.6f}" if seconds else "0.000000",
            "reads_per_sec": f"{mean:.1f}",
            "peak_rss_kb": f"{max(rss):.0f}" if rss else "",
            "cv": f"{(stdev / mean):.4f}" if mean else "0.0000",
            "exit_code": "0" if len(ok) == len(group) else "1",
        })
        out.append(merged)
    return out


def write_report(rows: list[dict[str, str]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    throughput = FIG_DIR / "barcode_demux_throughput.svg"
    memory = FIG_DIR / "barcode_demux_peak_memory.svg"
    assigned = FIG_DIR / "barcode_demux_assigned_reads.svg"
    verified = FIG_DIR / "barcode_demux_verified_per_read.svg"
    summary_rows = aggregate_rows(rows)
    svg_bar(summary_rows, "reads_per_sec", "Barcode demultiplexing throughput",
            "Inline FASTQ barcode demux; higher is better", throughput, " reads/s")
    svg_bar(summary_rows, "peak_rss_kb", "Barcode demultiplexing peak memory",
            "Peak resident memory; lower is better", memory, " KB", "#a45f3f")
    svg_bar(summary_rows, "assigned_reads", "Barcode demultiplexing assigned reads",
            "Reads assigned to a barcode-specific output FASTQ", assigned, " reads", "#476fb3")
    svg_bar([r for r in summary_rows if r.get("verified_per_read")], "verified_per_read",
            "DotMatch verified candidates per read",
            "Candidate verification work; lower is better", verified, "", "#7f5aa2")

    lines = [
        "# Barcode Demultiplexing Benchmark",
        "",
        "This report is the barcode-demultiplexing evidence track. It is separate from the CRISPR guide-counting report.",
        "",
        "Current status: DotMatch now has a native `demux` command for fixed-position inline barcodes. A state-of-the-art barcode claim requires real public barcode datasets plus competitor rows, not only the built-in fixture.",
        "",
        "The benchmark script can also emit a simple `hash_splitter_exact` row. This is a transparent exact-prefix baseline, not an edit-distance demultiplexer.",
        "",
        "## Figures",
        "",
    ]
    for label, path in [
        ("Throughput", throughput),
        ("Peak memory", memory),
        ("Assigned reads", assigned),
        ("Verified candidates/read", verified),
    ]:
        if path.exists():
            rel = os.path.relpath(path, OUT_DIR)
            lines.append(f"![{label}]({rel})")
            lines.append("")
    lines.extend([
        "## Raw Rows",
        "",
        "| tool | workflow | semantics | repeats | reads | barcodes | k | metric | mean seconds | mean reads/sec | peak RSS KB | assigned | ambiguous | unmatched | verified/read | cv | exit |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for row in summary_rows:
        lines.append(
            "| {tool} | {workflow} | {semantics} | {repeat} | {n_reads} | {n_barcodes} | {k} | {metric} | {seconds} | {reads_per_sec} | {peak_rss_kb} | {assigned_reads} | {ambiguous_reads} | {unmatched_reads} | {verified_per_read} | {cv} | {exit_code} |".format(**row)
        )
    lines.extend([
        "",
        "## Publication Gate",
        "",
        "Do not describe DotMatch as barcode state of the art until this table includes real public barcode workloads and fair competitor rows, at minimum Cutadapt plus a second comparator such as Ultraplex, Je, deML, sabre/fastx-style splitters, an exact hash splitter for the exact-prefix lane, and Illumina demux tools where their input model matches the benchmark.",
        "",
        "Suggested real-data starting point: SRP009896 / SRR391079-SRR391082, a maize GBS dataset described in public Cutadapt demultiplexing examples as 5-prime inline barcode reads with 96 demultiplexed outputs. `scripts/fetch_srp009896_barcode_demo.py --use-public-example-barcodes` extracts the first-member barcode sheet from the public Google Drive example archive with a ranged request instead of downloading the full 7.4 GB ZIP.",
        "",
        "Important boundary: the SRP009896 barcode sheet contains variable-length barcodes (`4-8 bp`), while the current DotMatch `demux` benchmark is fixed-length. Treat any fixed-length subset result as a subset claim unless DotMatch adds variable-length inline demux semantics.",
        "",
    ])
    (OUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = read_rows(RAW)
    if not rows:
        raise SystemExit(f"missing benchmark CSV: {RAW}")
    write_report(rows)
    print(OUT_DIR / "README.md")


if __name__ == "__main__":
    main()
