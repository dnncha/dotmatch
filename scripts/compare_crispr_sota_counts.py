#!/usr/bin/env python3
"""Compare CRISPR SOTA DotMatch count tables against fair competitors."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from compare_count_tables import ROOT, compare


RAW = ROOT / "benchmarks" / "raw"


def dataset_paths(dataset: str) -> dict[str, Path]:
    if dataset == "mageck_yusa":
        out = ROOT / "examples" / "crispr_guides" / "output"
        return {
            "dotmatch_hamming": out / "counts.hamming.mageck.tsv",
            "guide_counter": out / "guide_counter.counts.txt",
            "dotmatch_exact": out / "counts.exact.mageck.tsv",
            "mageck_exact": out / "mageck_exact_benchmark.count.txt",
        }
    out = ROOT / "examples" / dataset / "output"
    return {
        "dotmatch_hamming": out / "counts.hamming.mageck.tsv",
        "guide_counter": out / "guide_counter.counts.txt",
        "dotmatch_exact": out / "counts.exact.mageck.tsv",
        "mageck_exact": out / f"{dataset}_mageck_exact.count.txt",
    }


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", default="mageck_yusa,sanson_brunello")
    parser.add_argument("--summary-out", default=str(RAW / "crispr_sota_count_agreement_summary.csv"))
    parser.add_argument("--details-out", default=str(RAW / "crispr_sota_count_agreement_details.csv"))
    args = parser.parse_args()

    summaries: list[dict[str, str]] = []
    details: list[dict[str, str]] = []
    for dataset in [d.strip() for d in args.datasets.split(",") if d.strip()]:
        paths = dataset_paths(dataset)
        for name, left, right, left_label, right_label in [
            ("dotmatch_hamming_vs_guide_counter", paths["dotmatch_hamming"], paths["guide_counter"],
             "dotmatch_hamming", "guide_counter"),
            ("dotmatch_exact_vs_mageck_exact", paths["dotmatch_exact"], paths["mageck_exact"],
             "dotmatch_exact", "mageck_exact"),
        ]:
            summary, detail = compare(f"{dataset}:{name}", left, right, left_label, right_label)
            summary["dataset"] = dataset
            summaries.append(summary)
            for row in detail[:50]:
                row["dataset"] = dataset
                details.append(row)
    write_csv(Path(args.summary_out), summaries)
    write_csv(Path(args.details_out), details or [{"dataset": "", "comparison": "", "guide_id": "", "delta": ""}])
    print(args.summary_out)
    print(args.details_out)


if __name__ == "__main__":
    main()
