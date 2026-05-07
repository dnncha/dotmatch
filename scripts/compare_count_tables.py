#!/usr/bin/env python3
"""Compare DotMatch count tables against CRISPR workflow competitors."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples" / "crispr_guides" / "output"


def parse_counts(path: Path) -> dict[str, int]:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or len(reader.fieldnames) < 3:
            raise ValueError(f"{path} does not look like a count table")
        key_col = reader.fieldnames[0]
        count_cols = reader.fieldnames[2:]
        counts: dict[str, int] = {}
        for row in reader:
            key = row.get(key_col, "")
            if not key:
                continue
            total = 0
            for col in count_cols:
                value = row.get(col, "") or "0"
                try:
                    total += int(float(value))
                except ValueError:
                    pass
            counts[key] = total
    return counts


def pearson(xs: list[int], ys: list[int]) -> float:
    if len(xs) < 2:
        return float("nan")
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if den_x == 0 or den_y == 0:
        return float("nan")
    return num / (den_x * den_y)


def ranks(values: list[int]) -> list[float]:
    order = sorted(enumerate(values), key=lambda item: item[1])
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and order[j][1] == order[i][1]:
            j += 1
        rank = (i + 1 + j) / 2.0
        for idx in range(i, j):
            out[order[idx][0]] = rank
        i = j
    return out


def spearman(xs: list[int], ys: list[int]) -> float:
    return pearson([int(r * 1000000) for r in ranks(xs)], [int(r * 1000000) for r in ranks(ys)])


def compare(name: str, left_path: Path, right_path: Path, left_label: str, right_label: str) -> tuple[dict[str, str], list[dict[str, str]]]:
    if not left_path.exists() or not right_path.exists():
        return {
            "comparison": name,
            "left": left_label,
            "right": right_label,
            "status": "missing_input",
            "left_path": str(left_path),
            "right_path": str(right_path),
            "n_guides": "0",
            "total_left": "",
            "total_right": "",
            "total_delta": "",
            "differing_guides": "",
            "max_abs_delta": "",
            "pearson": "",
            "spearman": "",
        }, []

    left = parse_counts(left_path)
    right = parse_counts(right_path)
    keys = sorted(set(left) | set(right))
    detail: list[dict[str, str]] = []
    xs: list[int] = []
    ys: list[int] = []
    max_abs_delta = 0
    differing = 0
    for key in keys:
        lval = left.get(key, 0)
        rval = right.get(key, 0)
        delta = lval - rval
        xs.append(lval)
        ys.append(rval)
        if delta != 0:
            differing += 1
        max_abs_delta = max(max_abs_delta, abs(delta))
        detail.append({
            "comparison": name,
            "guide_id": key,
            left_label: str(lval),
            right_label: str(rval),
            "delta": str(delta),
            "abs_delta": str(abs(delta)),
        })
    detail.sort(key=lambda row: int(row["abs_delta"]), reverse=True)
    total_left = sum(xs)
    total_right = sum(ys)
    summary = {
        "comparison": name,
        "left": left_label,
        "right": right_label,
        "status": "ok",
        "left_path": str(left_path),
        "right_path": str(right_path),
        "n_guides": str(len(keys)),
        "total_left": str(total_left),
        "total_right": str(total_right),
        "total_delta": str(total_left - total_right),
        "differing_guides": str(differing),
        "max_abs_delta": str(max_abs_delta),
        "pearson": f"{pearson(xs, ys):.8f}",
        "spearman": f"{spearman(xs, ys):.8f}",
    }
    return summary, detail


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
    parser.add_argument("--dotmatch-hamming", default=str(OUT / "counts.hamming.mageck.tsv"))
    parser.add_argument("--guide-counter", default=str(OUT / "guide_counter.counts.txt"))
    parser.add_argument("--dotmatch-exact", default=str(OUT / "counts.exact.mageck.tsv"))
    parser.add_argument("--mageck-exact", default=str(OUT / "mageck_exact_benchmark.count.txt"))
    parser.add_argument("--summary-out", default=str(ROOT / "benchmarks" / "raw" / "count_agreement_summary.csv"))
    parser.add_argument("--details-out", default=str(ROOT / "benchmarks" / "raw" / "count_agreement_details.csv"))
    args = parser.parse_args()

    comparisons = [
        ("dotmatch_hamming_vs_guide_counter", Path(args.dotmatch_hamming), Path(args.guide_counter), "dotmatch_hamming", "guide_counter"),
        ("dotmatch_exact_vs_mageck_exact", Path(args.dotmatch_exact), Path(args.mageck_exact), "dotmatch_exact", "mageck_exact"),
    ]
    summaries: list[dict[str, str]] = []
    details: list[dict[str, str]] = []
    for item in comparisons:
        summary, detail = compare(*item)
        summaries.append(summary)
        details.extend(detail[:50])
    write_csv(Path(args.summary_out), summaries)
    write_csv(Path(args.details_out), details or [{"comparison": "", "guide_id": "", "delta": ""}])
    print(args.summary_out)
    print(args.details_out)


if __name__ == "__main__":
    main()
