#!/usr/bin/env python3
"""Compare DotMatch count tables against CRISPR workflow competitors."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples" / "crispr_guides" / "output"


def public_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def parse_matrix(path: Path) -> tuple[list[str], dict[str, tuple[str, dict[str, int]]]]:
    """Read raw integer counts without coercion, silent drops or duplicate keys."""
    import re

    with path.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, [])
        if len(header) < 3 or len(set(header)) != len(header) or any(not value for value in header):
            raise ValueError(f"{path} does not have distinct ID, gene and sample columns")
        samples = header[2:]
        rows: dict[str, tuple[str, dict[str, int]]] = {}
        for row in reader:
            if not row:
                continue
            if len(row) != len(header) or not row[0]:
                raise ValueError(f"{path}: malformed count row at line {reader.line_num}")
            if row[0] in rows:
                raise ValueError(f"{path}: duplicate guide ID {row[0]!r}")
            values: dict[str, int] = {}
            for sample, text in zip(samples, row[2:]):
                if not re.fullmatch(r"[0-9]+", text):
                    raise ValueError(f"{path}: {row[0]}/{sample} must be a non-negative integer count, got {text!r}")
                values[sample] = int(text)
            rows[row[0]] = (row[1], values)
        return samples, rows


def parse_counts(path: Path) -> dict[str, int]:
    """Legacy guide totals; not a claim of sample-by-sample matrix identity."""
    _samples, rows = parse_matrix(path)
    return {guide: sum(values.values()) for guide, (_gene, values) in rows.items()}


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
            "left_path": public_path(left_path),
            "right_path": public_path(right_path),
            "n_guides": "0",
            "total_left": "",
            "total_right": "",
            "total_delta": "",
            "differing_guides": "",
            "max_abs_delta": "",
            "pearson": "",
            "spearman": "",
        }, []

    left_samples, left_matrix = parse_matrix(left_path)
    right_samples, right_matrix = parse_matrix(right_path)
    left = {guide: sum(values.values()) for guide, (_gene, values) in left_matrix.items()}
    right = {guide: sum(values.values()) for guide, (_gene, values) in right_matrix.items()}
    same_axes = set(left_samples) == set(right_samples) and left_matrix.keys() == right_matrix.keys()
    same_annotations = same_axes and all(left_matrix[guide][0] == right_matrix[guide][0] for guide in left_matrix)
    matrix_identical = same_annotations and all(left_matrix[guide][1] == right_matrix[guide][1] for guide in left_matrix)
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
        "status": "ok",  # Legacy execution status, not a count-identity verdict.
        "execution_status": "completed",
        "aggregate_scope": "guide_totals_across_samples",
        "aggregate_counts_identical": str(left == right).lower(),
        "matrix_comparability": "same_named_axes_and_annotations" if same_annotations else "review_axes_or_annotations",
        "counts_identical": str(matrix_identical).lower() if same_annotations else "",
        "left_path": public_path(left_path),
        "right_path": public_path(right_path),
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
