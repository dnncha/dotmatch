#!/usr/bin/env python3
"""Compare complete MAGeCK-format count tables without depending on row order."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_table(path: Path, samples: list[str]) -> dict[str, tuple[str, tuple[int, ...]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = ["sgRNA", "Gene", *samples]
        missing = [name for name in required if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(missing)}")
        rows: dict[str, tuple[str, tuple[int, ...]]] = {}
        for line_number, row in enumerate(reader, start=2):
            guide = row["sgRNA"]
            if not guide:
                raise ValueError(f"{path}:{line_number}: empty sgRNA")
            if guide in rows:
                raise ValueError(f"{path}:{line_number}: duplicate sgRNA: {guide}")
            counts: list[int] = []
            for sample in samples:
                try:
                    count = int(row[sample])
                except ValueError as exc:
                    raise ValueError(
                        f"{path}:{line_number}: non-integer count for {sample}: {row[sample]!r}"
                    ) from exc
                if count < 0:
                    raise ValueError(f"{path}:{line_number}: negative count for {sample}: {count}")
                counts.append(count)
            rows[guide] = (row["Gene"], tuple(counts))
    return rows


def compare(reference_path: Path, candidate_path: Path, samples: list[str]) -> dict[str, object]:
    reference = read_table(reference_path, samples)
    candidate = read_table(candidate_path, samples)
    reference_ids = set(reference)
    candidate_ids = set(candidate)
    common_ids = reference_ids & candidate_ids
    gene_mismatches = sum(reference[guide][0] != candidate[guide][0] for guide in common_ids)
    count_mismatches = sum(
        expected != observed
        for guide in common_ids
        for expected, observed in zip(reference[guide][1], candidate[guide][1])
    )
    missing_ids = sorted(reference_ids - candidate_ids)
    extra_ids = sorted(candidate_ids - reference_ids)
    status = "pass" if not (missing_ids or extra_ids or gene_mismatches or count_mismatches) else "fail"
    return {
        "status": status,
        "samples": samples,
        "reference_rows": len(reference),
        "candidate_rows": len(candidate),
        "compared_count_cells": len(common_ids) * len(samples),
        "missing_guides": len(missing_ids),
        "extra_guides": len(extra_ids),
        "gene_mismatches": gene_mismatches,
        "count_mismatches": count_mismatches,
        "first_missing_guides": missing_ids[:10],
        "first_extra_guides": extra_ids[:10],
        "reference_sha256": file_sha256(reference_path),
        "candidate_sha256": file_sha256(candidate_path),
        "comparison_semantics": "sgRNA-keyed; exact Gene strings; integer counts; row-order independent",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path, help="Reference MAGeCK TSV")
    parser.add_argument("candidate", type=Path, help="Candidate MAGeCK TSV")
    parser.add_argument("--samples", required=True, help="Comma-separated sample columns to compare")
    parser.add_argument("--summary-out", type=Path)
    args = parser.parse_args()
    samples = [sample.strip() for sample in args.samples.split(",") if sample.strip()]
    if not samples or len(samples) != len(set(samples)):
        raise SystemExit("--samples must contain unique, non-empty sample names")
    try:
        summary = compare(args.reference, args.candidate, samples)
    except (OSError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc
    rendered = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.summary_out:
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        args.summary_out.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    if summary["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
