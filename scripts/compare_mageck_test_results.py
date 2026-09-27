#!/usr/bin/env python3
"""Compare MAGeCK test outputs while making tied-rank ordering explicit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


RANK_FIELDS = ("neg|rank", "pos|rank")
GENE_REQUIRED = {
    "id",
    "num",
    "neg|score",
    "neg|p-value",
    "neg|fdr",
    "neg|rank",
    "neg|goodsgrna",
    "neg|lfc",
    "pos|score",
    "pos|p-value",
    "pos|fdr",
    "pos|rank",
    "pos|goodsgrna",
    "pos|lfc",
}
SGRNA_REQUIRED = {"sgrna", "Gene"}
TOP_N_CUTOFFS = (10, 50, 100, 500, 1000, 5000, 10000)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_table(path: Path, key: str, required: set[str]) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        fields = reader.fieldnames or []
        missing = sorted(required - set(fields))
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(missing)}")
        rows: dict[str, dict[str, str]] = {}
        for line_number, row in enumerate(reader, start=2):
            identifier = row[key]
            if not identifier:
                raise ValueError(f"{path}:{line_number}: empty {key}")
            if identifier in rows:
                raise ValueError(f"{path}:{line_number}: duplicate {key}: {identifier}")
            if None in row:
                raise ValueError(f"{path}:{line_number}: extra fields")
            if any(value is None for value in row.values()):
                raise ValueError(f"{path}:{line_number}: missing fields")
            rows[identifier] = row
    return fields, rows


def _membership(reference: dict[str, object], candidate: dict[str, object]) -> tuple[list[str], list[str], list[str]]:
    reference_ids, candidate_ids = set(reference), set(candidate)
    return (
        sorted(reference_ids & candidate_ids),
        sorted(reference_ids - candidate_ids),
        sorted(candidate_ids - reference_ids),
    )


def _compare_sgrna(reference_path: Path, candidate_path: Path) -> dict[str, object]:
    reference_fields, reference = read_table(reference_path, "sgrna", SGRNA_REQUIRED)
    candidate_fields, candidate = read_table(candidate_path, "sgrna", SGRNA_REQUIRED)
    common, missing, extra = _membership(reference, candidate)
    header_match = reference_fields == candidate_fields
    fields = reference_fields if header_match else []
    mismatches = sum(
        reference[identifier][field] != candidate[identifier][field]
        for identifier in common
        for field in fields
    )
    return {
        "status": "pass" if header_match and not (missing or extra or mismatches) else "fail",
        "reference_rows": len(reference),
        "candidate_rows": len(candidate),
        "fields_per_row": len(fields),
        "compared_cells": len(common) * len(fields),
        "header_match": header_match,
        "missing_sgrnas": len(missing),
        "extra_sgrnas": len(extra),
        "value_mismatches": mismatches,
        "first_missing_sgrnas": missing[:10],
        "first_extra_sgrnas": extra[:10],
        "reference_sha256": file_sha256(reference_path),
        "candidate_sha256": file_sha256(candidate_path),
    }


def _tie_groups(rows: dict[str, dict[str, str]], direction: str) -> dict[tuple[str, str, str], list[str]]:
    groups: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    for identifier, row in rows.items():
        groups[(row[f"{direction}|score"], row[f"{direction}|p-value"], row[f"{direction}|fdr"])].append(identifier)
    return groups


def _compare_direction(
    reference: dict[str, dict[str, str]],
    candidate: dict[str, dict[str, str]],
    common: list[str],
    direction: str,
    fdr_threshold: float,
) -> dict[str, object]:
    rank_field = f"{direction}|rank"
    rank_changes = [
        identifier
        for identifier in common
        if reference[identifier][rank_field] != candidate[identifier][rank_field]
    ]
    rank_change_ids = set(rank_changes)
    invalid_rank_values: list[str] = []
    rank_delta: list[int] = []
    for identifier in rank_changes:
        try:
            rank_delta.append(int(candidate[identifier][rank_field]) - int(reference[identifier][rank_field]))
        except ValueError:
            invalid_rank_values.append(identifier)

    reference_groups = _tie_groups(reference, direction)
    candidate_groups = _tie_groups(candidate, direction)
    crossing_blocks: list[dict[str, object]] = []
    changed_blocks = 0
    for tie_key, identifiers in reference_groups.items():
        changed = [identifier for identifier in identifiers if identifier in rank_change_ids]
        if not changed:
            continue
        changed_blocks += 1
        reference_ranks = sorted(
            (reference[identifier][rank_field] for identifier in identifiers), key=int
        )
        candidate_ids = candidate_groups.get(tie_key, [])
        candidate_ranks = sorted(
            (candidate[identifier][rank_field] for identifier in candidate_ids), key=int
        )
        if set(identifiers) != set(candidate_ids) or reference_ranks != candidate_ranks:
            crossing_blocks.append(
                {
                    "score": tie_key[0],
                    "p_value": tie_key[1],
                    "fdr": tie_key[2],
                    "reference_genes": len(identifiers),
                    "candidate_genes": len(candidate_ids),
                }
            )

    def hit_set(rows: dict[str, dict[str, str]]) -> set[str]:
        hits: set[str] = set()
        for identifier, row in rows.items():
            try:
                if float(row[f"{direction}|fdr"]) <= fdr_threshold:
                    hits.add(identifier)
            except ValueError as exc:
                raise ValueError(f"invalid {direction}|fdr for {identifier}: {row[f'{direction}|fdr']!r}") from exc
        return hits

    reference_hits, candidate_hits = hit_set(reference), hit_set(candidate)
    top_n: dict[str, object] = {}
    for cutoff in TOP_N_CUTOFFS:
        reference_top = {identifier for identifier, row in reference.items() if int(row[rank_field]) <= cutoff}
        candidate_top = {identifier for identifier, row in candidate.items() if int(row[rank_field]) <= cutoff}
        top_n[str(cutoff)] = {
            "reference_genes": len(reference_top),
            "candidate_genes": len(candidate_top),
            "symmetric_difference": len(reference_top ^ candidate_top),
            "first_reference_only": sorted(reference_top - candidate_top)[:10],
            "first_candidate_only": sorted(candidate_top - reference_top)[:10],
        }

    return {
        "rank_mismatches": len(rank_changes),
        "genes_with_invalid_rank": invalid_rank_values[:10],
        "max_absolute_rank_delta": max(map(abs, rank_delta), default=0),
        "tie_blocks_with_rank_variation": changed_blocks,
        "rank_changes_crossing_tie_blocks": len(crossing_blocks),
        "first_crossing_tie_blocks": crossing_blocks[:10],
        "fdr_threshold": fdr_threshold,
        "reference_fdr_hits": len(reference_hits),
        "candidate_fdr_hits": len(candidate_hits),
        "missing_fdr_hits": sorted(reference_hits - candidate_hits)[:10],
        "extra_fdr_hits": sorted(candidate_hits - reference_hits)[:10],
        "fdr_hit_sets_identical": reference_hits == candidate_hits,
        "top_n": top_n,
    }


def _compare_genes(reference_path: Path, candidate_path: Path, fdr_threshold: float) -> dict[str, object]:
    reference_fields, reference = read_table(reference_path, "id", GENE_REQUIRED)
    candidate_fields, candidate = read_table(candidate_path, "id", GENE_REQUIRED)
    common, missing, extra = _membership(reference, candidate)
    header_match = reference_fields == candidate_fields
    fields = reference_fields if header_match else []
    non_rank_fields = [field for field in fields if field not in RANK_FIELDS]
    non_rank_mismatches = sum(
        reference[identifier][field] != candidate[identifier][field]
        for identifier in common
        for field in non_rank_fields
    )
    directions = {
        direction: _compare_direction(reference, candidate, common, direction, fdr_threshold)
        for direction in ("neg", "pos")
    }
    rank_mismatches = sum(item["rank_mismatches"] for item in directions.values())
    substantive_failure = bool(
        not header_match
        or missing
        or extra
        or non_rank_mismatches
        or any(item["genes_with_invalid_rank"] for item in directions.values())
        or any(item["rank_changes_crossing_tie_blocks"] for item in directions.values())
        or any(not item["fdr_hit_sets_identical"] for item in directions.values())
    )
    status = "fail" if substantive_failure else (
        "pass_with_tied_rank_order_variation" if rank_mismatches else "pass_exact"
    )
    return {
        "status": status,
        "reference_rows": len(reference),
        "candidate_rows": len(candidate),
        "fields_per_row": len(fields),
        "compared_non_rank_cells": len(common) * len(non_rank_fields),
        "header_match": header_match,
        "missing_genes": len(missing),
        "extra_genes": len(extra),
        "non_rank_value_mismatches": non_rank_mismatches,
        "total_rank_mismatches": rank_mismatches,
        "first_missing_genes": missing[:10],
        "first_extra_genes": extra[:10],
        "directions": directions,
        "reference_sha256": file_sha256(reference_path),
        "candidate_sha256": file_sha256(candidate_path),
    }


def compare(
    reference_gene: Path,
    candidate_gene: Path,
    reference_sgrna: Path,
    candidate_sgrna: Path,
    fdr_threshold: float = 0.05,
) -> dict[str, object]:
    if not 0.0 <= fdr_threshold <= 1.0:
        raise ValueError("FDR threshold must be between 0 and 1")
    genes = _compare_genes(reference_gene, candidate_gene, fdr_threshold)
    sgrnas = _compare_sgrna(reference_sgrna, candidate_sgrna)
    status = "fail" if genes["status"] == "fail" or sgrnas["status"] == "fail" else genes["status"]
    return {
        "schema_version": 1,
        "status": status,
        "gene_summary": genes,
        "sgrna_summary": sgrnas,
        "comparison_semantics": (
            "identifier-keyed exact strings; all sgRNA fields and all non-rank gene fields must match; "
            "gene ranks may vary only by permutation within identical score/p-value/FDR tie blocks"
        ),
        "claim_boundary": (
            "A passing result establishes output equivalence under the stated tie semantics. "
            "It does not establish that MAGeCK statistics or either counting workflow are biologically correct."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_gene", type=Path)
    parser.add_argument("candidate_gene", type=Path)
    parser.add_argument("reference_sgrna", type=Path)
    parser.add_argument("candidate_sgrna", type=Path)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--summary-out", type=Path)
    args = parser.parse_args()
    try:
        summary = compare(
            args.reference_gene,
            args.candidate_gene,
            args.reference_sgrna,
            args.candidate_sgrna,
            args.fdr_threshold,
        )
    except (OSError, ValueError, csv.Error) as exc:
        raise SystemExit(str(exc)) from exc
    rendered = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.summary_out:
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        args.summary_out.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    if summary["status"] == "fail":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
