#!/usr/bin/env python3
"""Regression tests for MAGeCK downstream-result comparison."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "compare_mageck_test_results", ROOT / "scripts/compare_mageck_test_results.py"
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("could not load MAGeCK downstream comparator")
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)

GENE_HEADER = (
    "id\tnum\tneg|score\tneg|p-value\tneg|fdr\tneg|rank\tneg|goodsgrna\tneg|lfc\t"
    "pos|score\tpos|p-value\tpos|fdr\tpos|rank\tpos|goodsgrna\tpos|lfc\n"
)
SGRNA_HEADER = "sgrna\tGene\tcontrol_count\ttreatment_count\n"


def gene(identifier: str, neg_score: str, neg_rank: int, pos_score: str, pos_rank: int) -> str:
    return (
        f"{identifier}\t2\t{neg_score}\t{neg_score}\t{neg_score}\t{neg_rank}\t1\t-1.0\t"
        f"{pos_score}\t{pos_score}\t{pos_score}\t{pos_rank}\t1\t1.0\n"
    )


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def compare(tmp: Path, reference_genes: str, candidate_genes: str, candidate_sgrnas: str | None = None):
    reference_gene = write(tmp / "reference.gene.txt", GENE_HEADER + reference_genes)
    candidate_gene = write(tmp / "candidate.gene.txt", GENE_HEADER + candidate_genes)
    reference_sgrna = write(tmp / "reference.sgrna.txt", SGRNA_HEADER + "g1\tA\t1\t2\n")
    candidate_sgrna = write(
        tmp / "candidate.sgrna.txt",
        SGRNA_HEADER + (candidate_sgrnas if candidate_sgrnas is not None else "g1\tA\t1\t2\n"),
    )
    return MODULE.compare(reference_gene, candidate_gene, reference_sgrna, candidate_sgrna)


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="dotmatch-mageck-test-") as raw_tmp:
        tmp = Path(raw_tmp)
        exact = gene("A", "0.01", 1, "0.9", 2) + gene("B", "0.5", 2, "0.1", 1)
        summary = compare(tmp, exact, exact)
        assert summary["status"] == "pass_exact"

        tied_reference = gene("A", "0.01", 1, "0.9", 2) + gene("B", "0.01", 2, "0.9", 1)
        tied_candidate = gene("B", "0.01", 1, "0.9", 2) + gene("A", "0.01", 2, "0.9", 1)
        summary = compare(tmp, tied_reference, tied_candidate)
        assert summary["status"] == "pass_with_tied_rank_order_variation"
        assert summary["gene_summary"]["directions"]["neg"]["rank_mismatches"] == 2
        assert summary["gene_summary"]["directions"]["neg"]["top_n"]["10"]["symmetric_difference"] == 0

        boundary_reference = "".join(
            gene(chr(65 + index), "0.01", index + 1, "0.9", index + 1)
            for index in range(11)
        )
        boundary_candidate = "".join(
            gene(chr(65 + index), "0.01", 11 if index == 9 else 10 if index == 10 else index + 1,
                 "0.9", index + 1)
            for index in range(11)
        )
        summary = compare(tmp, boundary_reference, boundary_candidate)
        assert summary["status"] == "pass_with_tied_rank_order_variation"
        assert summary["gene_summary"]["directions"]["neg"]["top_n"]["10"]["symmetric_difference"] == 2

        crossed_candidate = gene("A", "0.01", 2, "0.9", 2) + gene("B", "0.5", 1, "0.1", 1)
        summary = compare(tmp, exact, crossed_candidate)
        assert summary["status"] == "fail"
        assert summary["gene_summary"]["directions"]["neg"]["rank_changes_crossing_tie_blocks"] == 2

        changed_stat = gene("A", "0.02", 1, "0.9", 2) + gene("B", "0.5", 2, "0.1", 1)
        summary = compare(tmp, exact, changed_stat)
        assert summary["status"] == "fail"
        assert summary["gene_summary"]["non_rank_value_mismatches"] == 3

        summary = compare(tmp, exact, exact, candidate_sgrnas="g1\tA\t1\t3\n")
        assert summary["status"] == "fail"
        assert summary["sgrna_summary"]["value_mismatches"] == 1

        duplicate = write(tmp / "duplicate.gene.txt", GENE_HEADER + gene("A", "0.01", 1, "0.9", 1) * 2)
        try:
            MODULE.read_table(duplicate, "id", MODULE.GENE_REQUIRED)
        except ValueError as exc:
            assert "duplicate id" in str(exc)
        else:
            raise AssertionError("duplicate gene ID was accepted")

        malformed = write(tmp / "malformed.gene.txt", GENE_HEADER + "A\t2\t0.1\n")
        try:
            MODULE.read_table(malformed, "id", MODULE.GENE_REQUIRED)
        except ValueError as exc:
            assert "missing fields" in str(exc)
        else:
            raise AssertionError("truncated MAGeCK row was accepted")

    print("MAGeCK downstream comparator: PASS")


if __name__ == "__main__":
    main()
