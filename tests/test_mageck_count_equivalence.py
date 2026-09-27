#!/usr/bin/env python3
"""Regression tests for complete MAGeCK-table comparison."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "compare_mageck_count_tables", ROOT / "scripts/compare_mageck_count_tables.py"
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("could not load MAGeCK count-table comparator")
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="dotmatch-mageck-table-") as raw_tmp:
        tmp = Path(raw_tmp)
        reference = tmp / "reference.txt"
        candidate = tmp / "candidate.txt"
        write(reference, "sgRNA\tGene\ta\tb\ng1\tGENE1\t1\t2\ng2\tGENE2\t0\t3\n")
        write(candidate, "Gene\tb\tsgRNA\ta\nGENE2\t3\tg2\t0\nGENE1\t2\tg1\t1\n")
        summary = MODULE.compare(reference, candidate, ["a", "b"])
        assert summary["status"] == "pass"
        assert summary["reference_rows"] == 2
        assert summary["compared_count_cells"] == 4

        failures = {
            "missing": "sgRNA\tGene\ta\tb\ng1\tGENE1\t1\t2\n",
            "extra": "sgRNA\tGene\ta\tb\ng1\tGENE1\t1\t2\ng2\tGENE2\t0\t3\ng3\tGENE3\t0\t0\n",
            "gene": "sgRNA\tGene\ta\tb\ng1\tWRONG\t1\t2\ng2\tGENE2\t0\t3\n",
            "count": "sgRNA\tGene\ta\tb\ng1\tGENE1\t1\t2\ng2\tGENE2\t0\t4\n",
        }
        for name, contents in failures.items():
            path = tmp / f"{name}.txt"
            write(path, contents)
            assert MODULE.compare(reference, path, ["a", "b"])["status"] == "fail"

        invalid = tmp / "invalid.txt"
        write(invalid, "sgRNA\tGene\ta\tb\ng1\tGENE1\t1\t2\ng1\tGENE1\t1\t2\n")
        try:
            MODULE.read_table(invalid, ["a", "b"])
        except ValueError as exc:
            assert "duplicate sgRNA" in str(exc)
        else:
            raise AssertionError("duplicate sgRNA was accepted")

    print("MAGeCK complete-table comparator: PASS")


if __name__ == "__main__":
    main()
