#!/usr/bin/env python3
"""Check that the barcode-autopsy demo has explicit failure-mode fixtures."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REQUIRED_FAILURE_MODES = {
    "wrong_offset",
    "duplicate_barcode",
    "unsafe_one_edit_collision",
    "ambiguous_read",
    "unmatched_low_complexity",
    "low_quality_correction_rejected",
    "invalid_window",
    "reverse_complement_candidate",
}


class AuditResult:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def _read_findings(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def audit(root: Path) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    fixture_dir = root / "examples" / "barcode_autopsy" / "failure_modes"
    required_files = ["barcodes.tsv", "reads.fastq", "expected_findings.tsv", "README.md"]
    for name in required_files:
        path = fixture_dir / name
        if not path.is_file():
            result.failures.append(f"missing barcode failure fixture file: {path.relative_to(root) if path.is_absolute() else path}")

    findings_path = fixture_dir / "expected_findings.tsv"
    if findings_path.is_file():
        findings = _read_findings(findings_path)
        if not findings:
            result.failures.append("expected_findings.tsv must contain at least one finding")
        required_columns = {"failure_mode", "read_id", "observed", "diagnosis", "next_action"}
        present_columns = set(findings[0]) if findings else set()
        missing_columns = required_columns - present_columns
        for column in sorted(missing_columns):
            result.failures.append(f"expected_findings.tsv missing required column: {column}")
        modes = {row.get("failure_mode", "") for row in findings}
        for mode in sorted(REQUIRED_FAILURE_MODES - modes):
            result.failures.append(f"missing required failure mode: {mode}")
        for index, row in enumerate(findings, start=2):
            for column in required_columns:
                if not str(row.get(column) or "").strip():
                    result.failures.append(f"expected_findings.tsv:{index} has empty {column}")

    if result.ok:
        result.passed.append(f"{len(REQUIRED_FAILURE_MODES)} barcode failure modes documented")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="repository root")
    args = parser.parse_args()
    result = audit(Path(args.root))
    for item in result.passed:
        print(f"PASS: {item}")
    for failure in result.failures:
        print(f"FAIL: {failure}")
    if result.ok:
        print("BARCODE FAILURE FIXTURES: PASS")
        return 0
    print("BARCODE FAILURE FIXTURES: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
