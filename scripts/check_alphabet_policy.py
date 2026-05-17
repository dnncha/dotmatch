#!/usr/bin/env python3
"""Audit the DotMatch alphabet policy contract and benchmark provenance."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


POLICY = "literal-byte; A/C/G/T/N/IUPAC symbols are ordinary byte symbols; no wildcard expansion"


class AuditResult:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _header_policy(text: str) -> str:
    match = re.search(r'#define\s+QDALN_ALPHABET_POLICY\s+"([^"]+)"', text)
    return match.group(1) if match else ""


def _check_core(root: Path, result: AuditResult) -> None:
    header_path = root / "include" / "qdalign.h"
    source_path = root / "src" / "qdalign.c"
    try:
        header = _read(header_path)
        source = _read(source_path)
    except Exception as exc:
        result.failures.append(f"alphabet policy core files could not be read: {exc}")
        return

    declared = _header_policy(header)
    if declared != POLICY:
        result.failures.append("QDALN_ALPHABET_POLICY must declare the literal-byte N/IUPAC contract")
    if "return QDALN_ALPHABET_POLICY;" not in source:
        result.failures.append("qdaln_alphabet_policy must return QDALN_ALPHABET_POLICY")


def _is_dotmatch_tool(value: str) -> bool:
    return value.startswith("dotmatch")


def _check_raw_csvs(root: Path, result: AuditResult) -> None:
    raw_dir = root / "benchmarks" / "raw"
    if not raw_dir.exists():
        return
    for path in sorted(raw_dir.glob("*.csv")):
        try:
            with path.open(encoding="utf-8", newline="") as fh:
                reader = csv.DictReader(fh)
                if "alphabet_policy" not in (reader.fieldnames or []):
                    continue
                for line_no, row in enumerate(reader, start=2):
                    tool = str(row.get("tool") or "")
                    if not _is_dotmatch_tool(tool):
                        continue
                    policy = str(row.get("alphabet_policy") or "")
                    if not policy:
                        result.failures.append(
                            f"{tool} row in {path.relative_to(root).as_posix()} must record alphabet_policy"
                        )
                    elif policy != POLICY:
                        result.failures.append(
                            f"{tool} row in {path.relative_to(root).as_posix()}:{line_no} "
                            "must use literal-byte alphabet_policy"
                        )
        except Exception as exc:
            result.failures.append(f"{path.relative_to(root).as_posix()} could not be read: {exc}")


def audit(root: Path) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    _check_core(root, result)
    _check_raw_csvs(root, result)
    if result.ok:
        result.passed.append("alphabet policy contract exported")
        result.passed.append("dotmatch benchmark rows record literal-byte alphabet policy")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="repository root")
    args = parser.parse_args()

    result = audit(Path(args.root))
    for item in result.passed:
        print(f"PASS: {item}")
    for item in result.failures:
        print(f"FAIL: {item}")
    if result.ok:
        print("ALPHABET POLICY: PASS")
        return 0
    print("ALPHABET POLICY: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
