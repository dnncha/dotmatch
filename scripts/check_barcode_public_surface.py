#!/usr/bin/env python3
"""Check that barcode-autopsy docs are easy to use and conservatively worded."""

from __future__ import annotations

import argparse
from pathlib import Path


REQUIRED_README_PHRASES = [
    "dotmatch barcode autopsy",
    "make barcode-science-ready",
    "See [Barcode Science Readiness]",
    "Speed claims are",
]
REQUIRED_SCIENCE_PHRASES = [
    "make barcode-science-ready",
    "at least five public fixed-window evidence datasets",
    "explicit failure-mode fixtures",
    "These datasets are not interchangeable biological claims",
]
REQUIRED_DEMO_PHRASES = [
    "make barcode-autopsy-demo",
    "report.html",
    "offset_scan.tsv",
    "provenance.json",
]
BLOCKED_PHRASES = [
    "dominate",
    "will replace cutadapt",
    "replaces cutadapt",
    "will replace bcl convert",
    "replaces bcl convert",
    "100% scientifically accurate",
    "best demux tool",
    "clinical-grade",
]


class AuditResult:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def _read(path: Path, result: AuditResult) -> str:
    if not path.is_file():
        result.failures.append(f"missing required public barcode doc: {path}")
        return ""
    return path.read_text(encoding="utf-8")


def _require(text: str, phrase: str, label: str, result: AuditResult) -> None:
    if phrase not in text:
        result.failures.append(f"{label} must include: {phrase}")


def _reject_hype(text: str, label: str, result: AuditResult) -> None:
    lower = text.lower()
    for phrase in BLOCKED_PHRASES:
        if phrase in lower:
            result.failures.append(f"{label} contains overbroad or hype wording: {phrase}")


def audit(root: Path) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    readme = _read(root / "README.md", result)
    science = _read(root / "docs" / "barcode-science-readiness.md", result)
    demo = _read(root / "examples" / "barcode_autopsy" / "README.md", result)
    for phrase in REQUIRED_README_PHRASES:
        _require(readme, phrase, "README.md", result)
    for phrase in REQUIRED_SCIENCE_PHRASES:
        _require(science, phrase, "docs/barcode-science-readiness.md", result)
    for phrase in REQUIRED_DEMO_PHRASES:
        _require(demo, phrase, "examples/barcode_autopsy/README.md", result)
    _reject_hype("\n".join([readme, science, demo]), "barcode public surface", result)
    if result.ok:
        result.passed.append("one-command barcode autopsy is discoverable")
        result.passed.append("barcode public surface stays within evidence boundaries")
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
        print("BARCODE PUBLIC SURFACE: PASS")
        return 0
    print("BARCODE PUBLIC SURFACE: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
