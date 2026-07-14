#!/usr/bin/env python3
"""Fail closed when AssayCode platform surfaces drift or overclaim."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def require(path: str, needles: list[str]) -> list[str]:
    file_path = ROOT / path
    if not file_path.is_file():
        return [f"missing required AssayCode artifact: {path}"]
    text = file_path.read_text(encoding="utf-8")
    return [f"{path}: missing required text: {needle}" for needle in needles if needle not in text]


def main() -> int:
    errors: list[str] = []
    errors += require(
        "pyproject.toml",
        [
            'assaycode = "assaycode.cli:main"',
            'include = ["assaycode*", "dotmatch*", "quickdna*"]',
        ],
    )
    errors += require(
        "python/assaycode/cli.py",
        ["command_compile", "command_inspect", "command_watch"],
    )
    errors += require(
        "python/dotmatch/assayscript.py",
        ["class CompiledAssay", "source_sha256", "library_sha256", "safety_status"],
    )
    errors += require(
        "python/dotmatch/calibration.py",
        [
            "class ErrorModel",
            "def decode_joint",
            "def calibration_metrics",
            "def threshold_for_fdr",
            "This module is experimental",
        ],
    )
    errors += require(
        "python/dotmatch/assaywatch.py",
        ["class SequentialMonitor", "assignment_rate_interval95", "insufficient_data"],
    )
    for test in [
        "python/tests/test_assaycode_brand.py",
        "python/tests/test_assayscript.py",
        "python/tests/test_calibration.py",
        "python/tests/test_assaywatch.py",
    ]:
        errors += require(test, ["test_"])
    errors += require(
        "docs/assaycode.md",
        [
            "Compatibility Contract",
            "AssayScript v2 Compilation",
            "Experimental Calibration",
            "not yet claim",
        ],
    )
    errors += require(
        "docs/scientific-claims.md",
        [
            "AssayScript v2 compilation is experimental",
            "deterministic DotMatch assignment remains authoritative",
        ],
    )
    errors += require(
        "paper/paper.md",
        [
            "# AssayCode and AssayScript",
            "# Experimental uncertainty and run monitoring",
            "does not yet claim",
        ],
    )
    errors += require(
        "packaging/bioconda/meta.yaml",
        ["import assaycode, dotmatch", "assaycode --version"],
    )
    if errors:
        for error in errors:
            print(f"assaycode-readiness: {error}", file=sys.stderr)
        return 1
    print("assaycode-readiness: ready")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
