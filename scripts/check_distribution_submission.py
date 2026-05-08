#!/usr/bin/env python3
"""Validate the public distribution submission handoff dossier."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


DOSSIER = Path("docs") / "distribution-submission.md"
REQUIRED_CHECKS = [
    "make release-ready",
    "make python-package-test",
    "make bioconda-recipe-ready",
    "make distribution-record-ready",
    "make distribution-submission-ready",
    "make distribution-channels",
]
REQUIRED_TARGET_FRAGMENTS = {
    "PyPI": [
        "PyPI trusted publishing",
        "source distribution",
        "repaired manylinux/musllinux wheels",
        "dotmatch-linux-repaired-wheels",
        "scripts/check_python_wheel.py --sdist-only --out-dir dist",
    ],
    "GHCR": ["GHCR", "ghcr.io/dnncha/dotmatch", "OCI labels", "dotmatch --version", "dotmatch dist ACGT AGGT"],
    "Bioconda": ["Bioconda", "packaging/bioconda/meta.yaml", "REPLACE_WITH_RELEASE_TARBALL_SHA256", "bioconda-recipes"],
    "BioContainers": ["BioContainers", "quay.io/biocontainers/dotmatch"],
    "Zenodo": ["Zenodo", "CITATION.cff", "DOI"],
}
REQUIRED_RECORD_FRAGMENTS = [
    "docs/distribution-release.json",
    "public_url",
    "evidence_url",
    "verified_date",
    "status to `released`",
    "stable public HTTPS URLs",
    "`make distribution-channels` passes",
]
PUBLIC_CLAIM_PATTERNS = [
    r"\bis now available on\b",
    r"\bhas been published to\b",
    r"\bis published on\b",
    r"\bavailable through PyPI and Bioconda\b",
]


class AuditResult:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def _read(path: Path, result: AuditResult) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except Exception as exc:
        result.failures.append(f"{path.as_posix()} could not be read: {exc}")
        return ""


def _require(text: str, fragment: str, message: str, result: AuditResult) -> None:
    if fragment not in text:
        result.failures.append(message)


def audit(root: Path) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    text = _read(root / DOSSIER, result)
    if not text:
        return result

    _require(text, "Distribution Submission Dossier", "distribution submission dossier missing title", result)
    _require(
        text,
        "not a public distribution claim",
        "distribution submission dossier must state that it is not a public distribution claim",
        result,
    )
    for pattern in PUBLIC_CLAIM_PATTERNS:
        if re.search(pattern, text, flags=re.I):
            result.failures.append("distribution submission dossier must avoid public distribution claims")
            break

    for command in REQUIRED_CHECKS:
        _require(text, command, f"distribution submission dossier must include {command}", result)

    for target, fragments in REQUIRED_TARGET_FRAGMENTS.items():
        for fragment in fragments:
            _require(text, fragment, f"distribution submission dossier must include {target}: {fragment}", result)

    for fragment in REQUIRED_RECORD_FRAGMENTS:
        _require(text, fragment, f"distribution submission dossier must include {fragment}", result)

    if result.ok:
        result.passed.append("distribution submission dossier present")
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
        print("DISTRIBUTION SUBMISSION: PASS")
        return 0
    print("DISTRIBUTION SUBMISSION: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
