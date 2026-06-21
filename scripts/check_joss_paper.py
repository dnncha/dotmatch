#!/usr/bin/env python3
"""Validate the local JOSS paper draft enough to catch submission drift."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


REQUIRED_HEADINGS = [
    "Summary",
    "Statement of need",
    "State of the field",
    "Software design",
    "Research impact statement",
    "AI usage disclosure",
    "Acknowledgements",
    "References",
]
MIN_WORDS = 750
MAX_WORDS = 1750


class Result:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _word_count(markdown: str) -> int:
    body = re.sub(r"^---\n.*?\n---\n", "", markdown, flags=re.S)
    body = re.sub(r"\[@[A-Za-z0-9_:-]+\]", "", body)
    body = re.sub(r"https?://\S+", "", body)
    return len(re.findall(r"[A-Za-z0-9][A-Za-z0-9'/-]*", body))


def _metadata(markdown: str) -> str:
    match = re.match(r"^---\n(.*?)\n---\n", markdown, flags=re.S)
    return match.group(1) if match else ""


def _bib_keys(bibtex: str) -> set[str]:
    return set(re.findall(r"@\w+\{([^,\s]+)", bibtex))


def _used_keys(markdown: str) -> set[str]:
    return set(re.findall(r"@([A-Za-z0-9_:-]+)", markdown))


def audit(root: Path) -> Result:
    root = root.resolve()
    result = Result()
    paper_path = root / "paper" / "paper.md"
    bib_path = root / "paper" / "paper.bib"
    if not paper_path.exists():
        result.failures.append("paper/paper.md is missing")
        return result
    if not bib_path.exists():
        result.failures.append("paper/paper.bib is missing")
        return result

    paper = _read(paper_path)
    bib = _read(bib_path)
    metadata = _metadata(paper)
    if not metadata:
        result.failures.append("paper/paper.md must start with JOSS YAML metadata")
    for required in ["title:", "tags:", "authors:", "affiliations:", "date:", "bibliography: paper.bib"]:
        if required not in metadata:
            result.failures.append(f"paper metadata missing {required}")
    if "archive_doi: 10.5281/zenodo.20541628" not in metadata:
        result.failures.append("paper metadata must include the DotMatch Zenodo concept DOI")

    headings = set(re.findall(r"^# (.+)$", paper, flags=re.M))
    for heading in REQUIRED_HEADINGS:
        if heading not in headings:
            result.failures.append(f"paper missing required JOSS section: {heading}")

    words = _word_count(paper)
    if not MIN_WORDS <= words <= MAX_WORDS:
        result.failures.append(f"paper word count {words} outside JOSS range {MIN_WORDS}-{MAX_WORDS}")

    if re.search(r"\bTODO\b", paper, flags=re.I):
        result.failures.append("paper must not contain TODO markers")
    if "AI usage disclosure" in headings and "OpenAI Codex" not in paper:
        result.failures.append("paper AI usage disclosure must describe Codex usage")

    bib_keys = _bib_keys(bib)
    used_keys = _used_keys(paper)
    missing = sorted(used_keys - bib_keys)
    unused = sorted(bib_keys - used_keys)
    if missing:
        result.failures.append(f"paper cites missing bibliography keys: {', '.join(missing)}")
    if unused:
        result.failures.append(f"paper.bib has unused keys: {', '.join(unused)}")

    if result.ok:
        result.passed.append(f"JOSS paper draft valid ({words} words)")
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
        print("JOSS PAPER: PASS")
        return 0
    print("JOSS PAPER: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
