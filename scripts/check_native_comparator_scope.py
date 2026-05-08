#!/usr/bin/env python3

import argparse
from pathlib import Path


REQUIRED_FRAGMENTS = [
    "Edlib exhaustive global edit-distance assignment",
    "EDLIB_MODE_NW",
    "EDLIB_TASK_DISTANCE",
    "SeqAn and Parasail comparisons are not claimed as completed release evidence yet",
    "equivalent global edit-distance or documented semi-global scoring semantics",
    "fixed threshold `k`",
    "native dependency name, version, build flags, and platform",
    "raw CSV rows under `benchmarks/raw/`",
    "zero assignment mismatches",
    "gate script that fails when only scaffold, smoke, or unmatched-scoring rows are present",
    "limited to Edlib exhaustive global edit-distance assignment scans",
]


def check(root: Path) -> list[str]:
    path = root / "docs" / "native-comparator-scope.md"
    if not path.is_file():
        return [f"missing native comparator scope document: {path.relative_to(root).as_posix()}"]

    text = path.read_text(encoding="utf-8")
    failures = [
        f"docs/native-comparator-scope.md missing required evidence boundary: {fragment}"
        for fragment in REQUIRED_FRAGMENTS
        if fragment not in text
    ]
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description="Check native comparator scope and SeqAn/Parasail evidence boundaries.")
    parser.add_argument("--root", default=".", help="repository root")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    failures = check(root)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        print("NATIVE COMPARATOR SCOPE: FAIL")
        return 1
    print("PASS: native comparator scope and SeqAn/Parasail evidence boundaries documented")
    print("NATIVE COMPARATOR SCOPE: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
