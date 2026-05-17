#!/usr/bin/env python3

import argparse
from pathlib import Path


def check(root: Path) -> list[str]:
    path = root / "docs" / "native-comparator-scope.md"
    if not path.is_file():
        return [f"missing native comparator scope document: {path.relative_to(root).as_posix()}"]
    if not path.read_text(encoding="utf-8").strip():
        return [f"empty native comparator scope document: {path.relative_to(root).as_posix()}"]
    return []


def main() -> int:
    parser = argparse.ArgumentParser(description="Check native comparator scope document presence.")
    parser.add_argument("--root", default=".", help="repository root")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    failures = check(root)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        print("NATIVE COMPARATOR SCOPE: FAIL")
        return 1
    print("PASS: native comparator scope document present")
    print("NATIVE COMPARATOR SCOPE: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
