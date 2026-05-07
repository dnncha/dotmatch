#!/usr/bin/env python3
"""Fail when llvm-cov JSON line coverage is below a requested threshold."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coverage_json")
    parser.add_argument("--min-lines", type=float, required=True)
    args = parser.parse_args()

    data = json.loads(Path(args.coverage_json).read_text())
    percent = float(data["data"][0]["totals"]["lines"]["percent"])
    if percent < args.min_lines:
        raise SystemExit(f"line coverage {percent:.2f}% is below required {args.min_lines:.2f}%")
    print(f"line coverage {percent:.2f}% >= {args.min_lines:.2f}%")


if __name__ == "__main__":
    main()
