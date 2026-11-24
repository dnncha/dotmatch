#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path


def line_coverage_percent(coverage_json: Path) -> float:
    data = json.loads(coverage_json.read_text())
    return float(data["data"][0]["totals"]["lines"]["percent"])


if __name__ == "__main__":
    print(f"line coverage {line_coverage_percent(Path(sys.argv[1])):.2f}%")
