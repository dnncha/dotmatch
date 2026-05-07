#!/usr/bin/env python3
"""Fetch a pinned Edlib source tree for native benchmarks."""

from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EDLIB_DIR = ROOT / "build" / "edlib"
EDLIB_REPO = "https://github.com/Martinsos/edlib.git"
EDLIB_REF = "v1.2.7"


def main() -> None:
    EDLIB_DIR.parent.mkdir(parents=True, exist_ok=True)
    if not EDLIB_DIR.exists():
        subprocess.run(["git", "clone", "--depth", "1", "--branch", EDLIB_REF, EDLIB_REPO, str(EDLIB_DIR)], check=True)
    else:
        subprocess.run(["git", "fetch", "--depth", "1", "origin", EDLIB_REF], cwd=EDLIB_DIR, check=True)
        subprocess.run(["git", "checkout", EDLIB_REF], cwd=EDLIB_DIR, check=True)
    print(f"edlib_dir={EDLIB_DIR}")
    print(f"edlib_ref={EDLIB_REF}")


if __name__ == "__main__":
    main()
