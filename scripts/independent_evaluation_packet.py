"""Compatibility wrapper for the packaged independent-evaluation tool."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "python/dotmatch/evaluation_packet.py"
SPEC = importlib.util.spec_from_file_location("dotmatch_evaluation_packet", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(
        "could not load the packaged independent-evaluation implementation"
    )
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)
main = MODULE.main


if __name__ == "__main__":
    raise SystemExit(main())
