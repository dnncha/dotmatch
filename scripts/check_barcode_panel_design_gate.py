#!/usr/bin/env python3
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "barcode_panel_design.csv"


def main() -> int:
    subprocess.run([sys.executable, str(ROOT / "scripts" / "bench_barcode_panel_design.py")], cwd=ROOT, check=True)
    if not RAW.exists():
        raise SystemExit("barcode panel design benchmark did not write raw CSV")
    with RAW.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit("barcode panel design benchmark CSV is empty")
    by_workflow = {row["workflow"]: row for row in rows}
    required = {"small_symbolic_panel_smoke", "exact_k2_ambiguous_fixture", "contextual_flank_fixture"}
    missing = required - set(by_workflow)
    if missing:
        raise SystemExit(f"barcode panel design benchmark missing rows: {', '.join(sorted(missing))}")

    design = by_workflow["small_symbolic_panel_smoke"]
    if design["status"] not in {"pass", "warn"}:
        raise SystemExit(f"barcode panel design gate failed: status={design['status']}")
    if int(design["ambiguous_error_spheres"]) != 0:
        raise SystemExit("barcode panel design gate failed: ambiguous error spheres present")
    if int(design["silent_assignment_risk"]) != 0:
        raise SystemExit("barcode panel design gate failed: silent assignment risk present")

    k2 = by_workflow["exact_k2_ambiguous_fixture"]
    if k2["status"] != "fail" or int(k2["ambiguous_error_spheres"]) <= 0:
        raise SystemExit("barcode panel design gate failed: k=2 ambiguous fixture was not rejected")

    context = by_workflow["contextual_flank_fixture"]
    if context["status"] not in {"warn", "fail"} or int(context["context_warnings"]) <= 0:
        raise SystemExit("barcode panel design gate failed: contextual flank fixture was not reported")

    print(f"barcode panel design gate passed: {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
