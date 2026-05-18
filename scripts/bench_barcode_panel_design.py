#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_command(args: list[str], out_dir: Path, *, expected_returncodes: set[int] | None = None) -> tuple[float, subprocess.CompletedProcess[str]]:
    expected_returncodes = expected_returncodes or {0}
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "python")
    env["DOTMATCH_PYTHON_NO_DELEGATE"] = "1"
    start = time.perf_counter()
    completed = subprocess.run(
        [sys.executable, "-m", "dotmatch.cli", *args],
        cwd=ROOT,
        env=env,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    elapsed = time.perf_counter() - start
    if completed.returncode not in expected_returncodes:
        raise SystemExit(f"command failed in {out_dir}: {completed.stderr}")
    return elapsed, completed


def main() -> int:
    raw_dir = ROOT / "benchmarks" / "raw"
    report_dir = ROOT / "docs" / "benchmarks" / "barcode_panel_design"
    raw_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="dotmatch-panel-bench-") as tmp:
        work = Path(tmp)
        design_dir = work / "design"
        elapsed_design, _ = run_command(
            [
                "panel",
                "design",
                "--n",
                "8",
                "--length",
                "8",
                "--seed",
                "2026",
                "--candidate-pool-size",
                "2000",
                "--restarts",
                "2",
                "--min-hamming-distance",
                "4",
                "--min-levenshtein-distance",
                "4",
                "--out-dir",
                str(design_dir),
            ],
            design_dir,
        )
        check_dir = work / "check"
        elapsed_check, _ = run_command(
            ["panel", "check", str(design_dir / "barcodes.tsv"), "--k", "1", "--metric", "hamming", "--out-dir", str(check_dir)],
            check_dir,
        )
        summary = json.loads((check_dir / "panel_summary.json").read_text(encoding="utf-8"))
        k2_panel = work / "k2_unsafe.tsv"
        k2_panel.write_text("barcode_id\tsequence\nBC001\tAAAA\nBC002\tTTTT\n", encoding="utf-8")
        k2_dir = work / "k2_unsafe_check"
        elapsed_k2, _ = run_command(
            ["panel", "check", str(k2_panel), "--k", "2", "--metric", "hamming", "--out-dir", str(k2_dir)],
            k2_dir,
            expected_returncodes={1},
        )
        k2_summary = json.loads((k2_dir / "panel_summary.json").read_text(encoding="utf-8"))

        context_panel = work / "context.tsv"
        context_panel.write_text("barcode_id\tsequence\nBC001\tATGC\n", encoding="utf-8")
        context_dir = work / "context_check"
        elapsed_context, _ = run_command(
            [
                "panel",
                "check",
                str(context_panel),
                "--k",
                "1",
                "--metric",
                "levenshtein",
                "--left-flank",
                "AAA",
                "--right-flank",
                "TTT",
                "--context-window",
                "4",
                "--out-dir",
                str(context_dir),
            ],
            context_dir,
        )
        context_summary = json.loads((context_dir / "panel_summary.json").read_text(encoding="utf-8"))

    raw_path = raw_dir / "barcode_panel_design.csv"
    with raw_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "workflow",
                "status",
                "barcodes",
                "length",
                "design_seconds",
                "check_seconds",
                "minimum_hamming_distance",
                "ambiguous_error_spheres",
                "silent_assignment_risk",
                "context_warnings",
                "expected_outcome",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "workflow": "small_symbolic_panel_smoke",
                "status": summary["status"],
                "barcodes": summary["n_barcodes"],
                "length": ",".join(str(item) for item in summary["lengths"]),
                "design_seconds": f"{elapsed_design:.4f}",
                "check_seconds": f"{elapsed_check:.4f}",
                "minimum_hamming_distance": summary["minimum_hamming_distance"],
                "ambiguous_error_spheres": summary["ambiguous_error_spheres"],
                "silent_assignment_risk": summary["silent_assignment_risk"],
                "context_warnings": 0,
                "expected_outcome": "certified_safe_or_warn",
            }
        )
        writer.writerow(
            {
                "workflow": "exact_k2_ambiguous_fixture",
                "status": k2_summary["status"],
                "barcodes": k2_summary["n_barcodes"],
                "length": ",".join(str(item) for item in k2_summary["lengths"]),
                "design_seconds": "",
                "check_seconds": f"{elapsed_k2:.4f}",
                "minimum_hamming_distance": k2_summary["minimum_hamming_distance"],
                "ambiguous_error_spheres": k2_summary["ambiguous_error_spheres"],
                "silent_assignment_risk": k2_summary["silent_assignment_risk"],
                "context_warnings": 0,
                "expected_outcome": "must_fail",
            }
        )
        writer.writerow(
            {
                "workflow": "contextual_flank_fixture",
                "status": context_summary["status"],
                "barcodes": context_summary["n_barcodes"],
                "length": ",".join(str(item) for item in context_summary["lengths"]),
                "design_seconds": "",
                "check_seconds": f"{elapsed_context:.4f}",
                "minimum_hamming_distance": context_summary["minimum_hamming_distance"],
                "ambiguous_error_spheres": context_summary["ambiguous_error_spheres"],
                "silent_assignment_risk": context_summary["silent_assignment_risk"],
                "context_warnings": 1,
                "expected_outcome": "must_warn",
            }
        )
    (report_dir / "README.md").write_text(
        "# Barcode Panel Design Benchmark\n\n"
        "Smoke benchmark for DotMatch barcode panel design and machine-checkable safety certificates.\n\n"
        f"- Raw data: `{raw_path.relative_to(ROOT)}`\n"
        f"- Panel status: `{summary['status']}`\n"
        f"- Minimum Hamming distance: `{summary['minimum_hamming_distance']}`\n"
        f"- k=2 unsafe fixture status: `{k2_summary['status']}`\n"
        f"- Context fixture status: `{context_summary['status']}`\n",
        encoding="utf-8",
    )
    print(raw_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
