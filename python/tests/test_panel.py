from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LEGACY_ENV = {**os.environ, "DOTMATCH_PYTHON_NO_DELEGATE": "1", "PYTHONPATH": str(ROOT / "python")}


def _run_cli(args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "dotmatch.cli", *args],
        cwd=ROOT,
        env=LEGACY_ENV,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _rows(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def test_panel_check_fails_ambiguous_k1_hamming(tmp_path: Path) -> None:
    panel = tmp_path / "panel.tsv"
    panel.write_text("barcode_id\tsequence\nBC001\tAAAA\nBC002\tAATT\n", encoding="utf-8")
    out_dir = tmp_path / "check"

    rc = _run_cli(["panel", "check", str(panel), "--k", "1", "--metric", "hamming", "--out-dir", str(out_dir)])

    assert rc.returncode == 1, rc.stderr
    summary = json.loads((out_dir / "panel_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "fail"
    assert summary["assignment_metric"] == "hamming"
    assert summary["configured_assignment_k"] == 1
    assert summary["ambiguous_error_spheres"] > 0
    assert summary["unsafe_k1_variants"] > 0
    assert summary["safe_for_k1_hamming"] is False
    assert "dotmatch demux" in summary["certified_dotmatch_command"]
    assert _rows(out_dir / "ambiguous_error_spheres.tsv")
    assert (out_dir / "panel_report.html").read_text(encoding="utf-8").count("Do not use") >= 1


def test_panel_check_detects_silent_assignment_and_reverse_complement_risk(tmp_path: Path) -> None:
    panel = tmp_path / "panel.tsv"
    panel.write_text("barcode_id\tsequence\nBC001\tACGT\nBC002\tACGA\nBC003\tTCGT\n", encoding="utf-8")
    out_dir = tmp_path / "check"

    rc = _run_cli(
        [
            "panel",
            "check",
            str(panel),
            "--k",
            "1",
            "--metric",
            "hamming",
            "--reverse-complement-mode",
            "warn",
            "--out-dir",
            str(out_dir),
        ]
    )

    assert rc.returncode == 1, rc.stderr
    summary = json.loads((out_dir / "panel_summary.json").read_text(encoding="utf-8"))
    assert summary["silent_assignment_risk"] > 0
    assert summary["reverse_complement_warnings"] > 0
    safety = _rows(out_dir / "target_safety.tsv")
    assert any(row["status"] == "fail" for row in safety)


def test_panel_check_fails_ambiguous_k2_hamming_exact_sphere(tmp_path: Path) -> None:
    panel = tmp_path / "panel.tsv"
    panel.write_text("barcode_id\tsequence\nBC001\tAAAA\nBC002\tTTTT\n", encoding="utf-8")
    out_dir = tmp_path / "check"

    rc = _run_cli(["panel", "check", str(panel), "--k", "2", "--metric", "hamming", "--out-dir", str(out_dir)])

    assert rc.returncode == 1, rc.stderr
    summary = json.loads((out_dir / "panel_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "fail"
    assert summary["configured_assignment_k"] == 2
    assert summary["safe_for_k2_hamming"] is False
    rows = _rows(out_dir / "ambiguous_error_spheres.tsv")
    assert any(row["variant"] == "AATT" for row in rows)


def test_panel_check_refuses_uncertified_radius_above_two(tmp_path: Path) -> None:
    panel = tmp_path / "panel.tsv"
    panel.write_text("barcode_id\tsequence\nBC001\tAACCGGTT\nBC002\tTTGGCCAA\n", encoding="utf-8")
    out_dir = tmp_path / "check"

    rc = _run_cli(["panel", "check", str(panel), "--k", "3", "--metric", "hamming", "--out-dir", str(out_dir)])

    assert rc.returncode == 2
    assert "exact safety certificate supports k <= 2" in rc.stderr


def test_panel_check_writes_contextual_certificate_outputs(tmp_path: Path) -> None:
    panel = tmp_path / "panel.tsv"
    panel.write_text("barcode_id\tsequence\nBC001\tATGC\n", encoding="utf-8")
    out_dir = tmp_path / "check"

    rc = _run_cli(
        [
            "panel",
            "check",
            str(panel),
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
            str(out_dir),
        ]
    )

    assert rc.returncode == 0, rc.stderr
    context = _rows(out_dir / "context_risk.tsv")
    flanked = _rows(out_dir / "flanked_sequences.tsv")
    assert context[0]["status"] == "warn"
    assert "left_flank_homopolymer" in context[0]["risks"]
    assert flanked[0]["flanked_sequence"] == "AAAATGCTTT"
    assert "dotmatch demux" in context[0]["certified_command"]


def test_panel_design_is_seed_reproducible_and_certified(tmp_path: Path) -> None:
    out_a = tmp_path / "a"
    out_b = tmp_path / "b"
    args = [
        "panel",
        "design",
        "--n",
        "8",
        "--length",
        "8",
        "--seed",
        "42",
        "--candidate-pool-size",
        "5000",
        "--restarts",
        "3",
        "--min-hamming-distance",
        "4",
        "--min-levenshtein-distance",
        "4",
    ]

    assert _run_cli([*args, "--out-dir", str(out_a)]).returncode == 0
    assert _run_cli([*args, "--out-dir", str(out_b)]).returncode == 0

    assert (out_a / "barcodes.tsv").read_text(encoding="utf-8") == (out_b / "barcodes.tsv").read_text(encoding="utf-8")
    rows = _rows(out_a / "barcodes.tsv")
    assert len(rows) == 8
    assert {"barcode_id", "sequence", "gc", "min_hamming_neighbor", "certified_command"} <= set(rows[0])
    report = json.loads((out_a / "design_report.json").read_text(encoding="utf-8"))
    assert report["engine"] == "greedy"
    assert report["seed"] == 42
    assert report["certificate"]["status"] in {"pass", "warn"}
    assert (out_a / "panel_check" / "panel_summary.json").exists()
    assert (out_a / "sample_sheet_templates" / "SampleSheet.csv").exists()
    assert "machine-checkable safety certificate" in (out_a / "README_FOR_LAB.md").read_text(encoding="utf-8")


def test_panel_design_accepts_yaml_spec_and_graph_engine(tmp_path: Path) -> None:
    spec = tmp_path / "barcode_design.yml"
    spec.write_text(
        """
panel:
  name: "yaml_panel"
  count: 6
  length: 8
  seed: 11
assignment:
  metric: "hamming"
  k: 1
constraints:
  min_hamming_distance: 4
  min_levenshtein_distance: 4
  gc_min: 0.25
  gc_max: 0.75
  max_homopolymer: 3
cycle_balance:
  enabled: true
plate_layout:
  enabled: true
  format: 96
""".lstrip(),
        encoding="utf-8",
    )
    out_dir = tmp_path / "yaml"

    rc = _run_cli(["panel", "design", "--spec", str(spec), "--engine", "graph", "--candidate-pool-size", "2000", "--out-dir", str(out_dir)])

    assert rc.returncode == 0, rc.stderr
    report = json.loads((out_dir / "design_report.json").read_text(encoding="utf-8"))
    assert report["panel_name"] == "yaml_panel"
    assert report["engine"] == "graph"
    assert len(_rows(out_dir / "barcodes.tsv")) == 6


def test_panel_simulate_layout_export_optimize_compare_and_dual_design(tmp_path: Path) -> None:
    candidates = tmp_path / "candidates.tsv"
    candidates.write_text(
        "barcode_id\tsequence\n"
        "A\tTACGACAC\n"
        "B\tCCATTGGT\n"
        "C\tAGTCGTCG\n"
        "D\tGCTCACTA\n"
        "E\tTTAACCGG\n"
        "F\tCGAGAATG\n",
        encoding="utf-8",
    )
    opt = tmp_path / "opt"
    rc = _run_cli(
        [
            "panel",
            "optimize",
            str(candidates),
            "--n",
            "3",
            "--seed",
            "7",
            "--min-hamming-distance",
            "4",
            "--min-levenshtein-distance",
            "4",
            "--out-dir",
            str(opt),
        ]
    )
    assert rc.returncode == 0, rc.stderr
    panel = opt / "barcodes.tsv"
    assert len(_rows(panel)) == 3

    sim = tmp_path / "sim"
    rc = _run_cli(
        [
            "panel",
            "simulate",
            str(panel),
            "--reads",
            "200",
            "--seed",
            "7",
            "--substitution-rate",
            "0.02",
            "--insertion-rate",
            "0.001",
            "--deletion-rate",
            "0.001",
            "--out-dir",
            str(sim),
        ]
    )
    assert rc.returncode == 0, rc.stderr
    sim_summary = json.loads((sim / "simulation_summary.json").read_text(encoding="utf-8"))
    assert sim_summary["total_reads"] == 200
    assert {"unique_rate", "ambiguous_rate", "none_rate", "invalid_rate", "false_assignment_rate"} <= set(sim_summary)
    assert (sim / "per_barcode_confusion.tsv").exists()

    plate = tmp_path / "plate.tsv"
    rc = _run_cli(["panel", "layout", str(panel), "--plate", "96", "--out", str(plate)])
    assert rc.returncode == 0, rc.stderr
    assert _rows(plate)[0]["well"] == "A1"
    assert (tmp_path / "plate.svg").exists()
    assert (tmp_path / "neighbor_distance.tsv").exists()
    assert (tmp_path / "lab_picklist.csv").exists()

    export = tmp_path / "export"
    rc = _run_cli(["panel", "export", str(panel), "--format", "illumina-samplesheet", "--out-dir", str(export)])
    assert rc.returncode == 0, rc.stderr
    assert "[Data]" in (export / "SampleSheet.csv").read_text(encoding="utf-8")
    assert "dotmatch demux" in (export / "README_FOR_LAB.md").read_text(encoding="utf-8")

    compare = tmp_path / "compare"
    rc = _run_cli(["panel", "compare", str(candidates), str(panel), "--out-dir", str(compare)])
    assert rc.returncode == 0, rc.stderr
    compare_summary = json.loads((compare / "panel_compare.json").read_text(encoding="utf-8"))
    assert compare_summary["old_count"] == 6
    assert compare_summary["new_count"] == 3

    dual = tmp_path / "dual"
    rc = _run_cli(
        [
            "panel",
            "design-dual",
            "--samples",
            "6",
            "--i7-count",
            "6",
            "--i5-count",
            "6",
            "--i7-length",
            "8",
            "--i5-length",
            "8",
            "--unique-dual",
            "--min-i7-distance",
            "4",
            "--min-i5-distance",
            "4",
            "--min-pair-distance",
            "8",
            "--seed",
            "9",
            "--out-dir",
            str(dual),
        ]
    )
    assert rc.returncode == 0, rc.stderr
    dual_rows = _rows(dual / "dual_barcodes.tsv")
    assert len(dual_rows) == 6
    assert dual_rows[0]["index_hop_detectable"] == "true"
    assert (dual / "sample_sheet_templates" / "SampleSheet.csv").exists()


def test_panel_docs_scope_and_gate_exist() -> None:
    docs = (ROOT / "docs/barcode-panel-design.md").read_text(encoding="utf-8")
    assert "known-target assignment" in docs
    assert "not general genome alignment" in docs
    assert "not UMI entropy generation" in docs
    assert "not basecalling" in docs
    assert "DotMatch does not merely design barcodes" in docs
    assert (ROOT / "scripts/check_barcode_panel_design_gate.py").exists()
    assert (ROOT / "scripts/bench_barcode_panel_design.py").exists()
