import csv
import importlib.util
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "scripts" / "check_bcl_comparison_gate.py"


def _load_gate():
    spec = importlib.util.spec_from_file_location("check_bcl_comparison_gate", GATE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _passing_row(tool: str, workflow: str, fmt: str, speed: str, repeat: int) -> dict[str, str]:
    return {
        "tool": tool,
        "workflow": workflow,
        "format": fmt,
        "clusters_per_sec": speed,
        "validation_mismatches": "0",
        "validation_exit_code": "0",
        "exit_code": "0",
        "repeat": str(repeat),
    }


def _otherwise_passing_rows(workflow: str, repeats: int = 5) -> list[dict[str, str]]:
    rows = []
    for repeat in range(1, repeats + 1):
        rows.append(_passing_row("dotmatch_bcl_demux", workflow, "classic_bcl", "1000", repeat))
        rows.append(_passing_row("dotmatch_bcl_demux", workflow, "cbcl", "1000", repeat))
        rows.append(_passing_row("bcl2fastq", workflow, "classic_bcl", "50", repeat))
    return rows


def test_bcl_gate_rejects_tiny_demo_rows_even_when_other_conditions_pass(tmp_path, monkeypatch):
    gate = _load_gate()
    csv_path = tmp_path / "bcl.csv"
    _write_rows(csv_path, _otherwise_passing_rows("public_10x_tiny_bcl"))
    monkeypatch.setattr(sys, "argv", ["check_bcl_comparison_gate.py", "--csv", str(csv_path)])

    with pytest.raises(SystemExit) as exc:
        gate.main()

    assert "tiny" in str(exc.value)


def test_bcl_gate_requires_repeated_dotmatch_and_validated_competitor_rows(tmp_path, monkeypatch):
    gate = _load_gate()
    csv_path = tmp_path / "bcl.csv"
    _write_rows(csv_path, _otherwise_passing_rows("real_cbcl_and_classic_run", repeats=1))
    monkeypatch.setattr(sys, "argv", ["check_bcl_comparison_gate.py", "--csv", str(csv_path)])

    with pytest.raises(SystemExit) as exc:
        gate.main()

    assert "repeated" in str(exc.value)


def test_bcl_gate_requires_distinct_repeats_not_duplicated_rows(tmp_path, monkeypatch):
    gate = _load_gate()
    csv_path = tmp_path / "bcl.csv"
    rows = []
    for _ in range(5):
        rows.append(_passing_row("dotmatch_bcl_demux", "real_cbcl_and_classic_run", "classic_bcl", "1000", 1))
        rows.append(_passing_row("dotmatch_bcl_demux", "real_cbcl_and_classic_run", "cbcl", "1000", 1))
        rows.append(_passing_row("bcl2fastq", "real_cbcl_and_classic_run", "classic_bcl", "50", 1))
    _write_rows(csv_path, rows)
    monkeypatch.setattr(sys, "argv", ["check_bcl_comparison_gate.py", "--csv", str(csv_path)])

    with pytest.raises(SystemExit) as exc:
        gate.main()

    assert "distinct successful repeats" in str(exc.value)


def test_bcl_gate_requires_dotmatch_cbcl_not_competitor_only_cbcl(tmp_path, monkeypatch):
    gate = _load_gate()
    csv_path = tmp_path / "bcl.csv"
    rows = []
    for repeat in range(1, 6):
        rows.append(_passing_row("dotmatch_bcl_demux", "real_cbcl_and_classic_run", "classic_bcl", "1000", repeat))
        rows.append(_passing_row("bcl2fastq", "real_cbcl_and_classic_run", "cbcl", "50", repeat))
    _write_rows(csv_path, rows)
    monkeypatch.setattr(sys, "argv", ["check_bcl_comparison_gate.py", "--csv", str(csv_path)])

    with pytest.raises(SystemExit) as exc:
        gate.main()

    assert "DotMatch CBCL" in str(exc.value)
