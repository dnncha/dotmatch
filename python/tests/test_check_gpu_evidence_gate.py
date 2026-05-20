import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "scripts" / "check_gpu_evidence_gate.py"


def _load_gate():
    spec = importlib.util.spec_from_file_location("check_gpu_evidence_gate", GATE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _row(tool: str, **overrides) -> dict[str, str]:
    row = {
        "tool": tool,
        "backend": "metal" if tool == "dotmatch_gpu_metal" else "cpu",
        "status": "ok",
        "workload": "synthetic_hamming",
        "n_reads": "1000",
        "n_targets": "96",
        "len": "20",
        "k": "1",
        "error_rate": "0.01",
        "checksum": "123",
        "mismatches": "0",
        "device": "Apple test GPU",
    }
    row.update(overrides)
    return row


def test_gpu_gate_accepts_zero_mismatch_gpu_with_cpu_baseline():
    gate = _load_gate()
    failures = []

    gate.row_gate([_row("dotmatch_cpu_index"), _row("dotmatch_gpu_metal")], failures)

    assert failures == []


def test_gpu_gate_rejects_checksum_mismatch():
    gate = _load_gate()
    failures = []

    gate.row_gate([
        _row("dotmatch_cpu_index", checksum="123"),
        _row("dotmatch_gpu_metal", checksum="999"),
    ], failures)

    assert any("checksum differs" in failure for failure in failures)


def test_gpu_gate_accepts_explicit_unavailable_row():
    gate = _load_gate()
    failures = []

    gate.row_gate([
        _row("dotmatch_gpu_metal", status="unavailable", n_reads="0", n_targets="0", checksum="0", device="Darwin"),
    ], failures)

    assert failures == []


def test_gpu_gate_accepts_public_crispr_zero_mismatch_row():
    gate = _load_gate()
    failures = []
    cpu = {
        "tool": "dotmatch_cpu_index",
        "status": "ok",
        "packable_reads": "1000",
        "n_targets": "96",
        "target_length": "19",
        "k": "1",
        "checksum": "123",
    }
    gpu = {
        **cpu,
        "tool": "dotmatch_gpu_metal",
        "mismatches": "0",
        "count_delta": "0",
        "device": "Apple test GPU",
    }

    gate.real_row_gate([cpu, gpu], failures)

    assert failures == []


def test_gpu_gate_rejects_public_crispr_count_delta():
    gate = _load_gate()
    failures = []
    cpu = {
        "tool": "dotmatch_cpu_index",
        "status": "ok",
        "packable_reads": "1000",
        "n_targets": "96",
        "target_length": "19",
        "k": "1",
        "checksum": "123",
    }
    gpu = {
        **cpu,
        "tool": "dotmatch_gpu_metal",
        "mismatches": "0",
        "count_delta": "3",
        "device": "Apple test GPU",
    }

    gate.real_row_gate([cpu, gpu], failures)

    assert any("count delta" in failure for failure in failures)


def test_gpu_report_gate_requires_experimental_boundary(tmp_path):
    gate = _load_gate()
    report = tmp_path / "README.md"
    report.write_text("# Experimental GPU Acceleration Benchmark\n\nnot a production speed claim\n", encoding="utf-8")
    failures = []

    gate.report_gate(report, failures)

    assert failures == []
