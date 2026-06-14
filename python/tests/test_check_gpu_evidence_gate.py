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
        "workload": "public_crispr",
        "total_reads": "1200",
        "packable_reads": "1000",
        "n_targets": "96",
        "target_start": "23",
        "target_length": "19",
        "k": "1",
        "skipped_targets": "0",
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
        "workload": "public_crispr",
        "total_reads": "1200",
        "packable_reads": "1000",
        "n_targets": "96",
        "target_start": "23",
        "target_length": "19",
        "k": "1",
        "skipped_targets": "0",
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


def test_gpu_gate_rejects_public_crispr_mismatched_case_or_skipped_targets():
    gate = _load_gate()
    failures = []
    cpu = {
        "tool": "dotmatch_cpu_index",
        "status": "ok",
        "workload": "public_crispr",
        "total_reads": "1200",
        "packable_reads": "1000",
        "n_targets": "96",
        "target_start": "23",
        "target_length": "19",
        "k": "1",
        "skipped_targets": "0",
        "checksum": "123",
    }
    gpu = {
        **cpu,
        "tool": "dotmatch_gpu_metal",
        "target_start": "24",
        "mismatches": "0",
        "count_delta": "0",
        "device": "Apple test GPU",
    }

    gate.real_row_gate([cpu, gpu], failures)

    assert any("missing public CRISPR CPU baseline" in failure for failure in failures)

    failures = []
    gpu = {**cpu, "tool": "dotmatch_gpu_metal", "skipped_targets": "2", "mismatches": "0", "count_delta": "0", "device": "Apple test GPU"}
    gate.real_row_gate([cpu, gpu], failures)

    assert any("skipped targets" in failure for failure in failures)


def _production_cpu_row(**overrides) -> dict[str, str]:
    row = {
        "dataset": "mageck_yusa",
        "backend": "cpu",
        "metal_validate": "0",
        "status": "ok",
        "records_per_sample": "10000",
        "total_reads": "20000",
        "n_targets": "87437",
        "offset_mode": "fixed",
        "auto_offset": "0",
        "wall_seconds": "0.10",
        "reads_per_sec": "200000.0",
        "backend_effective": "cpu",
        "count_engine": "hamming_lookup_direct_single_offset",
        "metal_validation": "None",
        "count_match_cpu": "",
        "device": "test device",
    }
    row.update(overrides)
    return row


def _production_metal_row(**overrides) -> dict[str, str]:
    row = {
        **_production_cpu_row(),
        "backend": "gpu-metal-experimental",
        "wall_seconds": "0.08",
        "reads_per_sec": "250000.0",
        "backend_effective": "gpu-metal-experimental",
        "count_engine": "hamming_metal_seed_index",
        "count_match_cpu": "0",
    }
    row.update(overrides)
    return row


def _production_validation_failure(**overrides) -> dict[str, str]:
    row = {
        **_production_metal_row(),
        "metal_validate": "1",
        "status": "validation_failed",
        "total_reads": "0",
        "n_targets": "",
        "reads_per_sec": "",
        "backend_effective": "",
        "count_engine": "",
        "count_match_cpu": "",
        "device": "test device",
        "notes": "Metal validation failed against CPU authority checksum",
    }
    row.update(overrides)
    return row


def _production_sanson_ineligible(**overrides) -> dict[str, str]:
    row = {
        **_production_metal_row(),
        "dataset": "sanson_brunello",
        "status": "ineligible",
        "records_per_sample": "10000",
        "total_reads": "",
        "n_targets": "",
        "offset_mode": "multi",
        "auto_offset": "20",
        "wall_seconds": "",
        "reads_per_sec": "",
        "backend_effective": "",
        "count_engine": "",
        "count_match_cpu": "",
        "notes": "offset-mode multi blocks the production Metal path",
    }
    row.update(overrides)
    return row


def test_gpu_gate_accepts_production_crispr_metal_with_cpu_authority_and_boundaries():
    gate = _load_gate()
    failures = []

    gate.production_crispr_row_gate([
        _production_cpu_row(),
        _production_metal_row(),
        _production_validation_failure(),
        _production_sanson_ineligible(),
    ], failures)

    assert failures == []


def test_gpu_gate_rejects_production_crispr_metal_without_cpu_or_validation_evidence():
    gate = _load_gate()
    failures = []

    gate.production_crispr_row_gate([
        _production_metal_row(),
        _production_sanson_ineligible(),
    ], failures)

    assert any("missing production CRISPR CPU baseline" in failure for failure in failures)
    assert any("lacks validation-failure evidence" in failure for failure in failures)


def test_gpu_report_gate_requires_experimental_boundary(tmp_path):
    gate = _load_gate()
    report = tmp_path / "README.md"
    report.write_text("# Experimental GPU Acceleration Benchmark\n\nnot a production speed claim\n", encoding="utf-8")
    failures = []

    gate.report_gate(report, failures)

    assert failures == []


def test_production_crispr_report_gate_requires_cpu_authority_boundary(tmp_path):
    gate = _load_gate()
    report = tmp_path / "production_crispr_cpu_metal.md"
    report.write_text(
        "# Production CRISPR CPU vs Metal\n\n"
        "CPU remains the assignment authority. "
        "Metal is opt-in via `--backend gpu-metal-experimental`.\n\n"
        "Treat any Metal speedup as advisory until `metal_validation=passed` and guide-by-guide counts match the CPU shadow run.\n"
        "Do not use Metal for Sanson/Brunello-style multi-offset counting.\n"
        "The current Yusa rows are not count-identical to CPU. "
        "`--metal-validate` fails on Yusa today. "
        "Production `auto` staying on CPU is consistent with these measurements. "
        "Sanson/Brunello remains CPU-only.\n",
        encoding="utf-8",
    )
    failures = []

    gate.production_crispr_report_gate(report, failures)

    assert failures == []

    report.write_text("# Production CRISPR CPU vs Metal\n\nMetal is fast.\n", encoding="utf-8")
    failures = []

    gate.production_crispr_report_gate(report, failures)

    assert any("CPU remains the assignment authority" in failure for failure in failures)
