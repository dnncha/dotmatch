import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "scripts" / "generate_gpu_report.py"


def _load_report():
    spec = importlib.util.spec_from_file_location("generate_gpu_report", REPORT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_gpu_report_generates_relative_figure_and_speed_table(tmp_path, monkeypatch):
    report = _load_report()
    root = tmp_path / "repo"
    raw = root / "benchmarks" / "raw"
    out_dir = root / "docs" / "benchmarks" / "gpu"
    fig_dir = root / "benchmarks" / "figures"
    raw.mkdir(parents=True)
    csv_path = raw / "gpu_acceleration.csv"
    crispr_path = raw / "gpu_crispr.csv"
    fields = [
        "tool", "backend", "status", "workload", "n_reads", "n_targets", "len", "k",
        "error_rate", "prep_seconds", "seconds", "total_seconds", "reads_per_sec",
        "total_reads_per_sec", "pairs_per_sec", "checksum", "mismatches", "device", "notes",
    ]
    rows = [
        {
            "tool": "dotmatch_cpu_index", "backend": "cpu", "status": "ok",
            "workload": "synthetic_hamming", "n_reads": "1000", "n_targets": "96",
            "len": "20", "k": "1", "error_rate": "0.01", "prep_seconds": "0.1",
            "seconds": "1.0", "total_seconds": "1.1", "reads_per_sec": "1000",
            "total_reads_per_sec": "909.1", "pairs_per_sec": "96000", "checksum": "10",
            "mismatches": "0", "device": "cpu", "notes": "",
        },
        {
            "tool": "dotmatch_gpu_metal", "backend": "metal", "status": "ok",
            "workload": "synthetic_hamming", "n_reads": "1000", "n_targets": "96",
            "len": "20", "k": "1", "error_rate": "0.01", "prep_seconds": "0.1",
            "seconds": "0.5", "total_seconds": "0.6", "reads_per_sec": "2000",
            "total_reads_per_sec": "1666.7", "pairs_per_sec": "192000", "checksum": "10",
            "mismatches": "0", "device": "Apple test GPU", "notes": "",
        },
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    crispr_fields = [
        "tool", "backend", "status", "workload", "total_reads", "packable_reads",
        "n_targets", "target_start", "target_length", "k", "input_seconds",
        "prep_seconds", "seconds", "total_seconds", "reads_per_sec",
        "total_reads_per_sec", "assigned_unique", "assigned_exact",
        "assigned_corrected", "ambiguous", "unmatched", "invalid_windows",
        "non_acgt_windows", "skipped_targets", "checksum", "mismatches",
        "count_delta", "device", "notes",
    ]
    crispr_rows = [
        {
            "tool": "dotmatch_cpu_index", "backend": "cpu", "status": "ok",
            "workload": "public_crispr_yusa_hamming", "total_reads": "1000",
            "packable_reads": "990", "n_targets": "96", "target_start": "23",
            "target_length": "19", "k": "1", "input_seconds": "1.0",
            "prep_seconds": "0.1", "seconds": "1.0", "total_seconds": "2.1",
            "reads_per_sec": "1000", "total_reads_per_sec": "476.2",
            "assigned_unique": "900", "assigned_exact": "850",
            "assigned_corrected": "50", "ambiguous": "10", "unmatched": "80",
            "invalid_windows": "10", "non_acgt_windows": "0", "skipped_targets": "0",
            "checksum": "10", "mismatches": "0", "count_delta": "0",
            "device": "cpu", "notes": "",
        },
        {
            "tool": "dotmatch_gpu_metal", "backend": "metal", "status": "ok",
            "workload": "public_crispr_yusa_hamming", "total_reads": "1000",
            "packable_reads": "990", "n_targets": "96", "target_start": "23",
            "target_length": "19", "k": "1", "input_seconds": "1.0",
            "prep_seconds": "0.1", "seconds": "0.5", "total_seconds": "1.6",
            "reads_per_sec": "2000", "total_reads_per_sec": "625.0",
            "assigned_unique": "900", "assigned_exact": "850",
            "assigned_corrected": "50", "ambiguous": "10", "unmatched": "80",
            "invalid_windows": "10", "non_acgt_windows": "0", "skipped_targets": "0",
            "checksum": "10", "mismatches": "0", "count_delta": "0",
            "device": "Apple test GPU", "notes": "",
        },
    ]
    with crispr_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=crispr_fields)
        writer.writeheader()
        writer.writerows(crispr_rows)

    monkeypatch.setattr(report, "ROOT", root)
    monkeypatch.setattr(report, "RAW", csv_path)
    monkeypatch.setattr(report, "CRISPR_RAW", crispr_path)
    monkeypatch.setattr(report, "OUT_DIR", out_dir)
    monkeypatch.setattr(report, "FIG_DIR", fig_dir)

    report.main()

    text = (out_dir / "README.md").read_text(encoding="utf-8")
    assert "not a production speed claim" in text
    assert "](../../../benchmarks/figures/gpu_metal_speedup.svg)" in text
    assert "](../../../benchmarks/figures/gpu_crispr_metal_speedup.svg)" in text
    assert "| 1000 | 96 | 20 | 1 | 2000.0 | 1000.0 | 2.00x | 1.83x | 0 |" in text
    assert "| 1000 | 990 | 96 | 23 | 19 | 2000.0 | 1000.0 | 2.00x | 1.31x | 0 | 0 |" in text
