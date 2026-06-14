import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "scripts" / "generate_crispr_cpu_metal_report.py"


def _load_report():
    spec = importlib.util.spec_from_file_location("generate_crispr_cpu_metal_report", REPORT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_crispr_cpu_metal_report_documents_speedup_and_ineligible_rows(tmp_path, monkeypatch):
    report = _load_report()
    root = tmp_path / "repo"
    raw = root / "benchmarks" / "raw" / "crispr_cpu_metal.csv"
    out = root / "docs" / "benchmarks" / "gpu" / "production_crispr_cpu_metal.md"
    fig = root / "benchmarks" / "figures" / "crispr_cpu_metal_speedup.svg"
    raw.parent.mkdir(parents=True)
    fields = [
        "dataset", "backend", "metal_validate", "status", "records_per_sample",
        "n_samples", "total_reads", "n_targets", "guide_start", "guide_length", "k",
        "offset_mode", "auto_offset", "wall_seconds", "phase_total_seconds",
        "reads_per_sec", "backend_effective", "count_engine", "metal_validation",
        "exit_code", "count_match_cpu", "device", "notes",
    ]
    rows = [
        {
            "dataset": "mageck_yusa", "backend": "cpu", "metal_validate": "0", "status": "ok",
            "records_per_sample": "10000", "n_samples": "2", "total_reads": "20000",
            "n_targets": "87437", "guide_start": "23", "guide_length": "19", "k": "1",
            "offset_mode": "fixed", "auto_offset": "0", "wall_seconds": "0.40",
            "phase_total_seconds": "0.05", "reads_per_sec": "50000.0",
            "backend_effective": "cpu", "count_engine": "hamming_lookup_direct_single_offset",
            "metal_validation": "", "exit_code": "0", "count_match_cpu": "", "device": "test",
            "notes": "",
        },
        {
            "dataset": "mageck_yusa", "backend": "gpu-metal-experimental", "metal_validate": "0",
            "status": "ok", "records_per_sample": "10000", "n_samples": "2", "total_reads": "20000",
            "n_targets": "87437", "guide_start": "23", "guide_length": "19", "k": "1",
            "offset_mode": "fixed", "auto_offset": "0", "wall_seconds": "0.20",
            "phase_total_seconds": "0.15", "reads_per_sec": "100000.0",
            "backend_effective": "gpu-metal-experimental", "count_engine": "hamming_metal_seed_index",
            "metal_validation": "", "exit_code": "0", "count_match_cpu": "0", "device": "test",
            "notes": "",
        },
        {
            "dataset": "sanson_brunello", "backend": "gpu-metal-experimental", "metal_validate": "1",
            "status": "ineligible", "records_per_sample": "10000", "n_samples": "4", "total_reads": "",
            "n_targets": "", "guide_start": "20", "guide_length": "20", "k": "1",
            "offset_mode": "multi", "auto_offset": "20", "wall_seconds": "",
            "phase_total_seconds": "", "reads_per_sec": "", "backend_effective": "",
            "count_engine": "", "metal_validation": "", "exit_code": "", "count_match_cpu": "",
            "device": "test", "notes": "offset-mode multi blocks Metal",
        },
    ]
    with raw.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    monkeypatch.setattr(report, "ROOT", root)
    monkeypatch.setattr(report, "RAW", raw)
    monkeypatch.setattr(report, "OUT", out)
    monkeypatch.setattr(report, "FIG", fig)

    report.main()

    text = out.read_text(encoding="utf-8")
    assert "Production CRISPR CPU vs Metal" in text
    assert "2.00x" in text
    assert "count_match_cpu=0" in text
    assert "offset-mode multi blocks Metal" in text
    assert fig.exists()