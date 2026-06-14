import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "scripts" / "generate_public_crispr_report.py"


def _load_report():
    spec = importlib.util.spec_from_file_location("generate_public_crispr_report", REPORT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_public_crispr_report_includes_exact_and_hamming_speedups(tmp_path, monkeypatch):
    report = _load_report()
    root = tmp_path / "repo"
    raw = root / "benchmarks" / "raw" / "public_crispr_repeated.csv"
    out_dir = root / "docs" / "benchmarks" / "public_crispr"
    raw.parent.mkdir(parents=True)
    fields = [
        "tool", "version", "workflow", "semantics", "n_reads", "n_targets", "seconds",
        "reads_per_sec", "peak_rss_kb", "assigned_reads", "exact_reads", "corrected_reads",
        "ambiguous_reads", "rejected_reads", "overcount_reads", "candidates_per_read",
        "verified_per_read", "offset_mode", "hamming_index", "exit_code", "command",
        "dataset_id", "repeat", "requested_records_per_sample", "platform", "host", "python",
        "commit", "notes", "sample_id",
    ]
    rows = [
        {
            "tool": "dotmatch_exact_k0", "version": "local", "workflow": "public_crispr_yusa_small",
            "semantics": "exact_k0_no_errors", "n_reads": "20000", "n_targets": "87437",
            "seconds": "0.100000", "reads_per_sec": "200000.0", "peak_rss_kb": "120000",
            "assigned_reads": "17894", "exact_reads": "17894", "corrected_reads": "0",
            "ambiguous_reads": "0", "rejected_reads": "2106", "overcount_reads": "0",
            "candidates_per_read": "0.0", "verified_per_read": "0.9", "offset_mode": "best",
            "hamming_index": "query", "exit_code": "0", "command": "dotmatch exact",
            "dataset_id": "mageck_yusa", "repeat": "1", "requested_records_per_sample": "10000",
            "platform": "Darwin", "host": "arm64", "python": "3.11", "commit": "abc", "notes": "",
            "sample_id": "",
        },
        {
            "tool": "guide_counter_exact", "version": "0.1.3", "workflow": "public_crispr_yusa_small",
            "semantics": "exact_k0_no_errors", "n_reads": "20000", "n_targets": "87437",
            "seconds": "0.200000", "reads_per_sec": "100000.0", "peak_rss_kb": "500000",
            "assigned_reads": "17894", "exact_reads": "", "corrected_reads": "", "ambiguous_reads": "",
            "rejected_reads": "2106", "overcount_reads": "0", "candidates_per_read": "",
            "verified_per_read": "", "offset_mode": "", "hamming_index": "", "exit_code": "0",
            "command": "guide-counter exact", "dataset_id": "mageck_yusa", "repeat": "1",
            "requested_records_per_sample": "10000", "platform": "Darwin", "host": "arm64",
            "python": "3.11", "commit": "abc", "notes": "", "sample_id": "",
        },
        {
            "tool": "dotmatch_hamming_k1", "version": "local", "workflow": "public_crispr_yusa_small",
            "semantics": "hamming_k1_no_indels", "n_reads": "20000", "n_targets": "87437",
            "seconds": "0.090000", "reads_per_sec": "222222.2", "peak_rss_kb": "110000",
            "assigned_reads": "18376", "exact_reads": "17894", "corrected_reads": "482",
            "ambiguous_reads": "0", "rejected_reads": "1578", "overcount_reads": "0",
            "candidates_per_read": "0.9", "verified_per_read": "0.9", "offset_mode": "best",
            "hamming_index": "query", "exit_code": "0", "command": "dotmatch hamming",
            "dataset_id": "mageck_yusa", "repeat": "1", "requested_records_per_sample": "10000",
            "platform": "Darwin", "host": "arm64", "python": "3.11", "commit": "abc", "notes": "",
            "sample_id": "",
        },
    ]
    with raw.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    monkeypatch.setattr(report, "ROOT", root)
    monkeypatch.setattr(report, "RAW", raw)
    monkeypatch.setattr(report, "OUT_DIR", out_dir)

    report.main()

    text = (out_dir / "README.md").read_text(encoding="utf-8")
    assert "DotMatch Exact Count Speedup" in text
    assert "guide_counter_exact" in text
    assert "DotMatch Hamming Speedup" in text
