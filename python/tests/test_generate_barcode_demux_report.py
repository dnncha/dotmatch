import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "scripts" / "generate_barcode_demux_report.py"


def _load_report():
    spec = importlib.util.spec_from_file_location("generate_barcode_demux_report", REPORT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _row(tool: str, workflow: str, reads_per_sec: str, *, metric: str = "hamming", engine: str = "exact_lookup") -> dict[str, str]:
    return {
        "tool": tool,
        "workflow": workflow,
        "semantics": f"fixed_window_{metric}",
        "repeat": "1",
        "n_reads": "1000",
        "n_barcodes": "96",
        "barcode_length": "8",
        "k": "1" if metric != "exact" else "0",
        "metric": metric,
        "assignment_engine": engine,
        "seconds": "0.010000",
        "reads_per_sec": reads_per_sec,
        "peak_rss_kb": "1000",
        "assigned_reads": "900",
        "exact_reads": "850",
        "corrected_reads": "50",
        "ambiguous_reads": "0",
        "unmatched_reads": "100",
        "verified_per_read": "1.00",
        "exit_code": "0",
        "command": "fixture",
    }


def test_barcode_report_includes_gated_real_data_speedups(tmp_path, monkeypatch):
    report = _load_report()
    root = tmp_path / "repo"
    raw = root / "benchmarks" / "raw" / "barcode_demux.csv"
    out_dir = root / "docs" / "benchmarks" / "barcode_demux"
    fig_dir = root / "benchmarks" / "figures"
    raw.parent.mkdir(parents=True)
    rows = [
        _row("dotmatch_demux", "real_srp009896_inline_barcode", "3000.0", metric="exact", engine="exact_prefix_lookup"),
        _row("cutadapt_demux", "real_srp009896_inline_barcode", "500.0", metric="exact", engine="cutadapt"),
        _row("hash_splitter_exact", "real_srp009896_inline_barcode", "1000.0", metric="exact", engine="exact_hash"),
        _row("dotmatch_demux", "real_srp009896_inline_barcode_fixed8_k1", "5000.0", engine="hamming_k1_lookup_direct"),
        _row("cutadapt_demux", "real_srp009896_inline_barcode_fixed8_k1", "1000.0", engine="cutadapt"),
        _row("hamming_radius_splitter", "real_srp009896_inline_barcode_fixed8_k1", "400.0", engine="transparent_hamming"),
        _row("dotmatch_demux", "synthetic_levenshtein_one_edit_fixture", "2500.0", metric="levenshtein", engine="levenshtein_k1_lookup_direct"),
        _row("levenshtein_radius_splitter", "synthetic_levenshtein_one_edit_fixture", "250.0", metric="levenshtein", engine="transparent_levenshtein"),
    ]
    with raw.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    monkeypatch.setattr(report, "ROOT", root)
    monkeypatch.setattr(report, "RAW", raw)
    monkeypatch.setattr(report, "OUT_DIR", out_dir)
    monkeypatch.setattr(report, "FIG_DIR", fig_dir)

    report.main()

    text = (out_dir / "README.md").read_text(encoding="utf-8")
    assert "## Gated Real-Data Speedups" in text
    assert "| real_srp009896_inline_barcode | cutadapt_demux | 3000.0 | 500.0 | 6.00x | 5.00x |" in text
    assert "| real_srp009896_inline_barcode | hash_splitter_exact | 3000.0 | 1000.0 | 3.00x | 3.00x |" in text
    assert "| real_srp009896_inline_barcode_fixed8_k1 | hamming_radius_splitter | 5000.0 | 400.0 | 12.50x | 12.00x |" in text
    assert "levenshtein_k1_lookup_direct" in text
    assert "not public real-data evidence" in text
    assert "Levenshtein indel lane is fixture evidence" in text
    assert (fig_dir / "barcode_demux_throughput.svg").is_file()
