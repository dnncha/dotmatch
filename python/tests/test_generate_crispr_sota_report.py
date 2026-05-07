import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "scripts" / "generate_crispr_sota_report.py"


def _load_report():
    spec = importlib.util.spec_from_file_location("generate_crispr_sota_report", REPORT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_crispr_sota_report_uses_relative_figure_links(tmp_path, monkeypatch):
    report = _load_report()
    root = tmp_path / "repo"
    raw = root / "benchmarks" / "raw"
    out_dir = root / "docs" / "benchmarks" / "crispr_sota"
    fig_dir = root / "benchmarks" / "figures"
    raw.mkdir(parents=True)
    (raw / "crispr_sota_repeated.csv").write_text(
        "tool,dataset_id,requested_records_per_sample,exit_code,reads_per_sec,seconds,peak_rss_kb,verified_per_read\n"
        "dotmatch_exact_k0,mageck_yusa,100000,0,10,1,1024,1\n"
        "dotmatch_hamming_k1,mageck_yusa,full,0,50,1,1024,1\n"
        "guide_counter_one_mismatch,mageck_yusa,full,0,100,1,1024,\n",
        encoding="utf-8",
    )
    (raw / "crispr_sota_edlib_validation.csv").write_text(
        "dataset,sample,checked_reads,mismatches,oracle_strategy,edlib_alignments,bounded_windows,fallback_windows\n"
        "mageck_yusa,plasmid,10,0,bounded_edlib_candidates,12,3,0\n",
        encoding="utf-8",
    )
    (raw / "crispr_sota_count_agreement_summary.csv").write_text("dataset,comparison,status\n", encoding="utf-8")

    monkeypatch.setattr(report, "ROOT", root)
    monkeypatch.setattr(report, "RAW", raw)
    monkeypatch.setattr(report, "OUT_DIR", out_dir)
    monkeypatch.setattr(report, "FIG_DIR", fig_dir)

    report.main()

    text = (out_dir / "README.md").read_text(encoding="utf-8")
    assert "](../../../benchmarks/figures/crispr_sota_throughput.svg)" in text
    assert "|dataset|sample|checked_reads|mismatches|oracle_strategy|edlib_alignments|bounded_windows|fallback_windows|" in text
    assert "## Full Hamming Speed Check" in text
    assert "|mageck_yusa|50.0|100.0|0.50|blocked|" in text
    assert str(root) not in text


def test_crispr_sota_report_aggregates_complete_full_sample_rows(tmp_path, monkeypatch):
    report = _load_report()
    root = tmp_path / "repo"
    raw = root / "benchmarks" / "raw"
    out_dir = root / "docs" / "benchmarks" / "crispr_sota"
    fig_dir = root / "benchmarks" / "figures"
    raw.mkdir(parents=True)
    rows = [
        "tool,dataset_id,requested_records_per_sample,run_level,sample_id,repeat,exit_code,n_reads,reads_per_sec,seconds,peak_rss_kb,verified_per_read",
    ]
    sample_reads = {
        "plasmid": 9821128,
        "RepA": 76471324,
        "RepB": 85301059,
        "RepC": 75356900,
    }
    for sample_id, reads in sample_reads.items():
        rows.append(f"dotmatch_hamming_k1,sanson_brunello,full,full_sample,{sample_id},1,0,{reads},{reads / 1.0:.1f},1.0,1024,1.0")
        rows.append(f"guide_counter_one_mismatch,sanson_brunello,full,full_sample,{sample_id},1,0,{reads},{reads / 2.0:.1f},2.0,2048,")
    (raw / "crispr_sota_repeated.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    (raw / "crispr_sota_edlib_validation.csv").write_text("dataset,sample,checked_reads,mismatches\n", encoding="utf-8")
    (raw / "crispr_sota_count_agreement_summary.csv").write_text("dataset,comparison,status\n", encoding="utf-8")

    monkeypatch.setattr(report, "ROOT", root)
    monkeypatch.setattr(report, "RAW", raw)
    monkeypatch.setattr(report, "OUT_DIR", out_dir)
    monkeypatch.setattr(report, "FIG_DIR", fig_dir)

    report.main()

    text = (out_dir / "README.md").read_text(encoding="utf-8")
    assert "|sanson_brunello|dotmatch_hamming_k1|full|1|61737602.8|4.0000|1.0|1.000|" in text
    assert "|sanson_brunello|61737602.8|30868801.4|2.00|pass|" in text
