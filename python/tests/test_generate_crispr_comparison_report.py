import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "scripts" / "generate_crispr_comparison_report.py"


def _load_report():
    spec = importlib.util.spec_from_file_location("generate_crispr_comparison_report", REPORT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_crispr_comparison_report_uses_relative_figure_links(tmp_path, monkeypatch):
    report = _load_report()
    root = tmp_path / "repo"
    raw = root / "benchmarks" / "raw"
    out_dir = root / "docs" / "benchmarks" / "crispr_comparison"
    fig_dir = root / "benchmarks" / "figures"
    raw.mkdir(parents=True)
    (raw / "crispr_comparison_repeated.csv").write_text(
        "tool,dataset_id,requested_records_per_sample,exit_code,reads_per_sec,seconds,peak_rss_kb,verified_per_read\n"
        "dotmatch_exact_k0,mageck_yusa,100000,0,10,1,1024,1\n"
        "dotmatch_hamming_k1,mageck_yusa,full,0,50,1,1024,1\n"
        "guide_counter_one_mismatch,mageck_yusa,full,0,100,1,1024,\n",
        encoding="utf-8",
    )
    (raw / "crispr_comparison_edlib_validation.csv").write_text(
        "dataset,sample,checked_reads,mismatches,oracle_strategy,edlib_alignments,bounded_windows,fallback_windows\n"
        "mageck_yusa,plasmid,10,0,bounded_edlib_candidates,12,3,0\n",
        encoding="utf-8",
    )
    (raw / "crispr_comparison_count_agreement_summary.csv").write_text("dataset,comparison,status\n", encoding="utf-8")

    monkeypatch.setattr(report, "ROOT", root)
    monkeypatch.setattr(report, "RAW", raw)
    monkeypatch.setattr(report, "OUT_DIR", out_dir)
    monkeypatch.setattr(report, "FIG_DIR", fig_dir)

    report.main()

    text = (out_dir / "README.md").read_text(encoding="utf-8")
    assert "](../../../benchmarks/figures/crispr_comparison_throughput.svg)" in text
    assert "|dataset|sample|checked_reads|mismatches|oracle_strategy|edlib_alignments|bounded_windows|fallback_windows|" in text
    assert "## Full Hamming Guide-Counter Ratio" in text
    assert "|mageck_yusa|50.0|100.0|0.50|reported|" in text
    assert str(root) not in text


def test_crispr_comparison_report_aggregates_complete_full_sample_rows(tmp_path, monkeypatch):
    report = _load_report()
    root = tmp_path / "repo"
    raw = root / "benchmarks" / "raw"
    out_dir = root / "docs" / "benchmarks" / "crispr_comparison"
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
    (raw / "crispr_comparison_repeated.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    (raw / "crispr_comparison_edlib_validation.csv").write_text("dataset,sample,checked_reads,mismatches\n", encoding="utf-8")
    (raw / "crispr_comparison_count_agreement_summary.csv").write_text("dataset,comparison,status\n", encoding="utf-8")

    monkeypatch.setattr(report, "ROOT", root)
    monkeypatch.setattr(report, "RAW", raw)
    monkeypatch.setattr(report, "OUT_DIR", out_dir)
    monkeypatch.setattr(report, "FIG_DIR", fig_dir)

    report.main()

    text = (out_dir / "README.md").read_text(encoding="utf-8")
    assert "|sanson_brunello|dotmatch_hamming_k1|full|1|61737602.8|4.0000|1.0|1.000|" in text
    assert "|sanson_brunello|61737602.8|30868801.4|2.00|reported|" in text
    svg = (fig_dir / "crispr_comparison_throughput.svg").read_text(encoding="utf-8")
    assert "sanson_brunello dotmatch_hamming_k1 full FASTQs" in svg
    assert "bar-full" in svg


def test_crispr_comparison_report_has_explicit_guide_counter_style_public_data_lane(tmp_path, monkeypatch):
    report = _load_report()
    root = tmp_path / "repo"
    raw = root / "benchmarks" / "raw"
    out_dir = root / "docs" / "benchmarks" / "crispr_comparison"
    fig_dir = root / "benchmarks" / "figures"
    raw.mkdir(parents=True)
    (raw / "crispr_comparison_repeated.csv").write_text(
        "tool,dataset_id,requested_records_per_sample,exit_code,reads_per_sec,seconds,peak_rss_kb,verified_per_read\n"
        "dotmatch_hamming_k1,mageck_yusa,100000,0,200,1,1024,1\n"
        "guide_counter_one_mismatch,mageck_yusa,100000,0,100,1,2048,\n"
        "dotmatch_levenshtein_k1,mageck_yusa,100000,0,50,1,1024,2\n",
        encoding="utf-8",
    )
    (raw / "crispr_comparison_edlib_validation.csv").write_text("dataset,sample,checked_reads,mismatches\n", encoding="utf-8")
    (raw / "crispr_comparison_count_agreement_summary.csv").write_text(
        "dataset,comparison,status,total_left,total_right,total_delta,differing_guides,max_abs_delta,pearson,spearman\n"
        "mageck_yusa,mageck_yusa:dotmatch_hamming_vs_guide_counter,ok,200000,199000,1000,25,4,0.99,0.98\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(report, "ROOT", root)
    monkeypatch.setattr(report, "RAW", raw)
    monkeypatch.setattr(report, "OUT_DIR", out_dir)
    monkeypatch.setattr(report, "FIG_DIR", fig_dir)

    report.main()

    text = (out_dir / "README.md").read_text(encoding="utf-8")
    assert "## Guide-Counter-Style Public Paper-Data Lane" in text
    assert "DotMatch `dotmatch_hamming_k1` versus `guide_counter_one_mismatch`" in text
    assert "|dataset|records_per_sample|dotmatch_hamming_reads_per_sec|guide_counter_reads_per_sec|speedup|count_agreement_status|count_total_delta|semantics|" in text
    assert "|mageck_yusa|100000|200.0|100.0|2.00|ok|1000|one mismatch, no indels|" in text
    assert "dotmatch_levenshtein_k1" not in text.split("## Guide-Counter-Style Public Paper-Data Lane", 1)[1].split("## Full Hamming Guide-Counter Ratio", 1)[0]


def test_crispr_comparison_report_includes_backend_optimizer_artifact(tmp_path, monkeypatch):
    report = _load_report()
    root = tmp_path / "repo"
    raw = root / "benchmarks" / "raw"
    out_dir = root / "docs" / "benchmarks" / "crispr_comparison"
    fig_dir = root / "benchmarks" / "figures"
    raw.mkdir(parents=True)
    (raw / "crispr_comparison_repeated.csv").write_text(
        "tool,dataset_id,requested_records_per_sample,exit_code,reads_per_sec,seconds,peak_rss_kb\n"
        "dotmatch_hamming_k1,sanson_brunello,100000,0,200,1,1024\n",
        encoding="utf-8",
    )
    (raw / "crispr_comparison_edlib_validation.csv").write_text("dataset,sample,checked_reads,mismatches\n", encoding="utf-8")
    (raw / "crispr_comparison_count_agreement_summary.csv").write_text("dataset,comparison,status\n", encoding="utf-8")
    optimizer = raw / "crispr_sanson_brunello_backend_optimization_atlas.json"
    optimizer.write_text(
        """{
  "optimizer": "local_benchmark_informed_scorer_v1",
  "authority": "cpu",
  "selected_backend": "cpu",
  "candidate_backend": "gpu-metal-experimental",
  "recommendation": "gpu_candidate_requires_cpu_validation",
  "expected_speedup_band": "1.5-3x",
  "estimated_total_speedup": 2.6
}
""",
        encoding="utf-8",
    )

    monkeypatch.setattr(report, "ROOT", root)
    monkeypatch.setattr(report, "RAW", raw)
    monkeypatch.setattr(report, "OUT_DIR", out_dir)
    monkeypatch.setattr(report, "FIG_DIR", fig_dir)
    monkeypatch.setattr(report, "OPTIMIZER_ARTIFACTS", {"sanson_brunello": optimizer})

    report.main()

    text = (out_dir / "README.md").read_text(encoding="utf-8")
    assert "## Backend Optimizer" in text
    assert "|sanson_brunello|local_benchmark_informed_scorer_v1|cpu|cpu|gpu-metal-experimental|gpu_candidate_requires_cpu_validation|1.5-3x|2.6|" in text
