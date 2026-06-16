import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "scripts" / "generate_native_benchmark_report.py"


def _load_report():
    spec = importlib.util.spec_from_file_location("generate_native_benchmark_report", REPORT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_native_report_generates_without_pandas_or_matplotlib(tmp_path, monkeypatch):
    report = _load_report()
    root = tmp_path / "repo"
    out = root / "docs" / "benchmarks" / "native"
    raw = root / "benchmarks" / "raw"
    fig = root / "benchmarks" / "figures"
    sample = "\n".join([
        "tool,workload,error_mode,n_reads,n_targets,len,k,err,indel_rate,seconds,reads_per_sec,candidates_per_read,verified_per_read,peak_rss_kb,checksum,mismatches",
        "dotmatch_exact_batch,synthetic_barcode,exact,1000,4096,16,0,0.000,0.000,0.001,1000.0,1.00,1.00,10,1,0",
        "exact_hash_lookup,synthetic_barcode,exact,1000,4096,16,0,0.000,0.000,0.002,500.0,1.00,1.00,10,1,0",
        "dotmatch_indexed,synthetic_barcode,one_substitution,1000,4096,16,1,0.010,0.000,0.003,333.0,1.00,1.00,10,1,0",
        "edlib_native_scan,synthetic_barcode,one_substitution,1000,4096,16,1,0.010,0.000,0.300,3.0,4096.00,4096.00,10,1,0",
        "bk_tree,synthetic_barcode,one_substitution,1000,4096,16,1,0.010,0.000,0.010,100.0,20.00,20.00,10,1,0",
        "neighbor_lookup,synthetic_barcode,one_substitution,1000,4096,16,1,0.010,0.000,0.006,166.0,1.00,1.00,10,1,0",
        "dotmatch_indexed,synthetic_barcode,one_substitution,1000,4096,16,2,0.010,0.000,0.004,250.0,1.00,1.00,10,1,0",
        "edlib_native_scan,synthetic_barcode,one_substitution,1000,4096,16,2,0.010,0.000,0.300,20.0,4096.00,4096.00,10,1,0",
        "dotmatch_indexed,synthetic_barcode,one_insertion,1000,4096,16,2,0.000,0.062,0.004,250.0,1.00,8.00,10,1,0",
        "edlib_native_scan,synthetic_barcode,one_insertion,1000,4096,16,2,0.000,0.062,0.300,20.0,4096.00,4096.00,10,1,0",
        "dotmatch_indexed,synthetic_barcode,one_deletion,1000,4096,16,2,0.000,0.062,0.004,250.0,1.00,8.00,10,1,0",
        "edlib_native_scan,synthetic_barcode,one_deletion,1000,4096,16,2,0.000,0.062,0.300,20.0,4096.00,4096.00,10,1,0",
        "",
    ])

    def fake_run_command(cmd):
        if cmd == ["make", "build/bench_edlib_native"]:
            return ""
        return sample

    monkeypatch.setattr(report, "ROOT", root)
    monkeypatch.setattr(report, "OUT_DIR", out)
    monkeypatch.setattr(report, "RAW_DIR", raw)
    monkeypatch.setattr(report, "FIG_DIR", fig)
    monkeypatch.setattr(report, "REPORT_REPEATS", 1)
    monkeypatch.setattr(report, "REPORT_READS", 1000)
    monkeypatch.setattr(report, "plt", None)
    monkeypatch.setattr(report, "run_command", fake_run_command)

    report.main()

    text = (out / "README.md").read_text(encoding="utf-8")
    assert "Native Edlib Benchmark Report" in text
    assert "make native-exact-gate" in text
    assert "Gated Native Scaling Claims" in text
    assert "k=2 substitution indexed rows" in text
    assert "Levenshtein k=2 insertion/deletion rows" in text
    assert "| dotmatch_indexed | one_substitution | 4096 | 16 | 1 | 0.010 |" in text
    assert (raw / "native_edlib_assignment.csv").is_file()
    assert (raw / "native_edlib_assignment_summary.csv").is_file()
