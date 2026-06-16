import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_native_exact_gate.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_native_exact_gate", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


HEADER = (
    "tool,workload,error_mode,n_reads,n_targets,len,k,err,indel_rate,"
    "reads_per_sec,verified_per_read,mismatches\n"
)


def _row(
    tool: str,
    *,
    rps: float = 100.0,
    mismatches: str = "0",
    n_targets: int = 4096,
    k: int = 0,
    verified_per_read: float = 0.01,
) -> str:
    return (
        f"{tool},synthetic_barcode,exact,1000,{n_targets},16,{k},0.000,0.000,"
        f"{rps},{verified_per_read},{mismatches}\n"
    )


def _valid_k1_rows() -> str:
    return (
        _row("dotmatch_indexed", rps=5000.0, k=1, verified_per_read=0.02)
        + _row("edlib_native_scan", rps=100.0, k=1, verified_per_read=4096.0)
        + _row("bk_tree", rps=1000.0, k=1, verified_per_read=20.0)
        + _row("neighbor_lookup", rps=2500.0, k=1, verified_per_read=0.02)
    )


def _valid_k2_rows() -> str:
    return (
        _row("dotmatch_indexed", rps=4000.0, k=2, verified_per_read=8.0, n_targets=4096)
        .replace(",exact,", ",one_insertion,")
        + _row("edlib_native_scan", rps=100.0, k=2, verified_per_read=4096.0, n_targets=4096)
        .replace(",exact,", ",one_insertion,")
        + _row("dotmatch_indexed", rps=4200.0, k=2, verified_per_read=9.0, n_targets=4096)
        .replace(",exact,", ",one_deletion,")
        + _row("edlib_native_scan", rps=100.0, k=2, verified_per_read=4096.0, n_targets=4096)
        .replace(",exact,", ",one_deletion,")
    )


def _valid_k2_substitution_rows() -> str:
    return (
        _row("dotmatch_indexed", rps=2200.0, k=2, verified_per_read=1.0, n_targets=4096)
        .replace(",exact,", ",one_substitution,")
        + _row("edlib_native_scan", rps=100.0, k=2, verified_per_read=4096.0, n_targets=4096)
        .replace(",exact,", ",one_substitution,")
    )


def test_native_exact_gate_accepts_comparator_backed_rows(tmp_path):
    checker = _load_checker()
    csv = tmp_path / "native.csv"
    csv.write_text(
        HEADER
        + _row("dotmatch_exact_batch", rps=120.0)
        + _row("exact_hash_lookup", rps=100.0)
        + _row("edlib_native_scan", rps=1.0)
        + _valid_k1_rows()
        + _valid_k2_substitution_rows()
        + _valid_k2_rows(),
        encoding="utf-8",
    )

    result = checker.audit(csv)

    assert result.failures == []
    assert any("median ratio" in item for item in result.passed)
    assert any("large-library" in item for item in result.passed)


def test_native_exact_gate_requires_direct_rows(tmp_path):
    checker = _load_checker()
    csv = tmp_path / "native.csv"
    csv.write_text(
        HEADER + _row("exact_hash_lookup") + _row("edlib_native_scan") + _valid_k1_rows() + _valid_k2_rows(),
        encoding="utf-8",
    )

    result = checker.audit(csv)

    assert any("dotmatch_exact_direct or dotmatch_exact_batch" in failure for failure in result.failures)


def test_native_exact_gate_rejects_mismatches_and_missing_hash_pair(tmp_path):
    checker = _load_checker()
    csv = tmp_path / "native.csv"
    csv.write_text(
        HEADER
        + _row("dotmatch_exact_direct", mismatches="1")
        + _row("edlib_native_scan")
        + _valid_k1_rows()
        + _valid_k2_substitution_rows()
        + _valid_k2_rows(),
        encoding="utf-8",
    )

    result = checker.audit(csv)

    assert any("exact_hash_lookup k=0 rows" in failure for failure in result.failures)
    assert any("assignment mismatches" in failure for failure in result.failures)


def test_native_exact_gate_rejects_large_library_hash_regression(tmp_path):
    checker = _load_checker()
    csv = tmp_path / "native.csv"
    csv.write_text(
        HEADER
        + _row("dotmatch_exact_batch", rps=80.0)
        + _row("exact_hash_lookup", rps=100.0)
        + _row("edlib_native_scan", rps=1.0)
        + _valid_k1_rows()
        + _valid_k2_substitution_rows()
        + _valid_k2_rows(),
        encoding="utf-8",
    )

    result = checker.audit(csv)

    assert any("large-library native exact rows must beat exact_hash_lookup" in failure for failure in result.failures)


def test_native_gate_requires_k1_indexed_rows(tmp_path):
    checker = _load_checker()
    csv = tmp_path / "native.csv"
    csv.write_text(
        HEADER
        + _row("dotmatch_exact_batch", rps=120.0)
        + _row("exact_hash_lookup", rps=100.0)
        + _row("edlib_native_scan", rps=1.0),
        encoding="utf-8",
    )

    result = checker.audit(csv)

    assert any("dotmatch_indexed k=1 rows" in failure for failure in result.failures)


def test_native_gate_rejects_k1_edlib_or_baseline_regression(tmp_path):
    checker = _load_checker()
    csv = tmp_path / "native.csv"
    csv.write_text(
        HEADER
        + _row("dotmatch_exact_batch", rps=120.0)
        + _row("exact_hash_lookup", rps=100.0)
        + _row("edlib_native_scan", rps=1.0)
        + _row("dotmatch_indexed", rps=900.0, k=1, verified_per_read=0.02, n_targets=4096)
        + _row("edlib_native_scan", rps=100.0, k=1, verified_per_read=4096.0)
        + _row("neighbor_lookup", rps=1000.0, k=1, verified_per_read=0.02)
        + _valid_k2_substitution_rows()
        + _valid_k2_rows(),
        encoding="utf-8",
    )

    result = checker.audit(csv)

    assert any("exhaustive Edlib scan" in failure for failure in result.failures)
    assert any("large-library" in failure and "best BK-tree/neighbor baseline" in failure for failure in result.failures)


def test_native_gate_rejects_k1_high_candidate_verification(tmp_path):
    checker = _load_checker()
    csv = tmp_path / "native.csv"
    csv.write_text(
        HEADER
        + _row("dotmatch_exact_batch", rps=120.0)
        + _row("exact_hash_lookup", rps=100.0)
        + _row("edlib_native_scan", rps=1.0)
        + _row("dotmatch_indexed", rps=5000.0, k=1, verified_per_read=1.20)
        + _row("edlib_native_scan", rps=100.0, k=1, verified_per_read=4096.0)
        + _row("neighbor_lookup", rps=2500.0, k=1, verified_per_read=0.02)
        + _valid_k2_substitution_rows()
        + _valid_k2_rows(),
        encoding="utf-8",
    )

    result = checker.audit(csv)

    assert any("no more than 1.05 candidates/read" in failure for failure in result.failures)


def test_native_gate_requires_k2_indel_rows(tmp_path):
    checker = _load_checker()
    csv = tmp_path / "native.csv"
    csv.write_text(
        HEADER
        + _row("dotmatch_exact_batch", rps=120.0)
        + _row("exact_hash_lookup", rps=100.0)
        + _row("edlib_native_scan", rps=1.0)
        + _valid_k1_rows()
        + _valid_k2_substitution_rows(),
        encoding="utf-8",
    )

    result = checker.audit(csv)

    assert any("Levenshtein k=2 insertion/deletion rows" in failure for failure in result.failures)


def test_native_gate_rejects_k2_speed_or_candidate_regression(tmp_path):
    checker = _load_checker()
    csv = tmp_path / "native.csv"
    k2_bad = (
        _row("dotmatch_indexed", rps=400.0, k=2, verified_per_read=30.0, n_targets=4096)
        .replace(",exact,", ",one_insertion,")
        + _row("edlib_native_scan", rps=100.0, k=2, verified_per_read=4096.0, n_targets=4096)
        .replace(",exact,", ",one_insertion,")
        + _row("dotmatch_indexed", rps=4200.0, k=2, verified_per_read=9.0, n_targets=4096)
        .replace(",exact,", ",one_deletion,")
        + _row("edlib_native_scan", rps=100.0, k=2, verified_per_read=4096.0, n_targets=4096)
        .replace(",exact,", ",one_deletion,")
    )
    csv.write_text(
        HEADER
        + _row("dotmatch_exact_batch", rps=120.0)
        + _row("exact_hash_lookup", rps=100.0)
        + _row("edlib_native_scan", rps=1.0)
        + _valid_k1_rows()
        + _valid_k2_substitution_rows()
        + k2_bad,
        encoding="utf-8",
    )

    result = checker.audit(csv)

    assert any("k=2 indexed rows must beat exhaustive Edlib" in failure for failure in result.failures)
    assert any("no more than 25 candidates/read" in failure for failure in result.failures)


def test_native_gate_rejects_k2_substitution_speed_or_candidate_regression(tmp_path):
    checker = _load_checker()
    csv = tmp_path / "native.csv"
    k2_substitution_bad = (
        _row("dotmatch_indexed", rps=700.0, k=2, verified_per_read=1.20, n_targets=4096)
        .replace(",exact,", ",one_substitution,")
        + _row("edlib_native_scan", rps=100.0, k=2, verified_per_read=4096.0, n_targets=4096)
        .replace(",exact,", ",one_substitution,")
    )
    csv.write_text(
        HEADER
        + _row("dotmatch_exact_batch", rps=120.0)
        + _row("exact_hash_lookup", rps=100.0)
        + _row("edlib_native_scan", rps=1.0)
        + _valid_k1_rows()
        + k2_substitution_bad
        + _valid_k2_rows(),
        encoding="utf-8",
    )

    result = checker.audit(csv)

    assert any("k=2 substitution rows must beat exhaustive Edlib" in failure for failure in result.failures)
    assert any("k=2 substitution rows must verify no more than 1.05 candidates/read" in failure for failure in result.failures)


def test_native_report_gate_requires_gated_scaling_rows(tmp_path):
    checker = _load_checker()
    result = checker.AuditResult()
    report = tmp_path / "README.md"
    report.write_text(
        "# Native\n\n"
        "## Gated Native Scaling Claims\n\n"
        "| claim | large_library_rows | min_speedup_vs_edlib | median_speedup_vs_edlib | max_verified_per_read | min_speedup_required | max_verified_required |\n"
        "| --- | --- | --- | --- | --- | --- | --- |\n"
        "| k=1 substitution indexed rows | 36 | 531.68 | 698.44 | 1.00 | 10.00 | 1.05 |\n"
        "not end-to-end workflow speed claims\n"
        "This remains scoped to packed A/C/G/T fixed-window assignment up to 32 bases\n",
        encoding="utf-8",
    )

    checker.report_gate(
        report,
        [
            {
                "claim": "k=1 substitution indexed rows",
                "large_library_rows": "36",
                "min_speedup_vs_edlib": "531.68",
                "median_speedup_vs_edlib": "698.44",
                "max_verified_per_read": "1.00",
                "min_speedup_required": "8.00",
                "max_verified_required": "1.05",
            },
            {
                "claim": "k=2 substitution indexed rows",
                "large_library_rows": "36",
                "min_speedup_vs_edlib": "10.88",
                "median_speedup_vs_edlib": "20.89",
                "max_verified_per_read": "1.00",
                "min_speedup_required": "10.00",
                "max_verified_required": "1.05",
            },
        ],
        result,
    )

    assert any("k=2 substitution indexed rows" in failure for failure in result.failures)


def test_native_report_gate_accepts_matching_gated_scaling_rows(tmp_path):
    checker = _load_checker()
    result = checker.AuditResult()
    report = tmp_path / "README.md"
    report.write_text(
        "# Native\n\n"
        "## Gated Native Scaling Claims\n\n"
        "| claim | large_library_rows | min_speedup_vs_edlib | median_speedup_vs_edlib | max_verified_per_read | min_speedup_required | max_verified_required |\n"
        "| --- | --- | --- | --- | --- | --- | --- |\n"
        "| Levenshtein k=2 insertion/deletion rows | 18 | 9.76 | 15.83 | 1.00 | 8.00 | 25.00 |\n"
        "not end-to-end workflow speed claims\n"
        "This remains scoped to packed A/C/G/T fixed-window assignment up to 32 bases\n",
        encoding="utf-8",
    )

    checker.report_gate(
        report,
        [
            {
                "claim": "Levenshtein k=2 insertion/deletion rows",
                "large_library_rows": "18",
                "min_speedup_vs_edlib": "9.76",
                "median_speedup_vs_edlib": "15.83",
                "max_verified_per_read": "1.00",
                "min_speedup_required": "8.00",
                "max_verified_required": "25.00",
            },
        ],
        result,
    )

    assert result.failures == []
