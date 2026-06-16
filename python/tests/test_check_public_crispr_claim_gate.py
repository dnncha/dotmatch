import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "scripts" / "check_public_crispr_claim_gate.py"


def _load_gate():
    spec = importlib.util.spec_from_file_location("check_public_crispr_claim_gate", GATE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _row(tool, requested, repeat, reads_per_sec="100.0", verified_per_read="1.0"):
    return {
        "tool": tool,
        "requested_records_per_sample": str(requested),
        "repeat": str(repeat),
        "exit_code": "0",
        "reads_per_sec": str(reads_per_sec),
        "verified_per_read": verified_per_read,
    }


def test_public_gate_requires_guide_counter_rows_without_speed_superiority_claim():
    gate = _load_gate()
    rows = []
    for repeat in range(1, 3):
        rows.extend(
            [
                _row("dotmatch_exact_k0", 100000, repeat),
                _row("dotmatch_hamming_k1", 100000, repeat, reads_per_sec="90.0"),
                _row("dotmatch_levenshtein_k1", 100000, repeat, verified_per_read="2.0"),
                _row("mageck_count_exact", 100000, repeat),
                _row("guide_counter_one_mismatch", 100000, repeat, reads_per_sec="120.0"),
                _row("guide_counter_exact", 100000, repeat, reads_per_sec="80.0"),
            ]
        )
    failures = []

    gate.repeated_gate(rows, min_records=100000, min_repeats=2, require_guide_counter=True, failures=failures)

    assert failures == []


def test_public_gate_still_requires_guide_counter_when_requested():
    gate = _load_gate()
    rows = []
    for repeat in range(1, 3):
        rows.extend(
            [
                _row("dotmatch_exact_k0", 100000, repeat),
                _row("dotmatch_hamming_k1", 100000, repeat),
                _row("dotmatch_levenshtein_k1", 100000, repeat, verified_per_read="2.0"),
                _row("mageck_count_exact", 100000, repeat),
            ]
        )
    failures = []

    gate.repeated_gate(rows, min_records=100000, min_repeats=2, require_guide_counter=True, failures=failures)

    assert any("guide_counter_one_mismatch needs >= 2 successful repeats" in f for f in failures)


def test_public_validation_gate_requires_bounded_zero_fallback_when_recorded():
    gate = _load_gate()
    rows = [
        {
            "sample": "plasmid",
            "checked_reads": "1000",
            "mismatches": "0",
            "oracle_strategy": "bounded_edlib_candidates",
            "edlib_alignments": "100",
            "bounded_windows": "1000",
            "fallback_windows": "50",
        }
    ]
    failures = []

    gate.validation_gate(rows, min_checked=1000, failures=failures)

    assert failures == []


def test_public_validation_gate_rejects_fallback_or_unbounded_oracle_when_recorded():
    gate = _load_gate()
    rows = [
        {
            "sample": "plasmid",
            "checked_reads": "1000",
            "mismatches": "0",
            "oracle_strategy": "full_edlib_scan",
            "edlib_alignments": "100",
            "bounded_windows": "0",
            "fallback_windows": "51",
        }
    ]
    failures = []

    gate.validation_gate(rows, min_checked=1000, failures=failures)

    assert any("bounded_edlib_candidates" in failure for failure in failures)
    assert any("bounded_windows" in failure for failure in failures)
    assert any("fallback_windows exceeds 5%" in failure for failure in failures)


def test_public_report_gate_requires_speedup_rows(tmp_path):
    gate = _load_gate()
    failures = []
    report = tmp_path / "README.md"
    report.write_text(
        "# Public\n\n"
        "## DotMatch Hamming Speedup\n\n"
        "## DotMatch Exact Count Speedup\n\n"
        "## Edlib Oracle Validation\n\n"
        "bounded Edlib validation\n"
        "guide-counter is fast\n"
        "DotMatch assigns at most one target per read\n",
        encoding="utf-8",
    )
    rows = [
        _row("dotmatch_hamming_k1", 100000, 1, reads_per_sec="200.0"),
        _row("guide_counter_one_mismatch", 100000, 1, reads_per_sec="100.0"),
        _row("dotmatch_exact_k0", 100000, 1, reads_per_sec="300.0"),
        _row("mageck_count_exact", 100000, 1, reads_per_sec="30.0"),
        _row("guide_counter_exact", 100000, 1, reads_per_sec="150.0"),
    ]

    gate.report_gate(report, rows, require_guide_counter=True, failures=failures)

    assert any("Hamming speedup row" in failure for failure in failures)
    assert any("exact speedup row" in failure for failure in failures)


def test_public_report_gate_accepts_matching_speedup_rows(tmp_path):
    gate = _load_gate()
    failures = []
    report = tmp_path / "README.md"
    report.write_text(
        "# Public\n\n"
        "## DotMatch Hamming Speedup\n\n"
        "| baseline | records_per_sample | dotmatch_hamming_reads_per_sec | baseline_reads_per_sec | speedup |\n"
        "| --- | ---: | ---: | ---: | ---: |\n"
        "| guide_counter_one_mismatch | 100000 | 200.0 | 100.0 | 2.00x |\n\n"
        "## DotMatch Exact Count Speedup\n\n"
        "| baseline | records_per_sample | dotmatch_exact_reads_per_sec | baseline_reads_per_sec | speedup |\n"
        "| --- | ---: | ---: | ---: | ---: |\n"
        "| guide_counter_exact | 100000 | 300.0 | 150.0 | 2.00x |\n"
        "| mageck_count_exact | 100000 | 300.0 | 30.0 | 10.00x |\n\n"
        "## Edlib Oracle Validation\n\n"
        "bounded Edlib validation\n"
        "guide-counter is fast\n"
        "DotMatch assigns at most one target per read\n",
        encoding="utf-8",
    )
    rows = [
        _row("dotmatch_hamming_k1", 100000, 1, reads_per_sec="200.0"),
        _row("guide_counter_one_mismatch", 100000, 1, reads_per_sec="100.0"),
        _row("dotmatch_exact_k0", 100000, 1, reads_per_sec="300.0"),
        _row("mageck_count_exact", 100000, 1, reads_per_sec="30.0"),
        _row("guide_counter_exact", 100000, 1, reads_per_sec="150.0"),
    ]

    gate.report_gate(report, rows, require_guide_counter=True, failures=failures)

    assert failures == []
