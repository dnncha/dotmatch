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
