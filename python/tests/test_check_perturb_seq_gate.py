import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "scripts" / "check_perturb_seq_gate.py"


def _load_gate():
    spec = importlib.util.spec_from_file_location("check_perturb_seq_gate", GATE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _row(**overrides) -> dict[str, str]:
    row = {
        "tool": "dotmatch_pair_count",
        "workflow": "synthetic_perturb_seq_fixture",
        "status": "smoke",
        "exit_code": "0",
        "n_reads": "7",
        "assigned_pairs": "3",
        "pair_ambiguous": "1",
        "left_unmatched": "1",
        "right_unmatched": "1",
        "invalid_reads": "1",
        "validation_mismatches": "0",
        "command": "dotmatch pair-count --left-targets guides.tsv --right-targets features.tsv",
    }
    row.update(overrides)
    return row


def _public_row(**overrides) -> dict[str, str]:
    row = {
        "tool": "dotmatch_count",
        "workflow": "public_10x_crispr_guide_capture",
        "status": "supported",
        "exit_code": "0",
        "n_reads": "20000",
        "n_targets": "1",
        "target_start": "63",
        "target_length": "19",
        "k": "0",
        "metric": "hamming",
        "assigned_unique": "15000",
        "assigned_exact": "15000",
        "corrected_reads": "0",
        "ambiguous_reads": "0",
        "unmatched_reads": "5000",
        "validation_mismatches": "0",
        "command": "dotmatch count --targets crispr_guides.tsv --reads reads.fastq.gz --target-start 63 --target-length 19",
        "metadata": "examples/perturb_seq/data/metadata.json",
    }
    row.update(overrides)
    return row


def _write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def test_perturb_seq_gate_accepts_smoke_row_with_pair_diagnostics():
    gate = _load_gate()
    failures = []

    gate.row_gate([_row()], failures)

    assert failures == []


def test_perturb_seq_gate_requires_validation_and_pair_diagnostics():
    gate = _load_gate()
    failures = []

    gate.row_gate([
        _row(validation_mismatches="1"),
        _row(assigned_pairs="0"),
        _row(pair_ambiguous="0"),
        _row(left_unmatched="0"),
        _row(right_unmatched="0"),
        _row(invalid_reads="0"),
    ], failures)

    assert any("zero validation mismatches" in failure for failure in failures)
    assert any("assign at least one guide-feature pair" in failure for failure in failures)
    assert any("ambiguous" in failure for failure in failures)
    assert any("left-unmatched" in failure for failure in failures)
    assert any("right-unmatched" in failure for failure in failures)
    assert any("invalid-read" in failure for failure in failures)


def test_perturb_seq_gate_rejects_missing_rows_from_cli(tmp_path, monkeypatch):
    gate = _load_gate()
    csv_path = tmp_path / "perturb_seq.csv"
    _write_rows(csv_path, [])
    monkeypatch.setattr(gate, "RAW", csv_path)

    result = gate.main([])

    assert result == 1


def test_perturb_seq_public_gate_accepts_exact_baseline_agreement():
    gate = _load_gate()
    rows = [
        _row(),
        _public_row(k="0", assigned_unique="15000", assigned_exact="15000"),
        _public_row(k="1", assigned_unique="15500", assigned_exact="15000", corrected_reads="500"),
        _public_row(
            tool="exact_slice_hash",
            k="0",
            metric="exact",
            command="python3 scripts/bench_perturb_seq.py --include-public --metadata examples/perturb_seq/data/metadata.json",
            assigned_unique="15000",
            assigned_exact="15000",
        ),
    ]
    failures = []

    gate.public_row_gate(rows, failures)

    assert failures == []


def test_perturb_seq_public_gate_rejects_exact_baseline_mismatch():
    gate = _load_gate()
    rows = [
        _public_row(k="0", assigned_unique="15000"),
        _public_row(k="1", assigned_unique="15500"),
        _public_row(tool="exact_slice_hash", k="0", metric="exact", assigned_unique="14999"),
    ]
    failures = []

    gate.public_row_gate(rows, failures)

    assert any("exact-slice baseline" in failure for failure in failures)
