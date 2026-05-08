import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "scripts" / "check_feature_barcode_gate.py"


def _load_gate():
    spec = importlib.util.spec_from_file_location("check_feature_barcode_gate", GATE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _row(**overrides) -> dict[str, str]:
    row = {
        "tool": "dotmatch_count",
        "workflow": "synthetic_feature_barcode_fixture",
        "status": "smoke",
        "exit_code": "0",
        "n_reads": "6",
        "n_targets": "4",
        "assigned_unique": "4",
        "ambiguous_reads": "1",
        "unmatched_reads": "1",
        "validation_mismatches": "0",
        "command": "dotmatch count --targets feature_barcodes.tsv --reads feature_reads.fastq",
    }
    row.update(overrides)
    return row


def _public_row(**overrides) -> dict[str, str]:
    row = {
        "tool": "dotmatch_count",
        "workflow": "public_10x_totalseq_b_feature_barcode",
        "status": "supported",
        "exit_code": "0",
        "n_reads": "20000",
        "n_targets": "10",
        "target_start": "10",
        "target_length": "15",
        "k": "0",
        "metric": "hamming",
        "assigned_unique": "18000",
        "assigned_exact": "18000",
        "corrected_reads": "0",
        "ambiguous_reads": "0",
        "unmatched_reads": "2000",
        "validation_mismatches": "0",
        "command": "dotmatch count --targets feature_barcodes.tsv --reads reads.fastq.gz --target-start 10 --target-length 15",
        "metadata": "examples/feature_barcode/data/metadata.json",
    }
    row.update(overrides)
    return row


def _write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def test_feature_barcode_gate_accepts_smoke_row_with_validation_and_diagnostics():
    gate = _load_gate()
    failures = []

    gate.row_gate([_row()], failures)

    assert failures == []


def test_feature_barcode_gate_requires_validation_and_whitelist_coverage():
    gate = _load_gate()
    failures = []

    gate.row_gate([
        _row(validation_mismatches="1"),
        _row(assigned_unique="0"),
        _row(ambiguous_reads="0"),
    ], failures)

    assert any("zero validation mismatches" in failure for failure in failures)
    assert any("assign at least one feature barcode" in failure for failure in failures)
    assert any("ambiguous" in failure for failure in failures)


def test_feature_barcode_gate_rejects_missing_rows_from_cli(tmp_path, monkeypatch):
    gate = _load_gate()
    csv_path = tmp_path / "feature_barcode.csv"
    _write_rows(csv_path, [])
    monkeypatch.setattr(gate, "RAW", csv_path)

    result = gate.main([])

    assert result == 1


def test_feature_barcode_gate_accepts_public_rows_with_exact_baseline():
    gate = _load_gate()
    rows = [
        _row(),
        _public_row(k="0", assigned_unique="18000", assigned_exact="18000"),
        _public_row(k="1", assigned_unique="18100", assigned_exact="18000", corrected_reads="100"),
        _public_row(tool="exact_slice_hash", k="0", metric="exact", command="python3 scripts/bench_feature_barcode.py --metadata examples/feature_barcode/data/metadata.json", assigned_unique="18000", assigned_exact="18000"),
    ]
    failures = []

    gate.public_row_gate(rows, failures)

    assert failures == []


def test_feature_barcode_gate_rejects_public_rows_without_exact_agreement():
    gate = _load_gate()
    rows = [
        _public_row(k="0", assigned_unique="18000"),
        _public_row(k="1", assigned_unique="18100"),
        _public_row(tool="exact_slice_hash", k="0", metric="exact", assigned_unique="17999"),
    ]
    failures = []

    gate.public_row_gate(rows, failures)

    assert any("exact-slice baseline" in failure for failure in failures)
