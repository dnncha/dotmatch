import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "scripts" / "check_oligo_adapter_gate.py"


def _load_gate():
    spec = importlib.util.spec_from_file_location("check_oligo_adapter_gate", GATE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _row(**overrides) -> dict[str, str]:
    row = {
        "tool": "dotmatch_count",
        "workflow": "synthetic_oligo_adapter_fixture",
        "status": "smoke",
        "exit_code": "0",
        "n_reads": "6",
        "n_targets": "3",
        "target_start": "0",
        "target_length": "12",
        "k": "1",
        "metric": "hamming",
        "assigned_unique": "4",
        "assigned_exact": "3",
        "corrected_reads": "1",
        "ambiguous_reads": "1",
        "unmatched_reads": "1",
        "validation_mismatches": "0",
        "command": "dotmatch count --targets adapter_oligos.tsv --reads adapter_reads.fastq",
    }
    row.update(overrides)
    return row


def _public_row(**overrides) -> dict[str, str]:
    row = {
        "tool": "dotmatch_count",
        "workflow": "public_fast_adapter_truseq_r1",
        "status": "supported",
        "exit_code": "0",
        "n_reads": "10000",
        "n_targets": "10",
        "target_start": "229",
        "target_length": "20",
        "k": "0",
        "metric": "hamming",
        "assigned_unique": "76",
        "assigned_exact": "76",
        "corrected_reads": "0",
        "ambiguous_reads": "0",
        "unmatched_reads": "9924",
        "validation_mismatches": "0",
        "command": "dotmatch count --targets adapter_oligos.tsv --reads reads.fastq.gz --target-start 229 --target-length 20",
        "metadata": "examples/oligo_adapter/data/metadata.json",
    }
    row.update(overrides)
    return row


def _write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def test_oligo_adapter_gate_accepts_smoke_row_with_diagnostics():
    gate = _load_gate()
    failures = []

    gate.row_gate([_row()], failures)

    assert failures == []


def test_oligo_adapter_gate_requires_validation_correction_and_diagnostics():
    gate = _load_gate()
    failures = []

    gate.row_gate([
        _row(validation_mismatches="1"),
        _row(corrected_reads="0"),
        _row(ambiguous_reads="0"),
        _row(unmatched_reads="0"),
    ], failures)

    assert any("zero validation mismatches" in failure for failure in failures)
    assert any("corrected" in failure for failure in failures)
    assert any("ambiguous" in failure for failure in failures)
    assert any("unmatched" in failure for failure in failures)


def test_oligo_adapter_gate_rejects_missing_rows_from_cli(tmp_path, monkeypatch):
    gate = _load_gate()
    csv_path = tmp_path / "oligo_adapter.csv"
    _write_rows(csv_path, [])
    monkeypatch.setattr(gate, "RAW", csv_path)

    result = gate.main([])

    assert result == 1


def test_oligo_adapter_gate_accepts_public_rows_with_exact_baseline():
    gate = _load_gate()
    rows = [
        _row(),
        _public_row(k="0", assigned_unique="76", assigned_exact="76"),
        _public_row(k="1", assigned_unique="90", assigned_exact="76", corrected_reads="14"),
        _public_row(
            tool="exact_slice_hash",
            k="0",
            metric="exact",
            assigned_unique="76",
            assigned_exact="76",
            command="python3 scripts/bench_oligo_adapter.py --include-public --metadata examples/oligo_adapter/data/metadata.json",
        ),
    ]
    failures = []

    gate.public_row_gate(rows, failures)

    assert failures == []


def test_oligo_adapter_gate_rejects_public_rows_without_exact_agreement():
    gate = _load_gate()
    rows = [
        _public_row(k="0", assigned_unique="76"),
        _public_row(k="1", assigned_unique="90"),
        _public_row(tool="exact_slice_hash", k="0", metric="exact", assigned_unique="75"),
    ]
    failures = []

    gate.public_row_gate(rows, failures)

    assert any("exact-slice baseline" in failure for failure in failures)
