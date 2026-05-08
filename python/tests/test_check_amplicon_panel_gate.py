import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "scripts" / "check_amplicon_panel_gate.py"


def _load_gate():
    spec = importlib.util.spec_from_file_location("check_amplicon_panel_gate", GATE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _row(**overrides) -> dict[str, str]:
    row = {
        "tool": "dotmatch_count",
        "workflow": "synthetic_amplicon_panel_fixture",
        "status": "smoke",
        "exit_code": "0",
        "n_reads": "6",
        "n_targets": "4",
        "assigned_unique": "4",
        "ambiguous_reads": "1",
        "unmatched_reads": "1",
        "validation_mismatches": "0",
        "command": "dotmatch count --targets panel_targets.tsv --reads panel_reads.fastq",
    }
    row.update(overrides)
    return row


def _public_row(**overrides) -> dict[str, str]:
    row = {
        "tool": "dotmatch_count",
        "workflow": "public_nfcore_artic_v3_amplicon_primer",
        "status": "supported",
        "exit_code": "0",
        "n_reads": "20000",
        "n_targets": "80",
        "target_start": "0",
        "target_length": "22",
        "k": "0",
        "metric": "hamming",
        "assigned_unique": "6400",
        "assigned_exact": "6400",
        "corrected_reads": "0",
        "ambiguous_reads": "0",
        "unmatched_reads": "13600",
        "validation_mismatches": "0",
        "command": "dotmatch count --targets primers.tsv --reads sample1_R1.fastq.gz --target-start 0 --target-length 22",
        "metadata": "examples/amplicon_panel/data/metadata.json",
    }
    row.update(overrides)
    return row


def _write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def test_amplicon_panel_gate_accepts_smoke_row_with_validation_and_diagnostics():
    gate = _load_gate()
    failures = []

    gate.row_gate([_row()], failures)

    assert failures == []


def test_amplicon_panel_gate_requires_validation_and_ambiguity_coverage():
    gate = _load_gate()
    failures = []

    gate.row_gate([
        _row(validation_mismatches="1"),
        _row(ambiguous_reads="0"),
    ], failures)

    assert any("zero validation mismatches" in failure for failure in failures)
    assert any("ambiguous" in failure for failure in failures)


def test_amplicon_panel_gate_rejects_missing_rows_from_cli(tmp_path, monkeypatch):
    gate = _load_gate()
    csv_path = tmp_path / "amplicon_panel.csv"
    _write_rows(csv_path, [])
    monkeypatch.setattr(gate, "RAW", csv_path)

    result = gate.main([])

    assert result == 1


def test_amplicon_panel_public_gate_accepts_exact_baseline_agreement():
    gate = _load_gate()
    rows = [
        _row(),
        _public_row(k="0", assigned_unique="6400", assigned_exact="6400"),
        _public_row(k="1", assigned_unique="6500", assigned_exact="6400", corrected_reads="100"),
        _public_row(
            tool="exact_prefix_hash",
            k="0",
            metric="exact",
            command="python3 scripts/bench_amplicon_panel.py --include-public --metadata examples/amplicon_panel/data/metadata.json",
            assigned_unique="6400",
            assigned_exact="6400",
        ),
    ]
    failures = []

    gate.public_row_gate(rows, failures)

    assert failures == []


def test_amplicon_panel_public_gate_rejects_exact_baseline_mismatch():
    gate = _load_gate()
    rows = [
        _public_row(k="0", assigned_unique="6400"),
        _public_row(k="1", assigned_unique="6500"),
        _public_row(tool="exact_prefix_hash", k="0", metric="exact", assigned_unique="6399"),
    ]
    failures = []

    gate.public_row_gate(rows, failures)

    assert any("exact-prefix baseline" in failure for failure in failures)
