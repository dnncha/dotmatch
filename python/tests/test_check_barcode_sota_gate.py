import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "scripts" / "check_barcode_sota_gate.py"


def _load_gate():
    spec = importlib.util.spec_from_file_location("check_barcode_sota_gate", GATE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _metadata(path: Path, barcode_length: int, barcode_lengths: list[int]) -> None:
    path.write_text(
        json.dumps(
            {
                "claim_grade_ready": True,
                "barcode_count": 192,
                "barcode_length": barcode_length,
                "barcode_lengths": barcode_lengths,
                "runs": [
                    {
                        "accession": "SRR391079",
                        "ena": {"fastq_md5": "remote-md5"},
                    }
                ],
            }
        ),
        encoding="utf-8",
    )


def test_barcode_metadata_rejects_variable_length_sheet_without_fixed_benchmark_length(tmp_path):
    gate = _load_gate()
    metadata = tmp_path / "metadata.json"
    _metadata(metadata, barcode_length=0, barcode_lengths=[4, 5, 6, 7, 8])
    failures = []

    gate.metadata_gate(metadata, failures)

    assert any("fixed barcode length" in failure for failure in failures)


def test_barcode_metadata_accepts_declared_fixed_benchmark_length(tmp_path):
    gate = _load_gate()
    metadata = tmp_path / "metadata.json"
    _metadata(metadata, barcode_length=8, barcode_lengths=[8])
    failures = []

    gate.metadata_gate(metadata, failures)

    assert not any("fixed barcode length" in failure for failure in failures)


def _row(tool: str, k: str = "0") -> dict[str, str]:
    return {
        "tool": tool,
        "workflow": "real_public_inline_barcode",
        "exit_code": "0",
        "n_reads": "100",
        "n_barcodes": "4",
        "k": k,
    }


def test_hash_splitter_counts_as_second_comparator_for_exact_lane_only():
    gate = _load_gate()
    failures = []

    gate.row_gate([
        _row("dotmatch_demux", k="1"),
        _row("cutadapt_demux", k="1"),
        _row("hash_splitter_exact", k="1"),
    ], min_repeats=1, require_cutadapt=True, require_second_comparator=True,
        allow_fixture=False, failures=failures)

    assert any("second successful comparator" in failure for failure in failures)


def test_hash_splitter_counts_as_second_comparator_for_k0_exact_lane():
    gate = _load_gate()
    failures = []

    gate.row_gate([
        _row("dotmatch_demux", k="0"),
        _row("cutadapt_demux", k="0"),
        _row("hash_splitter_exact", k="0"),
    ], min_repeats=1, require_cutadapt=True, require_second_comparator=True,
        allow_fixture=False, failures=failures)

    assert not any("second successful comparator" in failure for failure in failures)
