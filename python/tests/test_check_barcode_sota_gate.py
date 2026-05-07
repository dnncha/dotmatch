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
