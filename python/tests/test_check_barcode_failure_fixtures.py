import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "check_barcode_failure_fixtures.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_barcode_failure_fixtures", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_barcode_failure_fixture_catalog_is_complete():
    checker = _load_checker()

    result = checker.audit(ROOT)

    assert result.ok, result.failures
    assert any("failure modes" in item for item in result.passed)


def test_barcode_failure_fixture_rejects_missing_required_mode(tmp_path):
    checker = _load_checker()
    root = tmp_path / "repo"
    fixture_dir = root / "examples" / "barcode_autopsy" / "failure_modes"
    fixture_dir.mkdir(parents=True)
    (fixture_dir / "barcodes.tsv").write_text("barcode_id\tbarcode_seq\ns1\tACGT\n", encoding="utf-8")
    (fixture_dir / "reads.fastq").write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")
    with (fixture_dir / "expected_findings.tsv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["failure_mode", "read_id", "observed", "diagnosis", "next_action"], delimiter="\t")
        writer.writeheader()
        writer.writerow({
            "failure_mode": "wrong_offset",
            "read_id": "wrong_offset",
            "observed": "NACG",
            "diagnosis": "Barcode window starts one base early.",
            "next_action": "Scan offsets and rerun with the best-supported start.",
        })

    result = checker.audit(root)

    assert not result.ok
    assert any("missing required failure mode" in failure for failure in result.failures)
