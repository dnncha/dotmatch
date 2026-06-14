import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "scripts" / "check_bcl_tiny_public_gate.py"


def _load_gate():
    spec = importlib.util.spec_from_file_location("check_bcl_tiny_public_gate", GATE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _dotmatch_row(**overrides) -> dict[str, str]:
    row = {
        "tool": "dotmatch_bcl_demux",
        "workflow": "public_10x_tiny_bcl",
        "format": "classic_bcl",
        "clusters": "2136539",
        "cycles": "132",
        "samples": "1",
        "total_clusters": "2136539",
        "assigned_reads": "2079501",
        "undetermined_reads": "57038",
        "filtered_clusters": "0",
        "tiles": "1",
        "exit_code": "0",
        "output_sha256": "a" * 64,
        "fastq_content_sha256": "b" * 64,
        "command": "dotmatch bcl-demux --run-folder examples/bcl_demux/data/cellranger-tiny-bcl-1.2.0 --sample-sheet examples/bcl_demux/data/cellranger-tiny-bcl-samplesheet.normalized.csv",
    }
    row.update(overrides)
    return row


def _bcl2fastq_row(**overrides) -> dict[str, str]:
    row = {
        "tool": "bcl2fastq",
        "workflow": "public_10x_tiny_bcl",
        "format": "classic_bcl",
        "clusters": "2136539",
        "total_clusters": "2136539",
        "assigned_reads": "2079501",
        "undetermined_reads": "57038",
        "exit_code": "0",
        "validation_mismatches": "0",
        "validation_exit_code": "0",
        "validation_mode": "count_totals",
        "fastq_content_sha256": "c" * 64,
        "command": "bcl2fastq --runfolder-dir examples/bcl_demux/data/cellranger-tiny-bcl-1.2.0",
    }
    row.update(overrides)
    return row


def _bcl2fastq_unavailable_row(**overrides) -> dict[str, str]:
    row = {
        "tool": "bcl2fastq",
        "version": "not_installed",
        "workflow": "public_10x_tiny_bcl",
        "format": "classic_bcl_or_cbcl",
        "clusters": "2136539",
        "total_clusters": "2136539",
        "exit_code": "127",
        "command": "bcl2fastq not-run",
    }
    row.update(overrides)
    return row


def _write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def test_bcl_tiny_public_gate_accepts_public_dotmatch_and_validated_bcl2fastq_rows(tmp_path):
    gate = _load_gate()
    csv_path = tmp_path / "bcl_demux.csv"
    report = tmp_path / "README.md"
    _write_rows(csv_path, [_dotmatch_row(), _bcl2fastq_row()])
    report.write_text(
        "# BCL benchmark report\n\n"
        "Current status: DotMatch supports a first lane-1 classic per-cycle BCL milestone. "
        "CBCL/NovaSeq-style input and broader multi-lane BCL operation are still gated as future work.\n"
        "Rows with exit code `127` are environment records for missing competitors, not runtime comparisons.\n"
        "Run `make bcl-tiny-public-gate` to verify the narrow public 10x tiny-BCL classic per-cycle milestone. "
        "This gate checks bcl2fastq count-total validation when bcl2fastq is available; "
        "an explicit `not_installed` environment row keeps the boundary visible without turning it into a false comparison.\n"
        "Broader raw-BCL evaluation needs real classic-BCL and CBCL run-folder rows, a successful DotMatch CBCL row, "
        "and `dotmatch bcl-validate` zero-mismatch evidence.\n",
        encoding="utf-8",
    )

    failures = []
    gate.row_gate(gate.read_rows(csv_path), failures)
    gate.report_gate(report, failures)

    assert failures == []


def test_bcl_tiny_public_gate_accepts_explicitly_unavailable_bcl2fastq(tmp_path):
    gate = _load_gate()
    csv_path = tmp_path / "bcl_demux.csv"
    _write_rows(csv_path, [_dotmatch_row(), _bcl2fastq_unavailable_row()])

    failures = []
    gate.row_gate(gate.read_rows(csv_path), failures)

    assert failures == []


def test_bcl_tiny_public_gate_rejects_installed_but_invalid_comparator(tmp_path):
    gate = _load_gate()
    csv_path = tmp_path / "bcl_demux.csv"
    _write_rows(csv_path, [_dotmatch_row(), _bcl2fastq_row(validation_mismatches="2")])

    failures = []
    gate.row_gate(gate.read_rows(csv_path), failures)

    assert any("validated bcl2fastq" in failure for failure in failures)


def test_bcl_tiny_public_gate_rejects_count_disagreement(tmp_path):
    gate = _load_gate()
    csv_path = tmp_path / "bcl_demux.csv"
    _write_rows(csv_path, [_dotmatch_row(), _bcl2fastq_row(assigned_reads="1")])

    failures = []
    gate.row_gate(gate.read_rows(csv_path), failures)

    assert any("assigned read count" in failure for failure in failures)


def test_bcl_tiny_public_gate_rejects_empty_report(tmp_path):
    gate = _load_gate()
    report = tmp_path / "README.md"
    report.write_text("", encoding="utf-8")

    failures = []
    gate.report_gate(report, failures)

    assert any("empty BCL benchmark report" in failure for failure in failures)


def test_bcl_tiny_public_gate_rejects_report_without_boundary(tmp_path):
    gate = _load_gate()
    report = tmp_path / "README.md"
    report.write_text("# BCL benchmark report\n\nFast BCL demux.\n", encoding="utf-8")

    failures = []
    gate.report_gate(report, failures)

    assert any("classic per-cycle BCL milestone" in failure for failure in failures)
    assert any("not runtime comparisons" in failure for failure in failures)
    assert any("DotMatch CBCL row" in failure for failure in failures)
