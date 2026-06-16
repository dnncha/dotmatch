import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "scripts" / "check_barcode_comparison_gate.py"


def _load_gate():
    spec = importlib.util.spec_from_file_location("check_barcode_comparison_gate", GATE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _metadata(path: Path, barcode_length: int, barcode_lengths: list[int], barcode_length_mode: str = "fixed") -> None:
    path.write_text(
        json.dumps(
            {
                "evidence_ready": True,
                "barcode_count": 192,
                "barcode_length": barcode_length,
                "barcode_length_mode": barcode_length_mode,
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


def test_barcode_metadata_rejects_variable_length_sheet_without_length_mode(tmp_path):
    gate = _load_gate()
    metadata = tmp_path / "metadata.json"
    _metadata(metadata, barcode_length=0, barcode_lengths=[4, 5, 6, 7, 8], barcode_length_mode="")
    failures = []

    gate.metadata_gate(metadata, failures)

    assert any("barcode length mode" in failure for failure in failures)


def test_barcode_metadata_accepts_auto_length_sheet(tmp_path):
    gate = _load_gate()
    metadata = tmp_path / "metadata.json"
    _metadata(metadata, barcode_length=0, barcode_lengths=[4, 5, 6, 7, 8], barcode_length_mode="auto")
    failures = []

    gate.metadata_gate(metadata, failures)

    assert not any("barcode length" in failure for failure in failures)


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
        "barcode_length": "8",
        "metric": "hamming",
        "assigned_reads": "80",
        "corrected_reads": "0",
        "ambiguous_reads": "0",
        "reads_per_sec": "1000.0",
        "assignment_engine": "hamming_k1_lookup_direct" if k == "1" else "hamming_exact_lookup_direct",
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


def test_hamming_radius_splitter_counts_as_second_comparator_for_k1_lane():
    gate = _load_gate()
    failures = []
    dotmatch = _row("dotmatch_demux", k="1")
    dotmatch["corrected_reads"] = "25"
    cutadapt = _row("cutadapt_demux", k="1")
    cutadapt["reads_per_sec"] = "100.0"
    oracle = _row("hamming_radius_splitter", k="1")
    oracle["corrected_reads"] = "25"
    oracle["reads_per_sec"] = "80.0"

    gate.row_gate([
        dotmatch,
        cutadapt,
        oracle,
    ], min_repeats=1, require_cutadapt=True, require_second_comparator=True,
        allow_fixture=False, failures=failures)

    assert not any("second successful comparator" in failure for failure in failures)
    assert not any("corrected_reads" in failure for failure in failures)


def test_hamming_radius_splitter_must_match_dotmatch_counts():
    gate = _load_gate()
    failures = []
    dotmatch = _row("dotmatch_demux", k="1")
    dotmatch["corrected_reads"] = "25"
    oracle = _row("hamming_radius_splitter", k="1")
    oracle["corrected_reads"] = "24"

    gate.row_gate([
        dotmatch,
        _row("cutadapt_demux", k="1"),
        oracle,
    ], min_repeats=1, require_cutadapt=True, require_second_comparator=True,
        allow_fixture=False, failures=failures)

    assert any("corrected_reads" in failure for failure in failures)


def test_hash_splitter_counts_as_second_comparator_for_k0_exact_lane():
    gate = _load_gate()
    failures = []
    cutadapt = _row("cutadapt_demux", k="0")
    cutadapt["reads_per_sec"] = "100.0"
    hash_splitter = _row("hash_splitter_exact", k="0")
    hash_splitter["reads_per_sec"] = "300.0"

    gate.row_gate([
        _row("dotmatch_demux", k="0"),
        cutadapt,
        hash_splitter,
    ], min_repeats=1, require_cutadapt=True, require_second_comparator=True,
        allow_fixture=False, failures=failures)

    assert not any("second successful comparator" in failure for failure in failures)


def test_real_barcode_rows_must_assign_reads():
    gate = _load_gate()
    failures = []
    dotmatch = _row("dotmatch_demux", k="0")
    dotmatch["assigned_reads"] = "0"

    gate.row_gate([
        dotmatch,
        _row("cutadapt_demux", k="0"),
        _row("hash_splitter_exact", k="0"),
    ], min_repeats=1, require_cutadapt=True, require_second_comparator=True,
        allow_fixture=False, failures=failures)

    assert any("assigned zero reads" in failure for failure in failures)


def test_fixed_hamming_dotmatch_rows_must_record_direct_engine():
    gate = _load_gate()
    failures = []
    dotmatch = _row("dotmatch_demux", k="1")
    dotmatch["assignment_engine"] = "generic_indexed"
    oracle = _row("hamming_radius_splitter", k="1")

    gate.row_gate([
        dotmatch,
        _row("cutadapt_demux", k="1"),
        oracle,
    ], min_repeats=1, require_cutadapt=True, require_second_comparator=True,
        allow_fixture=False, failures=failures)

    assert any("assignment_engine" in failure for failure in failures)


def test_levenshtein_splitter_must_match_dotmatch_fixture_counts():
    gate = _load_gate()
    failures = []
    dotmatch = _row("dotmatch_demux", k="1")
    dotmatch.update({
        "workflow": "synthetic_levenshtein_one_edit_fixture",
        "metric": "levenshtein",
        "assignment_engine": "levenshtein_k1_lookup_direct",
        "assigned_reads": "80",
        "corrected_reads": "60",
        "ambiguous_reads": "0",
    })
    oracle = _row("levenshtein_radius_splitter", k="1")
    oracle.update({
        "workflow": "synthetic_levenshtein_one_edit_fixture",
        "metric": "levenshtein",
        "assigned_reads": "80",
        "corrected_reads": "60",
        "ambiguous_reads": "0",
    })

    gate.row_gate([dotmatch, oracle], min_repeats=1, require_cutadapt=False,
                  require_second_comparator=False, allow_fixture=True, failures=failures)

    assert not any("Levenshtein splitter" in failure for failure in failures)
    assert not any("corrected_reads" in failure for failure in failures)


def test_levenshtein_fixture_rows_must_record_direct_engine():
    gate = _load_gate()
    failures = []
    dotmatch = _row("dotmatch_demux", k="1")
    dotmatch.update({
        "workflow": "synthetic_levenshtein_one_edit_fixture",
        "metric": "levenshtein",
        "assignment_engine": "generic_indexed",
        "corrected_reads": "60",
    })
    oracle = _row("levenshtein_radius_splitter", k="1")
    oracle.update({
        "workflow": "synthetic_levenshtein_one_edit_fixture",
        "metric": "levenshtein",
        "corrected_reads": "60",
    })

    gate.row_gate([dotmatch, oracle], min_repeats=1, require_cutadapt=False,
                  require_second_comparator=False, allow_fixture=True, failures=failures)

    assert any("Levenshtein" in failure and "assignment_engine" in failure for failure in failures)


def test_real_barcode_rows_must_clear_speedup_floor():
    gate = _load_gate()
    failures = []
    dotmatch = _row("dotmatch_demux", k="0")
    dotmatch["reads_per_sec"] = "1000.0"
    cutadapt = _row("cutadapt_demux", k="0")
    cutadapt["reads_per_sec"] = "900.0"
    hash_splitter = _row("hash_splitter_exact", k="0")
    hash_splitter["reads_per_sec"] = "800.0"

    gate.row_gate([
        dotmatch,
        cutadapt,
        hash_splitter,
    ], min_repeats=1, require_cutadapt=True, require_second_comparator=True,
        allow_fixture=False, failures=failures)

    assert any("speedup below" in failure for failure in failures)


def test_barcode_report_gate_requires_published_gated_speedup_rows(tmp_path):
    gate = _load_gate()
    failures = []
    dotmatch = _row("dotmatch_demux", k="0")
    dotmatch.update({"workflow": "real_srp009896_inline_barcode", "reads_per_sec": "3000.0"})
    cutadapt = _row("cutadapt_demux", k="0")
    cutadapt.update({"workflow": "real_srp009896_inline_barcode", "reads_per_sec": "500.0"})
    hash_splitter = _row("hash_splitter_exact", k="0")
    hash_splitter.update({"workflow": "real_srp009896_inline_barcode", "reads_per_sec": "1000.0"})
    speedups = gate.gated_speedups([dotmatch, cutadapt, hash_splitter], failures)
    report = tmp_path / "README.md"
    report.write_text(
        "# Barcode\n\n"
        "## Gated Real-Data Speedups\n\n"
        "| workflow | comparator | DotMatch reads/sec | comparator reads/sec | speedup | gate floor |\n"
        "| --- | --- | ---: | ---: | ---: | ---: |\n"
        "| real_srp009896_inline_barcode | cutadapt_demux | 3000.0 | 500.0 | 6.00x | 5.00x |\n"
        "not public real-data evidence\n"
        "Levenshtein indel lane is fixture evidence\n",
        encoding="utf-8",
    )

    gate.report_gate(report, speedups, failures)

    assert any("hash_splitter_exact" in failure for failure in failures)


def test_barcode_report_gate_accepts_current_speedup_rows(tmp_path):
    gate = _load_gate()
    failures = []
    dotmatch = _row("dotmatch_demux", k="1")
    dotmatch.update({
        "workflow": "real_srp009896_inline_barcode_fixed8_k1",
        "reads_per_sec": "5000.0",
        "corrected_reads": "25",
    })
    cutadapt = _row("cutadapt_demux", k="1")
    cutadapt.update({"workflow": "real_srp009896_inline_barcode_fixed8_k1", "reads_per_sec": "1000.0"})
    oracle = _row("hamming_radius_splitter", k="1")
    oracle.update({
        "workflow": "real_srp009896_inline_barcode_fixed8_k1",
        "reads_per_sec": "400.0",
        "corrected_reads": "25",
    })
    speedups = gate.gated_speedups([dotmatch, cutadapt, oracle], failures)
    report = tmp_path / "README.md"
    report.write_text(
        "# Barcode\n\n"
        "## Gated Real-Data Speedups\n\n"
        "| workflow | comparator | DotMatch reads/sec | comparator reads/sec | speedup | gate floor |\n"
        "| --- | --- | ---: | ---: | ---: | ---: |\n"
        "| real_srp009896_inline_barcode_fixed8_k1 | cutadapt_demux | 5000.0 | 1000.0 | 5.00x | 5.00x |\n"
        "| real_srp009896_inline_barcode_fixed8_k1 | hamming_radius_splitter | 5000.0 | 400.0 | 12.50x | 12.00x |\n"
        "not public real-data evidence\n"
        "Levenshtein indel lane is fixture evidence\n",
        encoding="utf-8",
    )

    gate.report_gate(report, speedups, failures)

    assert failures == []
