import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "scripts" / "check_crispr_sota_gate.py"


def _load_gate():
    spec = importlib.util.spec_from_file_location("check_crispr_sota_gate", GATE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _agreement_row(dataset, comparison, status="ok", total_left="200000", total_right="200000"):
    return {
        "dataset": dataset,
        "comparison": f"{dataset}:{comparison}",
        "status": status,
        "total_left": total_left,
        "total_right": total_right,
    }


def test_strict_gate_rejects_shallow_guide_counter_count_agreement():
    gate = _load_gate()
    rows = [
        _agreement_row("mageck_yusa", "dotmatch_exact_vs_mageck_exact"),
        _agreement_row("mageck_yusa", "dotmatch_hamming_vs_guide_counter"),
        _agreement_row("sanson_brunello", "dotmatch_exact_vs_mageck_exact"),
        _agreement_row("sanson_brunello", "dotmatch_hamming_vs_guide_counter", total_left="16", total_right="17"),
    ]
    failures = []

    gate.agreement_gate(rows, require_guide_counter=True, failures=failures)

    assert any("sanson_brunello Hamming count agreement is below publication depth" in f for f in failures)


def _repeated_row(dataset, tool, verified_per_read=""):
    return {
        "dataset_id": dataset,
        "tool": tool,
        "requested_records_per_sample": "100000",
        "repeat": "1",
        "run_level": "subsample",
        "exit_code": "0",
        "n_targets": "77441",
        "verified_per_read": verified_per_read,
    }


def _full_row(dataset, tool, n_reads):
    row = _repeated_row(dataset, tool, verified_per_read="1.0")
    row["requested_records_per_sample"] = "full"
    row["run_level"] = "full"
    row["n_reads"] = str(n_reads)
    return row


def _full_sample_row(dataset, tool, sample_id, n_reads, seconds="1.0"):
    row = _full_row(dataset, tool, n_reads)
    row["run_level"] = "full_sample"
    row["sample_id"] = sample_id
    row["seconds"] = seconds
    row["reads_per_sec"] = str(float(n_reads) / float(seconds))
    return row


def test_repeated_gate_accepts_multi_offset_levenshtein_candidate_collapse():
    gate = _load_gate()
    rows = []
    for dataset in gate.DATASETS:
        rows.extend([
            _repeated_row(dataset, "dotmatch_exact_k0"),
            _repeated_row(dataset, "dotmatch_hamming_k1"),
            _repeated_row(dataset, "dotmatch_levenshtein_k1", verified_per_read="8.9290"),
        ])
    failures = []

    gate.repeated_gate(rows, min_records=100000, min_repeats=1, require_full=False,
                       require_mageck=False, require_guide_counter=False, failures=failures)

    assert not any("candidate collapse" in f for f in failures)


def test_repeated_gate_rejects_weak_levenshtein_candidate_collapse():
    gate = _load_gate()
    rows = []
    for dataset in gate.DATASETS:
        rows.extend([
            _repeated_row(dataset, "dotmatch_exact_k0"),
            _repeated_row(dataset, "dotmatch_hamming_k1"),
            _repeated_row(dataset, "dotmatch_levenshtein_k1", verified_per_read="100.0"),
        ])
    failures = []

    gate.repeated_gate(rows, min_records=100000, min_repeats=1, require_full=False,
                       require_mageck=False, require_guide_counter=False, failures=failures)

    assert any("candidate collapse" in f for f in failures)


def _speed_row(dataset, tool, reads_per_sec):
    row = _repeated_row(dataset, tool)
    row["reads_per_sec"] = str(reads_per_sec)
    return row


def test_repeated_gate_accepts_hamming_speedup_that_beats_guide_counter():
    gate = _load_gate()
    rows = []
    for dataset in gate.DATASETS:
        rows.extend([
            _repeated_row(dataset, "dotmatch_exact_k0"),
            _speed_row(dataset, "dotmatch_hamming_k1", "212.0"),
            _repeated_row(dataset, "dotmatch_levenshtein_k1", verified_per_read="1.0"),
            _speed_row(dataset, "guide_counter_one_mismatch", "100.0"),
        ])
    failures = []

    gate.repeated_gate(rows, min_records=100000, min_repeats=1, require_full=False,
                       require_mageck=False, require_guide_counter=True, failures=failures)

    assert not any("speedup vs guide-counter" in f for f in failures)


def test_repeated_gate_rejects_hamming_slower_than_guide_counter():
    gate = _load_gate()
    rows = []
    for dataset in gate.DATASETS:
        rows.extend([
            _repeated_row(dataset, "dotmatch_exact_k0"),
            _speed_row(dataset, "dotmatch_hamming_k1", "99.0"),
            _repeated_row(dataset, "dotmatch_levenshtein_k1", verified_per_read="1.0"),
            _speed_row(dataset, "guide_counter_one_mismatch", "100.0"),
        ])
    failures = []

    gate.repeated_gate(rows, min_records=100000, min_repeats=1, require_full=False,
                       require_mageck=False, require_guide_counter=True, failures=failures)

    assert any("speedup vs guide-counter" in f for f in failures)


def test_repeated_gate_no_full_speedup_ignores_full_rows():
    gate = _load_gate()
    rows = []
    for dataset in gate.DATASETS:
        rows.extend([
            _speed_row(dataset, "dotmatch_exact_k0", "100.0"),
            _speed_row(dataset, "dotmatch_hamming_k1", "200.0"),
            _speed_row(dataset, "dotmatch_levenshtein_k1", "50.0"),
            _speed_row(dataset, "guide_counter_one_mismatch", "100.0"),
            _full_sample_row(dataset, "dotmatch_hamming_k1", next(iter(gate.FULL_FASTQ_SAMPLE_READS[dataset])), 1, seconds="100.0"),
            _full_sample_row(dataset, "guide_counter_one_mismatch", next(iter(gate.FULL_FASTQ_SAMPLE_READS[dataset])), 1000, seconds="1.0"),
        ])
    failures = []

    gate.repeated_gate(rows, min_records=100000, min_repeats=1, require_full=False,
                       require_mageck=False, require_guide_counter=True, failures=failures)

    assert not any("DotMatch Hamming mean speedup vs guide-counter" in f for f in failures)


def test_repeated_gate_rejects_mislabeled_full_fastq_rows():
    gate = _load_gate()
    rows = []
    for dataset in gate.DATASETS:
        rows.extend([
            _repeated_row(dataset, "dotmatch_exact_k0"),
            _repeated_row(dataset, "dotmatch_hamming_k1"),
            _repeated_row(dataset, "dotmatch_levenshtein_k1", verified_per_read="1.0"),
        ])
    rows.extend([
        _full_row("mageck_yusa", "dotmatch_exact_k0", 200000),
        _full_row("mageck_yusa", "dotmatch_hamming_k1", 200000),
        _full_row("mageck_yusa", "dotmatch_levenshtein_k1", 200000),
        _full_row("sanson_brunello", "dotmatch_exact_k0", 246950411),
        _full_row("sanson_brunello", "dotmatch_hamming_k1", 246950411),
        _full_row("sanson_brunello", "dotmatch_levenshtein_k1", 246950411),
    ])
    failures = []

    gate.repeated_gate(rows, min_records=1, min_repeats=1, require_full=True,
                       require_mageck=False, require_guide_counter=False, failures=failures)

    assert any("mageck_yusa:dotmatch_exact_k0 full FASTQ row has too few reads" in f for f in failures)


def test_repeated_gate_accepts_full_fastq_rows_at_dataset_depth():
    gate = _load_gate()
    rows = []
    for dataset in gate.DATASETS:
        rows.extend([
            _repeated_row(dataset, "dotmatch_exact_k0"),
            _repeated_row(dataset, "dotmatch_hamming_k1"),
            _repeated_row(dataset, "dotmatch_levenshtein_k1", verified_per_read="1.0"),
        ])
    rows.extend([
        _full_row("mageck_yusa", "dotmatch_exact_k0", 20394663),
        _full_row("mageck_yusa", "dotmatch_hamming_k1", 20394663),
        _full_row("mageck_yusa", "dotmatch_levenshtein_k1", 20394663),
        _full_row("sanson_brunello", "dotmatch_exact_k0", 246950411),
        _full_row("sanson_brunello", "dotmatch_hamming_k1", 246950411),
        _full_row("sanson_brunello", "dotmatch_levenshtein_k1", 246950411),
    ])
    failures = []

    gate.repeated_gate(rows, min_records=1, min_repeats=1, require_full=True,
                       require_mageck=False, require_guide_counter=False, failures=failures)

    assert not any("full FASTQ row has too few reads" in f for f in failures)


def test_repeated_gate_accepts_complete_full_sample_rows_at_dataset_depth():
    gate = _load_gate()
    rows = []
    for dataset in gate.DATASETS:
        rows.extend([
            _speed_row(dataset, "dotmatch_exact_k0", "100.0"),
            _speed_row(dataset, "dotmatch_hamming_k1", "200.0"),
            _speed_row(dataset, "dotmatch_levenshtein_k1", "50.0"),
            _speed_row(dataset, "guide_counter_one_mismatch", "100.0"),
        ])
    for dataset, samples in gate.FULL_FASTQ_SAMPLE_READS.items():
        for tool in ["dotmatch_exact_k0", "dotmatch_hamming_k1", "dotmatch_levenshtein_k1", "guide_counter_one_mismatch"]:
            for sample_id, reads in samples.items():
                rows.append(_full_sample_row(dataset, tool, sample_id, reads, seconds="1.0"))
    failures = []

    gate.repeated_gate(rows, min_records=1, min_repeats=1, require_full=True,
                       require_mageck=False, require_guide_counter=True, failures=failures)

    assert not any("needs at least one full FASTQ timing row" in f for f in failures)
    assert not any("missing full rows for DotMatch-vs-guide-counter" in f for f in failures)


def test_repeated_gate_rejects_incomplete_full_sample_rows():
    gate = _load_gate()
    rows = []
    for dataset in gate.DATASETS:
        rows.extend([
            _repeated_row(dataset, "dotmatch_exact_k0"),
            _repeated_row(dataset, "dotmatch_hamming_k1"),
            _repeated_row(dataset, "dotmatch_levenshtein_k1", verified_per_read="1.0"),
        ])
    for sample_id, reads in gate.FULL_FASTQ_SAMPLE_READS["sanson_brunello"].items():
        if sample_id == "RepC":
            continue
        rows.append(_full_sample_row("sanson_brunello", "dotmatch_exact_k0", sample_id, reads))
    rows.extend([
        _full_row("mageck_yusa", "dotmatch_exact_k0", gate.FULL_FASTQ_MIN_READS["mageck_yusa"]),
        _full_row("mageck_yusa", "dotmatch_hamming_k1", gate.FULL_FASTQ_MIN_READS["mageck_yusa"]),
        _full_row("mageck_yusa", "dotmatch_levenshtein_k1", gate.FULL_FASTQ_MIN_READS["mageck_yusa"]),
    ])
    failures = []

    gate.repeated_gate(rows, min_records=1, min_repeats=1, require_full=True,
                       require_mageck=False, require_guide_counter=False, failures=failures)

    assert any("sanson_brunello:dotmatch_exact_k0 needs at least one full FASTQ timing row" in f for f in failures)


def test_repeated_gate_rejects_full_hamming_slower_than_full_guide_counter():
    gate = _load_gate()
    rows = []
    for dataset in gate.DATASETS:
        rows.extend([
            _speed_row(dataset, "dotmatch_exact_k0", "100.0"),
            _speed_row(dataset, "dotmatch_hamming_k1", "200.0"),
            _speed_row(dataset, "dotmatch_levenshtein_k1", "50.0"),
            _speed_row(dataset, "guide_counter_one_mismatch", "100.0"),
        ])
    rows.extend([
        _full_row("mageck_yusa", "dotmatch_exact_k0", 20394663),
        _full_row("mageck_yusa", "dotmatch_hamming_k1", 20394663),
        _full_row("mageck_yusa", "dotmatch_levenshtein_k1", 20394663),
        _full_row("mageck_yusa", "guide_counter_one_mismatch", 20394663),
        _full_row("sanson_brunello", "dotmatch_exact_k0", 246950411),
        _full_row("sanson_brunello", "dotmatch_hamming_k1", 246950411),
        _full_row("sanson_brunello", "dotmatch_levenshtein_k1", 246950411),
        _full_row("sanson_brunello", "guide_counter_one_mismatch", 246950411),
    ])
    for row in rows:
        if row["requested_records_per_sample"] == "full":
            row["reads_per_sec"] = "100.0"
    for row in rows:
        if (
            row["dataset_id"] == "mageck_yusa"
            and row["requested_records_per_sample"] == "full"
            and row["tool"] == "guide_counter_one_mismatch"
        ):
            row["reads_per_sec"] = "120.0"
    failures = []

    gate.repeated_gate(rows, min_records=1, min_repeats=1, require_full=True,
                       require_mageck=False, require_guide_counter=True, failures=failures)

    assert any("mageck_yusa full DotMatch Hamming speedup vs guide-counter" in f for f in failures)
