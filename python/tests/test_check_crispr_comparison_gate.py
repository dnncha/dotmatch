import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "scripts" / "check_crispr_comparison_gate.py"


def _load_gate():
    spec = importlib.util.spec_from_file_location("check_crispr_comparison_gate", GATE)
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
        "total_delta": "0",
        "differing_guides": "0",
        "pearson": "1.00000000",
        "spearman": "1.00000000",
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

    assert any("sanson_brunello Hamming count agreement is below evidence threshold" in f for f in failures)


def test_agreement_gate_rejects_ok_exact_rows_with_deltas_or_nan():
    gate = _load_gate()
    rows = [
        _agreement_row("mageck_yusa", "dotmatch_exact_vs_mageck_exact"),
        _agreement_row("mageck_yusa", "dotmatch_hamming_vs_guide_counter"),
        _agreement_row(
            "sanson_brunello",
            "dotmatch_exact_vs_mageck_exact",
            total_left="321536",
            total_right="0",
        ) | {"total_delta": "321536", "differing_guides": "67253", "pearson": "nan", "spearman": "nan"},
        _agreement_row("sanson_brunello", "dotmatch_hamming_vs_guide_counter"),
    ]
    failures = []

    gate.agreement_gate(rows, require_guide_counter=False, failures=failures)

    assert any("sanson_brunello exact total differs" in f for f in failures)
    assert any("sanson_brunello exact guide-level differences" in f for f in failures)
    assert any("finite Pearson" in f for f in failures)


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


def test_repeated_gate_accepts_guide_counter_rows_without_speed_superiority_claim():
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

    assert not any("guide-counter" in f for f in failures)


def test_repeated_gate_still_requires_guide_counter_when_requested():
    gate = _load_gate()
    rows = []
    for dataset in gate.DATASETS:
        rows.extend([
            _repeated_row(dataset, "dotmatch_exact_k0"),
            _speed_row(dataset, "dotmatch_hamming_k1", "99.0"),
            _repeated_row(dataset, "dotmatch_levenshtein_k1", verified_per_read="1.0"),
        ])
    failures = []

    gate.repeated_gate(rows, min_records=100000, min_repeats=1, require_full=False,
                       require_mageck=False, require_guide_counter=True, failures=failures)

    assert any("guide_counter_one_mismatch needs >= 1 repeats" in f for f in failures)


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


def test_repeated_gate_accepts_full_hamming_slower_than_full_guide_counter_when_full_rows_exist():
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

    assert not any("guide-counter" in f for f in failures)


def test_repeated_gate_requires_named_full_guide_counter_dataset_only():
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
        _full_row("sanson_brunello", "dotmatch_hamming_k1", gate.FULL_FASTQ_MIN_READS["sanson_brunello"]),
        _full_row("sanson_brunello", "guide_counter_one_mismatch", gate.FULL_FASTQ_MIN_READS["sanson_brunello"]),
    ])
    for row in rows:
        if row["requested_records_per_sample"] == "full":
            row["reads_per_sec"] = "100.0"
            row["seconds"] = "1.0"
    failures = []

    gate.repeated_gate(
        rows,
        min_records=1,
        min_repeats=1,
        require_full=False,
        require_mageck=False,
        require_guide_counter=True,
        failures=failures,
        required_full_guide_counter_datasets=["sanson_brunello"],
    )

    assert not failures


def test_repeated_gate_rejects_missing_named_full_guide_counter_dataset():
    gate = _load_gate()
    rows = []
    for dataset in gate.DATASETS:
        rows.extend([
            _speed_row(dataset, "dotmatch_exact_k0", "100.0"),
            _speed_row(dataset, "dotmatch_hamming_k1", "200.0"),
            _speed_row(dataset, "dotmatch_levenshtein_k1", "50.0"),
            _speed_row(dataset, "guide_counter_one_mismatch", "100.0"),
        ])
    failures = []

    gate.repeated_gate(
        rows,
        min_records=1,
        min_repeats=1,
        require_full=False,
        require_mageck=False,
        require_guide_counter=True,
        failures=failures,
        required_full_guide_counter_datasets=["sanson_brunello"],
    )

    assert any("sanson_brunello:dotmatch_hamming_k1 needs a full FASTQ guide-counter comparison row" in f for f in failures)
    assert any("sanson_brunello:guide_counter_one_mismatch needs a full FASTQ guide-counter comparison row" in f for f in failures)


def _validation_row(dataset, sample, **overrides):
    row = {
        "dataset": dataset,
        "sample": sample,
        "checked_reads": "1000",
        "mismatches": "0",
        "oracle_strategy": "bounded_edlib_candidates",
        "edlib_alignments": "12000",
        "bounded_windows": "3000",
        "fallback_windows": "50",
    }
    row.update({key: str(value) for key, value in overrides.items()})
    return row


def test_validation_gate_accepts_bounded_edlib_oracle_rows():
    gate = _load_gate()
    rows = [
        _validation_row("mageck_yusa", "plasmid"),
        _validation_row("sanson_brunello", "RepA"),
    ]
    failures = []

    gate.validation_gate(rows, min_checked=1000, failures=failures)

    assert failures == []


def test_validation_gate_rejects_weak_oracle_metadata():
    gate = _load_gate()
    rows = [
        _validation_row(
            "mageck_yusa",
            "plasmid",
            oracle_strategy="full_edlib_scan",
            edlib_alignments="0",
            bounded_windows="0",
            fallback_windows="51",
        ),
        _validation_row("sanson_brunello", "RepA"),
    ]
    failures = []

    gate.validation_gate(rows, min_checked=1000, failures=failures)

    assert any("bounded_edlib_candidates" in f for f in failures)
    assert any("edlib_alignments missing" in f for f in failures)
    assert any("bounded_windows missing" in f for f in failures)
    assert any("fallback_windows exceeds 5%" in f for f in failures)


def test_hamming_k23_comparator_gate_rejects_missing_bowtie1_rows():
    gate = _load_gate()
    rows = [
        {
            "dataset": "mageck_yusa",
            "k": "2",
            "records_per_sample": "100000",
            "comparison": "dotmatch_hamming_k2_vs_internal_baseline",
            "status": "ok",
        }
    ]
    failures = []

    gate.hamming_k23_comparator_gate(
        rows,
        required_ks=["2", "3"],
        failures=failures,
        min_records=100000,
        datasets=["mageck_yusa"],
    )

    assert any("missing DotMatch-vs-Bowtie 1 Hamming k2 comparator row for mageck_yusa" in f for f in failures)
    assert any("missing DotMatch-vs-Bowtie 1 Hamming k3 comparator row for mageck_yusa" in f for f in failures)


def test_hamming_k23_comparator_gate_accepts_minimal_passing_rows():
    gate = _load_gate()
    rows = [
        {
            "dataset": "mageck_yusa",
            "k": "2",
            "records_per_sample": "100000",
            "comparison": "dotmatch_hamming_k2_vs_bowtie1",
            "dotmatch_tool": "dotmatch_hamming_k2",
            "bowtie1_tool": "bowtie1",
            "status": "ok",
            "exit_code": "0",
            "semantics": "Hamming k=2, no indels, same-strand fixed guide window",
            "n_targets": "77441",
            "dotmatch_reads_per_sec": "120000.0",
            "bowtie1_reads_per_sec": "10000.0",
            "speedup": "12.0",
            "dotmatch_assigned_reads": "5000",
            "bowtie1_assigned_reads": "5000",
        },
        {
            "dataset": "mageck_yusa",
            "k": "3",
            "records_per_sample": "100000",
            "comparison": "dotmatch_hamming_k3_vs_bowtie1",
            "dotmatch_tool": "dotmatch_hamming_k3",
            "bowtie1_tool": "bowtie1",
            "status": "ok",
            "exit_code": "0",
            "semantics": "Hamming k=3, no indels, same-strand fixed guide window",
            "n_targets": "77441",
            "dotmatch_reads_per_sec": "12000.0",
            "bowtie1_reads_per_sec": "5000.0",
            "speedup": "2.4",
            "dotmatch_assigned_reads": "17000",
            "bowtie1_assigned_reads": "17000",
        },
    ]
    failures = []

    gate.hamming_k23_comparator_gate(
        rows,
        required_ks=["2", "3"],
        failures=failures,
        min_records=100000,
        datasets=["mageck_yusa"],
    )

    assert failures == []


def test_hamming_k23_comparator_gate_rejects_weak_speed_or_disagreement():
    gate = _load_gate()
    rows = [
        {
            "dataset": "mageck_yusa",
            "k": "2",
            "records_per_sample": "100000",
            "comparison": "dotmatch_hamming_k2_vs_bowtie1",
            "dotmatch_tool": "dotmatch_hamming_k2",
            "bowtie1_tool": "bowtie1",
            "status": "ok",
            "exit_code": "0",
            "semantics": "Hamming k=2, no indels, same-strand fixed guide window",
            "n_targets": "1000",
            "dotmatch_reads_per_sec": "11000.0",
            "bowtie1_reads_per_sec": "10000.0",
            "speedup": "1.1",
            "dotmatch_assigned_reads": "41",
            "bowtie1_assigned_reads": "40",
        }
    ]
    failures = []

    gate.hamming_k23_comparator_gate(
        rows,
        required_ks=["2"],
        failures=failures,
        min_records=100000,
        datasets=["mageck_yusa"],
    )

    assert any("speedup below" in f for f in failures)
    assert any("target library is too small" in f for f in failures)
    assert any("assigned reads differ" in f for f in failures)


def test_crispr_report_gate_requires_hamming_k23_rows(tmp_path):
    gate = _load_gate()
    failures = []
    report = tmp_path / "README.md"
    report.write_text(
        "# CRISPR\n\n"
        "Broad comparisons require `make crispr-comparison-gate` to pass.\n\n"
        "## Hamming k2/k3 External Comparator Rows\n\n"
        "compare DotMatch directly with Bowtie 1\n\n"
        "Hamming k=2 must clear >=8x vs Bowtie 1; Hamming k=3 must clear >=2x vs Bowtie 1.\n\n"
        "|dataset|k|records_per_sample|dotmatch_tool|bowtie1_tool|dotmatch_reads_per_sec|bowtie1_reads_per_sec|speedup|status|semantics|artifact|\n"
        "|---|---|---|---|---|---|---|---|---|---|---|\n"
        "## Edlib Oracle Validation\n\n"
        "bounded_edlib_candidates\n",
        encoding="utf-8",
    )
    rows = [
        {
            "dataset": "mageck_yusa",
            "k": "2",
            "records_per_sample": "100000",
            "dotmatch_tool": "dotmatch_hamming_k2",
            "bowtie1_tool": "bowtie1",
            "dotmatch_reads_per_sec": "300.0",
            "bowtie1_reads_per_sec": "150.0",
            "speedup": "2.00",
            "status": "ok",
            "semantics": "Hamming k=2, no indels",
            "artifact": "benchmarks/raw/k2.csv",
        }
    ]

    gate.report_gate(report, rows, failures)

    assert any("Hamming k2 comparator row" in failure for failure in failures)


def test_crispr_report_gate_accepts_matching_hamming_k23_rows(tmp_path):
    gate = _load_gate()
    failures = []
    report = tmp_path / "README.md"
    report.write_text(
        "# CRISPR\n\n"
        "Broad comparisons require `make crispr-comparison-gate` to pass.\n\n"
        "## Hamming k2/k3 External Comparator Rows\n\n"
        "compare DotMatch directly with Bowtie 1\n\n"
        "Hamming k=2 must clear >=8x vs Bowtie 1; Hamming k=3 must clear >=2x vs Bowtie 1.\n\n"
        "|dataset|k|records_per_sample|dotmatch_tool|bowtie1_tool|dotmatch_reads_per_sec|bowtie1_reads_per_sec|speedup|status|semantics|artifact|\n"
        "|---|---|---|---|---|---|---|---|---|---|---|\n"
        "|mageck_yusa|2|100000|dotmatch_hamming_k2|bowtie1|300.0|150.0|2.00|ok|Hamming k=2, no indels|benchmarks/raw/k2.csv|\n\n"
        "## Edlib Oracle Validation\n\n"
        "bounded_edlib_candidates\n",
        encoding="utf-8",
    )
    rows = [
        {
            "dataset": "mageck_yusa",
            "k": "2",
            "records_per_sample": "100000",
            "dotmatch_tool": "dotmatch_hamming_k2",
            "bowtie1_tool": "bowtie1",
            "dotmatch_reads_per_sec": "300.0",
            "bowtie1_reads_per_sec": "150.0",
            "speedup": "2.00",
            "status": "ok",
            "semantics": "Hamming k=2, no indels",
            "artifact": "benchmarks/raw/k2.csv",
        }
    ]

    gate.report_gate(report, rows, failures)

    assert failures == []
