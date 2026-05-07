import csv
import importlib.util
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
VALIDATOR = ROOT / "scripts" / "validate_crispr_sota_edlib.py"


def _load_validator():
    spec = importlib.util.spec_from_file_location("validate_crispr_sota_edlib", VALIDATOR)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_validation_main_uses_persistent_work_dir(tmp_path, monkeypatch):
    validator = _load_validator()
    out = tmp_path / "raw" / "validation.csv"
    work = tmp_path / "work"
    library = tmp_path / "library.tsv"
    reads = tmp_path / "reads.fastq.gz"
    library.write_text("id\tseq\nsg1\tACGT\n", encoding="utf-8")
    reads.write_bytes(b"")

    def fake_run(*args, **kwargs):
        return subprocess.CompletedProcess(args=args[0], returncode=0)

    def fake_dataset(records):
        return (
            "mageck_yusa",
            library,
            0,
            4,
            0,
            records,
            "best",
            0.005,
            [("sample", reads)],
        )

    seen_work_dirs = []

    def fake_validate_one(dataset_id, sample, targets, source_reads, target_start, target_length,
                          auto_offset, auto_offset_sample, records, offset_mode,
                          offset_min_fraction, sample_size, tmp, edlib_threads=1):
        seen_work_dirs.append(tmp)
        tmp.mkdir(parents=True, exist_ok=True)
        assignments = tmp / f"{dataset_id}.{sample}.assignments.tsv"
        assignments.write_text("read_id\tstatus\n", encoding="utf-8")
        return {
            "dataset": dataset_id,
            "sample": sample,
            "oracle": "edlib_native",
            "checked_reads": str(sample_size),
            "mismatches": "0",
            "assignments_path": str(assignments),
        }

    monkeypatch.setattr(validator.subprocess, "run", fake_run)
    monkeypatch.setattr(validator, "yusa_dataset", fake_dataset)
    monkeypatch.setattr(validator, "validate_one", fake_validate_one)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "validate_crispr_sota_edlib",
            "--datasets",
            "mageck_yusa",
            "--records-per-sample",
            "3",
            "--sample",
            "2",
            "--out",
            str(out),
            "--work-dir",
            str(work),
        ],
    )

    validator.main()

    assert seen_work_dirs == [work]
    with out.open() as fh:
        row = next(csv.DictReader(fh))
    assert Path(row["assignments_path"]).exists()


def test_validation_main_accepts_parallel_jobs(tmp_path, monkeypatch):
    validator = _load_validator()
    out = tmp_path / "validation.csv"
    work = tmp_path / "work"
    library = tmp_path / "library.tsv"
    reads_a = tmp_path / "a.fastq.gz"
    reads_b = tmp_path / "b.fastq.gz"
    library.write_text("id\tseq\nsg1\tACGT\n", encoding="utf-8")
    reads_a.write_bytes(b"")
    reads_b.write_bytes(b"")

    def fake_run(*args, **kwargs):
        return subprocess.CompletedProcess(args=args[0], returncode=0)

    def fake_dataset(records):
        return (
            "mageck_yusa",
            library,
            0,
            4,
            0,
            records,
            "best",
            0.005,
            [("sample_a", reads_a), ("sample_b", reads_b)],
        )

    def fake_validate_one(dataset_id, sample, targets, source_reads, target_start, target_length,
                          auto_offset, auto_offset_sample, records, offset_mode,
                          offset_min_fraction, sample_size, tmp, edlib_threads=1):
        tmp.mkdir(parents=True, exist_ok=True)
        assignments = tmp / f"{dataset_id}.{sample}.assignments.tsv"
        assignments.write_text("read_id\tstatus\n", encoding="utf-8")
        return {
            "dataset": dataset_id,
            "sample": sample,
            "oracle": "edlib_native",
            "checked_reads": str(sample_size),
            "mismatches": "0",
            "assignments_path": str(assignments),
        }

    monkeypatch.setattr(validator.subprocess, "run", fake_run)
    monkeypatch.setattr(validator, "yusa_dataset", fake_dataset)
    monkeypatch.setattr(validator, "validate_one", fake_validate_one)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "validate_crispr_sota_edlib",
            "--datasets",
            "mageck_yusa",
            "--records-per-sample",
            "3",
            "--sample",
            "2",
            "--jobs",
            "2",
            "--out",
            str(out),
            "--work-dir",
            str(work),
        ],
    )

    validator.main()

    with out.open() as fh:
        rows = list(csv.DictReader(fh))
    assert [row["sample"] for row in rows] == ["sample_a", "sample_b"]


def test_validation_resume_skips_complete_rows(tmp_path, monkeypatch):
    validator = _load_validator()
    out = tmp_path / "validation.csv"
    work = tmp_path / "work"
    library = tmp_path / "library.tsv"
    reads_a = tmp_path / "a.fastq.gz"
    reads_b = tmp_path / "b.fastq.gz"
    library.write_text("id\tseq\nsg1\tACGT\n", encoding="utf-8")
    reads_a.write_bytes(b"")
    reads_b.write_bytes(b"")
    out.write_text(
        "dataset,sample,checked_reads,mismatches,records_available_for_validation,assignments_path\n"
        "mageck_yusa,sample_a,2,0,3,/tmp/sample_a.assignments.tsv\n"
        "mageck_yusa,sample_b,1,0,3,/tmp/sample_b.assignments.tsv\n",
        encoding="utf-8",
    )

    def fake_run(*args, **kwargs):
        return subprocess.CompletedProcess(args=args[0], returncode=0)

    def fake_dataset(records):
        return (
            "mageck_yusa",
            library,
            0,
            4,
            0,
            records,
            "best",
            0.005,
            [("sample_a", reads_a), ("sample_b", reads_b)],
        )

    called = []

    def fake_validate_one(dataset_id, sample, targets, source_reads, target_start, target_length,
                          auto_offset, auto_offset_sample, records, offset_mode,
                          offset_min_fraction, sample_size, tmp, edlib_threads=1):
        called.append(sample)
        return {
            "dataset": dataset_id,
            "sample": sample,
            "oracle": "edlib_native",
            "checked_reads": str(sample_size),
            "mismatches": "0",
            "records_available_for_validation": str(records),
            "assignments_path": str(tmp / f"{sample}.assignments.tsv"),
        }

    monkeypatch.setattr(validator.subprocess, "run", fake_run)
    monkeypatch.setattr(validator, "yusa_dataset", fake_dataset)
    monkeypatch.setattr(validator, "validate_one", fake_validate_one)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "validate_crispr_sota_edlib",
            "--datasets",
            "mageck_yusa",
            "--records-per-sample",
            "3",
            "--sample",
            "2",
            "--resume",
            "--out",
            str(out),
            "--work-dir",
            str(work),
        ],
    )

    validator.main()

    assert called == ["sample_b"]
    with out.open() as fh:
        rows = list(csv.DictReader(fh))
    assert [(row["sample"], row["checked_reads"]) for row in rows] == [("sample_a", "2"), ("sample_b", "2")]


def test_validation_resume_checkpoints_completed_rows_before_later_failure(tmp_path, monkeypatch):
    validator = _load_validator()
    out = tmp_path / "validation.csv"
    work = tmp_path / "work"
    library = tmp_path / "library.tsv"
    reads_a = tmp_path / "a.fastq.gz"
    reads_b = tmp_path / "b.fastq.gz"
    library.write_text("id\tseq\nsg1\tACGT\n", encoding="utf-8")
    reads_a.write_bytes(b"")
    reads_b.write_bytes(b"")

    def fake_run(*args, **kwargs):
        return subprocess.CompletedProcess(args=args[0], returncode=0)

    def fake_dataset(records):
        return (
            "mageck_yusa",
            library,
            0,
            4,
            0,
            records,
            "best",
            0.005,
            [("sample_a", reads_a), ("sample_b", reads_b)],
        )

    def fake_validate_one(dataset_id, sample, targets, source_reads, target_start, target_length,
                          auto_offset, auto_offset_sample, records, offset_mode,
                          offset_min_fraction, sample_size, tmp, edlib_threads=1):
        if sample == "sample_b":
            raise SystemExit("forced validation failure")
        return {
            "dataset": dataset_id,
            "sample": sample,
            "oracle": "edlib_native",
            "checked_reads": str(sample_size),
            "mismatches": "0",
            "records_available_for_validation": str(records),
            "assignments_path": str(tmp / f"{sample}.assignments.tsv"),
        }

    monkeypatch.setattr(validator.subprocess, "run", fake_run)
    monkeypatch.setattr(validator, "yusa_dataset", fake_dataset)
    monkeypatch.setattr(validator, "validate_one", fake_validate_one)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "validate_crispr_sota_edlib",
            "--datasets",
            "mageck_yusa",
            "--records-per-sample",
            "3",
            "--sample",
            "2",
            "--resume",
            "--out",
            str(out),
            "--work-dir",
            str(work),
        ],
    )

    try:
        validator.main()
    except SystemExit:
        pass
    else:
        raise AssertionError("expected validation failure")

    with out.open() as fh:
        rows = list(csv.DictReader(fh))
    assert [(row["sample"], row["checked_reads"]) for row in rows] == [("sample_a", "2")]


def test_validate_one_records_bounded_edlib_oracle_stats(tmp_path, monkeypatch):
    validator = _load_validator()
    targets = tmp_path / "library.tsv"
    source_reads = tmp_path / "source.fastq.gz"
    targets.write_text("id\tseq\nsg1\tACGT\n", encoding="utf-8")
    source_reads.write_bytes(b"")

    def fake_write_fastq_prefix(src, dest, records):
        dest.write_bytes(b"")

    def fake_assignment_strata(dataset_id, targets_path, reads_path, sample, target_start,
                               target_length, auto_offset, auto_offset_sample, offset_mode,
                               offset_min_fraction, records, tmp):
        assignments = tmp / f"{dataset_id}.{sample}.assignments.tsv"
        assignments.write_text("read_id\tstatus\n", encoding="utf-8")
        return assignments, 0, "[0]", {
            "exact": 1,
            "corrected": 0,
            "ambiguous": 0,
            "unmatched": 0,
            "contains_n": 0,
            "offset_shift_candidate": 0,
            "indel_window_candidate": 0,
        }

    def fake_run_json(cmd):
        return {
            "oracle": "edlib_native",
            "checked_reads": 3,
            "mismatches": 0,
            "k": 1,
            "oracle_strategy": "bounded_edlib_candidates",
            "edlib_alignments": 9,
            "bounded_windows": 12,
            "fallback_windows": 0,
        }

    monkeypatch.setattr(validator, "write_fastq_prefix", fake_write_fastq_prefix)
    monkeypatch.setattr(validator, "dotmatch_assignment_strata", fake_assignment_strata)
    monkeypatch.setattr(validator, "run_json", fake_run_json)

    row = validator.validate_one(
        "dataset",
        "sample",
        targets,
        source_reads,
        0,
        4,
        0,
        3,
        3,
        "best",
        0.005,
        3,
        tmp_path,
        edlib_threads=2,
    )

    assert row["oracle_strategy"] == "bounded_edlib_candidates"
    assert row["edlib_alignments"] == "9"
    assert row["bounded_windows"] == "12"
    assert row["fallback_windows"] == "0"
