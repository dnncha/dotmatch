import gzip
import json

import pytest

import dotmatch
from dotmatch.calibration_io import _dict_rows, _read_targets
from dotmatch.count_io import read_count_table
from dotmatch.fastq_io import iter_fastq


def test_posterior_min_ignores_library_targets_of_another_length(tmp_path, monkeypatch):
    from dotmatch.cli import main

    monkeypatch.setenv("DOTMATCH_PYTHON_NO_DELEGATE", "1")
    targets = tmp_path / "targets.tsv"
    targets.write_text("target_id\ttarget_seq\nlong\tACGTACGT\nshort\tGGGGCCC\n", encoding="utf-8")
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGTACGT\n+\nIIIIIIII\n", encoding="utf-8")
    summary = tmp_path / "summary.json"
    main([
        "count", "--targets", str(targets), "--reads", str(reads), "--target-start", "0", "--target-length", "8",
        "--k", "1", "--posterior-min", "0.9", "--out", str(tmp_path / "counts.tsv"), "--summary", str(summary),
    ])
    assert json.loads(summary.read_text(encoding="utf-8"))["assigned_unique"] == 1


def test_python_count_output_is_one_sample(tmp_path):
    path = tmp_path / "counts.tsv"
    path.write_text(
        "target_id\ttarget_seq\tgene\tcount_exact\tcount_corrected_substitution\tcount_corrected_insertion\t"
        "count_corrected_deletion\tcount_corrected_other\tcount_total\tambiguous_nearby\n"
        "t0\tACGT\tG\t3\t1\t0\t0\t0\t4\t0\n",
        encoding="utf-8",
    )
    table = read_count_table(path)
    assert table.sample_names == ("count_total",)
    assert table.counts == ((4,),)


def test_truncated_gzip_fastq_and_counts_raise_value_error(tmp_path):
    fastq = tmp_path / "reads.fastq.gz"
    fastq.write_bytes(gzip.compress(b"@r\nACGT\n+\nIIII\n" * 200)[:40])
    with pytest.raises(ValueError, match="truncated gzip"):
        list(iter_fastq(fastq))
    counts = tmp_path / "counts.tsv.gz"
    payload = b"id\ts\n" + b"".join(f"t{i}\t{i}\n".encode() for i in range(2000))
    compressed = gzip.compress(payload)
    counts.write_bytes(compressed[: len(compressed) // 2])
    with pytest.raises(ValueError, match="truncated gzip"):
        read_count_table(counts)


def test_calibration_rows_must_be_complete(tmp_path):
    path = tmp_path / "reads.tsv"
    path.write_text("read_id\tobserved\tquality\nr1\tACGT\n", encoding="utf-8")
    with pytest.raises(ValueError, match="different number of fields"):
        _dict_rows(path)


def test_calibration_targets_honour_bom_and_named_columns(tmp_path):
    path = tmp_path / "targets.tsv"
    path.write_text("﻿target_id\ttarget_seq\tgene\nt0\tAAAA\tCCCC\n", encoding="utf-8")
    assert _read_targets(path) == ["AAAA"]
