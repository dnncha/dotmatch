import csv
import importlib.util
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "bench_bowtie1_crispr.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("bench_bowtie1_crispr", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_extracts_fixed_guide_window_and_reads_dotmatch_targets(tmp_path):
    bench = _load_module()
    reads = tmp_path / "reads.fastq"
    reads.write_text(
        "@r1\nNNACGTAA\n+\nIIIIIIII\n"
        "@r2\nTTACGAAA\n+\nJJJJJJJJ\n",
        encoding="utf-8",
    )
    extracted = tmp_path / "guide.fastq"

    n_reads = bench.extract_guide_fastq([reads], extracted, start=2, length=4)
    targets_csv = tmp_path / "guides.csv"
    targets_csv.write_text(
        "id,gRNA.sequence,Gene\n"
        "guide_a,ACGT,GENEA\n"
        "guide_b,ACGA,GENEB\n",
        encoding="utf-8",
    )

    assert n_reads == 2
    assert extracted.read_text(encoding="utf-8") == (
        "@r1\nACGT\n+\nIIII\n"
        "@r2\nACGA\n+\nJJJJ\n"
    )
    assert bench.read_guides(targets_csv) == [
        ("guide_a", "ACGT"),
        ("guide_b", "ACGA"),
    ]


def test_counts_unique_ambiguous_and_unmatched_bowtie_best_strata_rows(tmp_path):
    bench = _load_module()
    bowtie_out = tmp_path / "bowtie.out"
    bowtie_out.write_text(
        "r1\t+\tguide_a\t0\tACGT\tIIII\t0\t\n"
        "r2\t+\tguide_a\t0\tACGA\tIIII\t1\t3:A>G\n"
        "r2\t+\tguide_b\t0\tACGA\tIIII\t1\t\n"
        "r3\t+\tguide_c\t0\tTTTT\tIIII\t0\t\n"
        "r3\t-\tguide_c\t0\tAAAA\tIIII\t0\t\n",
        encoding="utf-8",
    )

    stats = bench.parse_bowtie_assignments(bowtie_out, n_reads=4)

    assert stats == {
        "assigned_reads": "2",
        "ambiguous_reads": "1",
        "rejected_reads": "1",
    }


def test_extract_prefixes_repeated_read_ids_across_multiple_fastqs(tmp_path):
    bench = _load_module()
    read_a = tmp_path / "a.fastq"
    read_b = tmp_path / "b.fastq"
    read_a.write_text("@read1\nNNACGT\n+\nIIIIII\n", encoding="utf-8")
    read_b.write_text("@read1\nTTACGA\n+\nJJJJJJ\n", encoding="utf-8")
    extracted = tmp_path / "combined.fastq"

    n_reads = bench.extract_guide_fastq([read_a, read_b], extracted, start=2, length=4)

    assert n_reads == 2
    assert extracted.read_text(encoding="utf-8") == (
        "@0:read1\nACGT\n+\nIIII\n"
        "@1:read1\nACGA\n+\nJJJJ\n"
    )


def test_allow_missing_writes_not_installed_row(tmp_path):
    guides = tmp_path / "guides.tsv"
    guides.write_text("sgRNAID\tSeq\tgene\ng1\tACGT\tGENE\n", encoding="utf-8")
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r1\nNNACGT\n+\nIIIIII\n", encoding="utf-8")
    out = tmp_path / "rows.csv"
    env = os.environ.copy()
    env["PATH"] = str(tmp_path / "empty-path")

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--guides",
            str(guides),
            "--reads",
            str(reads),
            "--target-start",
            "2",
            "--target-length",
            "4",
            "--k",
            "2",
            "--out",
            str(out),
            "--workflow",
            "unit_smoke",
            "--dataset-id",
            "toy",
            "--sample-id",
            "s1",
            "--allow-missing",
        ],
        check=False,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert result.returncode == 0, result.stderr
    rows = list(csv.DictReader(out.open()))
    assert len(rows) == 1
    row = rows[0]
    assert row["tool"] == "bowtie1_crispr_hamming_k2"
    assert row["version"] == "not_installed"
    assert row["semantics"] == "hamming_k2_no_indels_bowtie1_v"
    assert row["n_reads"] == "1"
    assert row["n_targets"] == "1"
    assert row["seconds"] == "0.000000"
    assert row["reads_per_sec"] == "0.0"
    assert row["exit_code"] == "127"
    assert row["assigned_reads"] == ""
    assert row["ambiguous_reads"] == ""
    assert row["rejected_reads"] == ""
    assert row["workflow"] == "unit_smoke"
    assert row["dataset_id"] == "toy"
    assert row["sample_id"] == "s1"
    assert "-v 2 --best --strata --norc -a" in row["command"]
