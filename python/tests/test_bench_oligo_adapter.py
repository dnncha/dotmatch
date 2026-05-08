import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "scripts" / "bench_oligo_adapter.py"


def _load_bench():
    spec = importlib.util.spec_from_file_location("bench_oligo_adapter", BENCH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_make_fixture_includes_adapter_oligos_and_diagnostics(tmp_path):
    bench = _load_bench()

    fixture = bench.make_fixture(tmp_path)

    targets = fixture.targets.read_text(encoding="utf-8")
    reads = fixture.reads.read_text(encoding="utf-8")
    assert "ADAPTER_A\tACGTACGTACGT\tadapter_A" in targets
    assert "ADAPTER_B\tTTTTCCCCAAAA\tadapter_B" in targets
    assert "ADAPTER_D\tACGTACGTACGA\tadapter_D" in targets
    assert "@exact_adapter_a\nACGTACGTACGT" in reads
    assert "@corrected_adapter_b\nTTTTCCCCAAAT" in reads
    assert "@ambiguous_adapter_ad\nACGTACGTACGC" in reads
    assert "@unmatched\nCCCCCCCCCCCC" in reads
    assert fixture.expected == {
        "total_reads": 6,
        "assigned_unique": 4,
        "assigned_exact": 3,
        "assigned_corrected": 1,
        "ambiguous": 1,
        "unmatched": 1,
    }


def test_validation_detects_oligo_adapter_summary_mismatches(tmp_path):
    bench = _load_bench()
    summary = tmp_path / "summary.json"
    summary.write_text(
        json.dumps(
            {
                "samples": [
                    {
                        "total_reads": 6,
                        "assigned_unique": 4,
                        "assigned_exact": 3,
                        "assigned_corrected": 1,
                        "ambiguous": 1,
                        "unmatched": 1,
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    assert bench.validation_mismatches(summary, {"assigned_unique": 4, "assigned_corrected": 1}) == []
    assert bench.validation_mismatches(summary, {"assigned_unique": 5}) == [
        "assigned_unique expected 5 observed 4"
    ]


def test_exact_slice_hash_counts_public_adapter_prefixes(tmp_path):
    bench = _load_bench()
    targets = tmp_path / "adapter_oligos.tsv"
    reads = tmp_path / "reads.fastq"
    targets.write_text(
        "target_id\ttarget_seq\tgene\n"
        "TruSeq_I7_full_length\tAGATCGGAAGAGCACACGTC\tTruSeq_I7_full_length\n"
        "TruSeq_I7_threeprime\tATCTCGTATGCCGTCTTCTG\tTruSeq_I7_threeprime\n",
        encoding="utf-8",
    )
    reads.write_text(
        "@r0\n" + ("A" * 229) + "AGATCGGAAGAGCACACGTC" + "GG\n+\n" + ("I" * 251) + "\n"
        "@r1\n" + ("C" * 229) + "ATCTCGTATGCCGTCTTCTG" + "GG\n+\n" + ("I" * 251) + "\n"
        "@r2\n" + ("G" * 229) + "NNNNNNNNNNNNNNNNNNNN" + "GG\n+\n" + ("I" * 251) + "\n",
        encoding="utf-8",
    )

    stats = bench.exact_slice_hash_stats(targets, reads, 229, 20)

    assert stats == {
        "n_reads": "3",
        "n_targets": "2",
        "assigned_unique": "2",
        "assigned_exact": "2",
        "corrected_reads": "0",
        "ambiguous_reads": "0",
        "unmatched_reads": "1",
    }


def test_public_benchmark_rows_require_metadata_and_include_exact_baseline(tmp_path, monkeypatch):
    bench = _load_bench()
    metadata = tmp_path / "metadata.json"
    targets = tmp_path / "adapter_oligos.tsv"
    reads = tmp_path / "reads.fastq"
    targets.write_text(
        "target_id\ttarget_seq\tgene\n"
        "TruSeq_I7_full_length\tAGATCGGAAGAGCACACGTC\tTruSeq_I7_full_length\n",
        encoding="utf-8",
    )
    reads.write_text(
        "@r0\n" + ("A" * 229) + "AGATCGGAAGAGCACACGTC" + "GG\n+\n" + ("I" * 251) + "\n",
        encoding="utf-8",
    )
    metadata.write_text(
        json.dumps(
            {
                "evidence_ready": True,
                "targets": str(targets),
                "local_fastq": str(reads),
                "target_start": 229,
                "target_length": 20,
                "target_count": 1,
                "written_records": 1,
                "dataset_id": "fast_adapter_trimming_truseq_r1",
            }
        ),
        encoding="utf-8",
    )

    def fake_run(cmd, cwd, stdout, stderr, text, check):
        assert "--target-start" in cmd
        assert "229" in cmd
        out_summary = tmp_path / (
            "public_fast_adapter_truseq_r1_k1_summary.json"
            if "--k" in cmd and cmd[cmd.index("--k") + 1] == "1"
            else "public_fast_adapter_truseq_r1_k0_summary.json"
        )
        out_summary.write_text(
            json.dumps(
                {
                    "n_targets": 1,
                    "alphabet_policy": "literal-byte; A/C/G/T/N/IUPAC symbols are ordinary byte symbols; no wildcard expansion",
                    "samples": [
                        {
                            "total_reads": 1,
                            "assigned_unique": 1,
                            "assigned_exact": 1,
                            "assigned_corrected": 0,
                            "ambiguous": 0,
                            "unmatched": 0,
                            "invalid": 0,
                            "candidates_verified": 1,
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )
        return type("Proc", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    monkeypatch.setattr(bench.subprocess, "run", fake_run)

    rows = bench.run_public_benchmark(Path("dotmatch"), metadata, tmp_path)

    assert [row["tool"] for row in rows] == ["dotmatch_count", "dotmatch_count", "exact_slice_hash"]
    assert {row["workflow"] for row in rows} == {"public_fast_adapter_truseq_r1"}
    assert rows[0]["k"] == "0"
    assert rows[1]["k"] == "1"
    assert rows[2]["metric"] == "exact"
    assert rows[2]["assigned_unique"] == "1"
    assert all(row["status"] == "supported" for row in rows)


def test_write_report_scopes_oligo_adapter_claims(tmp_path):
    bench = _load_bench()
    report = tmp_path / "README.md"
    rows = [
        {
            "tool": "dotmatch_count",
            "workflow": "synthetic_oligo_adapter_fixture",
            "status": "smoke",
            "n_targets": "3",
            "n_reads": "6",
            "target_start": "0",
            "target_length": "12",
            "k": "1",
            "metric": "hamming",
            "assigned_unique": "4",
            "assigned_exact": "3",
            "corrected_reads": "1",
            "ambiguous_reads": "1",
            "unmatched_reads": "1",
            "validation_mismatches": "0",
            "command": "dotmatch count --targets adapter_oligos.tsv",
        }
    ]

    bench.write_report(rows, report)

    text = report.read_text(encoding="utf-8")
    assert "Oligo/Adapter Assignment Evidence" in text
    assert "synthetic_oligo_adapter_fixture" in text
    assert "not adapter trimming" in text
    assert "make oligo-adapter-smoke-gate" in text


def test_write_report_handles_public_adapter_rows(tmp_path):
    bench = _load_bench()
    report = tmp_path / "README.md"
    rows = [
        {
            "tool": "dotmatch_count",
            "workflow": "synthetic_oligo_adapter_fixture",
            "status": "smoke",
            "n_targets": "3",
            "n_reads": "6",
            "target_start": "0",
            "target_length": "12",
            "k": "1",
            "metric": "hamming",
            "assigned_unique": "4",
            "assigned_exact": "3",
            "corrected_reads": "1",
            "ambiguous_reads": "1",
            "unmatched_reads": "1",
            "validation_mismatches": "0",
            "command": "dotmatch count --targets adapter_oligos.tsv",
        },
        {
            "tool": "dotmatch_count",
            "workflow": "public_fast_adapter_truseq_r1",
            "status": "supported",
            "n_targets": "10",
            "n_reads": "10000",
            "target_start": "229",
            "target_length": "20",
            "k": "0",
            "metric": "hamming",
            "assigned_unique": "76",
            "assigned_exact": "76",
            "corrected_reads": "0",
            "ambiguous_reads": "0",
            "unmatched_reads": "9924",
            "validation_mismatches": "0",
            "command": "dotmatch count --targets adapter_oligos.tsv --target-start 229",
        },
    ]

    bench.write_report(rows, report)

    text = report.read_text(encoding="utf-8")
    assert "public_fast_adapter_truseq_r1" in text
    assert "exact-slice hash baseline" in text
    assert "not adapter trimming" in text
