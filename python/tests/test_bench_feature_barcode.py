import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "scripts" / "bench_feature_barcode.py"


def _load_bench():
    spec = importlib.util.spec_from_file_location("bench_feature_barcode", BENCH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_make_fixture_includes_feature_whitelist_and_diagnostic_reads(tmp_path):
    bench = _load_bench()

    fixture = bench.make_fixture(tmp_path)

    targets = fixture.targets.read_text(encoding="utf-8")
    reads = fixture.reads.read_text(encoding="utf-8")
    assert "HTO_A\tACGTACGTAA\tcell_hash_A" in targets
    assert "ADT_CD3\tGGGGAAAACC\tCD3" in targets
    assert "HTO_D\tACGTACGTAT\tcell_hash_D" in targets
    assert "@exact_hto_a\nACGTACGTAA" in reads
    assert "@ambiguous_hto_ad\nACGTACGTAC" in reads
    assert "@unmatched\nCCCCCCCCCC" in reads
    assert fixture.expected == {
        "total_reads": 6,
        "assigned_unique": 4,
        "assigned_exact": 4,
        "ambiguous": 1,
        "unmatched": 1,
    }


def test_validation_detects_feature_summary_mismatches(tmp_path):
    bench = _load_bench()
    summary = tmp_path / "summary.json"
    summary.write_text(
        json.dumps(
            {
                "samples": [
                    {
                        "total_reads": 6,
                        "assigned_unique": 4,
                        "assigned_exact": 4,
                        "ambiguous": 1,
                        "unmatched": 1,
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    assert bench.validation_mismatches(summary, {"assigned_unique": 4, "ambiguous": 1}) == []
    assert bench.validation_mismatches(summary, {"assigned_unique": 3, "ambiguous": 1}) == [
        "assigned_unique expected 3 observed 4"
    ]


def test_exact_slice_hash_counts_public_feature_barcodes(tmp_path):
    bench = _load_bench()
    targets = tmp_path / "feature_barcodes.tsv"
    reads = tmp_path / "reads.fastq"
    targets.write_text(
        "target_id\ttarget_seq\tgene\n"
        "CD3\tCTCATTGTAACTCCT\tCD3\n"
        "CD4\tTGTTCCCGCTCAACT\tCD4\n",
        encoding="utf-8",
    )
    reads.write_text(
        "@r0\nAAAAAAAAAACTCATTGTAACTCCTGG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        "@r1\nAAAAAAAAAATGTTCCCGCTCAACTGG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        "@r2\nAAAAAAAAAANNNNNNNNNNNNNNNGG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        encoding="utf-8",
    )

    stats = bench.exact_slice_hash_stats(targets, reads, 10, 15)

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
    targets = tmp_path / "feature_barcodes.tsv"
    reads = tmp_path / "reads.fastq"
    summary = tmp_path / "public_k0_summary.json"
    targets.write_text(
        "target_id\ttarget_seq\tgene\n"
        "CD3\tCTCATTGTAACTCCT\tCD3\n",
        encoding="utf-8",
    )
    reads.write_text(
        "@r0\nAAAAAAAAAACTCATTGTAACTCCTGG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        encoding="utf-8",
    )
    metadata.write_text(
        json.dumps(
            {
                "evidence_ready": True,
                "targets": str(targets),
                "local_fastq": str(reads),
                "target_start": 10,
                "target_length": 15,
                "feature_count": 1,
                "written_records": 1,
                "dataset_id": "10x_1k_pbmc_totalseq_b_3p",
            }
        ),
        encoding="utf-8",
    )

    def fake_run(cmd, cwd, stdout, stderr, text, check):
        assert "--target-start" in cmd
        assert "10" in cmd
        out_summary = tmp_path / ("public_k1_summary.json" if "--k" in cmd and cmd[cmd.index("--k") + 1] == "1" else "public_k0_summary.json")
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
    assert {row["workflow"] for row in rows} == {"public_10x_totalseq_b_feature_barcode"}
    assert rows[0]["k"] == "0"
    assert rows[1]["k"] == "1"
    assert rows[2]["metric"] == "exact"
    assert rows[2]["assigned_unique"] == "1"
    assert all(row["status"] == "supported" for row in rows)


def test_write_report_handles_public_rows_with_target_start(tmp_path):
    bench = _load_bench()
    report = tmp_path / "README.md"
    rows = [
        {
            "tool": "dotmatch_count",
            "workflow": "synthetic_feature_barcode_fixture",
            "status": "smoke",
            "n_targets": "4",
            "n_reads": "6",
            "target_start": "0",
            "target_length": "10",
            "k": "1",
            "metric": "hamming",
            "assigned_unique": "4",
            "assigned_exact": "4",
            "corrected_reads": "0",
            "ambiguous_reads": "1",
            "unmatched_reads": "1",
            "validation_mismatches": "0",
            "command": "dotmatch count --targets feature_barcodes.tsv",
        },
        {
            "tool": "dotmatch_count",
            "workflow": "public_10x_totalseq_b_feature_barcode",
            "status": "supported",
            "n_targets": "10",
            "n_reads": "20000",
            "target_start": "10",
            "target_length": "15",
            "k": "0",
            "metric": "hamming",
            "assigned_unique": "18000",
            "assigned_exact": "18000",
            "corrected_reads": "0",
            "ambiguous_reads": "0",
            "unmatched_reads": "2000",
            "validation_mismatches": "0",
            "command": "dotmatch count --targets feature_barcodes.tsv --target-start 10",
        },
    ]

    bench.write_report(rows, report)

    text = report.read_text(encoding="utf-8")
    assert "public_10x_totalseq_b_feature_barcode" in text
    assert "| dotmatch_count | public_10x_totalseq_b_feature_barcode | supported | 10 | 20000 | 10 | 15 |" in text
