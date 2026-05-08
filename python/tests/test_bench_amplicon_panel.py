import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "scripts" / "bench_amplicon_panel.py"


def _load_bench():
    spec = importlib.util.spec_from_file_location("bench_amplicon_panel", BENCH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_make_fixture_includes_panel_targets_and_diagnostic_reads(tmp_path):
    bench = _load_bench()

    fixture = bench.make_fixture(tmp_path)

    targets = fixture.targets.read_text(encoding="utf-8")
    reads = fixture.reads.read_text(encoding="utf-8")
    assert "AMP_A\tACGTACGTACGT\tGENE_A" in targets
    assert "AMP_D\tACGTACGTACGA\tGENE_D" in targets
    assert "@exact_a\nACGTACGTACGT" in reads
    assert "@ambiguous_ad\nACGTACGTACGG" in reads
    assert "@unmatched\nCCCCCCCCCCCC" in reads
    assert fixture.expected == {
        "total_reads": 6,
        "assigned_unique": 4,
        "assigned_exact": 4,
        "ambiguous": 1,
        "unmatched": 1,
    }


def test_validation_detects_summary_mismatches(tmp_path):
    bench = _load_bench()
    summary = tmp_path / "summary.json"
    summary.write_text(
        json.dumps(
            {
                "total_reads": 6,
                "assigned_unique": 4,
                "assigned_exact": 4,
                "ambiguous": 1,
                "unmatched": 1,
            }
        ),
        encoding="utf-8",
    )

    assert bench.validation_mismatches(summary, {"assigned_unique": 4, "ambiguous": 1}) == []
    assert bench.validation_mismatches(summary, {"assigned_unique": 5, "ambiguous": 1}) == [
        "assigned_unique expected 5 observed 4"
    ]


def test_public_exact_prefix_baseline_counts_amplicon_primers(tmp_path):
    bench = _load_bench()
    targets = tmp_path / "targets.tsv"
    reads = tmp_path / "reads.fastq"
    targets.write_text(
        "target_id\ttarget_seq\tgene\n"
        "primer_1_LEFT\tACGT\tprimer_1\n"
        "primer_2_LEFT\tTTTT\tprimer_2\n",
        encoding="utf-8",
    )
    reads.write_text(
        "@r0\nACGTAAAA\n+\nIIIIIIII\n"
        "@r1\nTTTTCCCC\n+\nIIIIIIII\n"
        "@r2\nGGGGCCCC\n+\nIIIIIIII\n",
        encoding="utf-8",
    )

    stats = bench.exact_prefix_hash_stats(targets, reads, target_length=4)

    assert stats["n_reads"] == "3"
    assert stats["n_targets"] == "2"
    assert stats["assigned_unique"] == "2"
    assert stats["assigned_exact"] == "2"
    assert stats["unmatched_reads"] == "1"
