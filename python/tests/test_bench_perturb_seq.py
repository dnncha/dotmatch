import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "scripts" / "bench_perturb_seq.py"


def _load_bench():
    spec = importlib.util.spec_from_file_location("bench_perturb_seq", BENCH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_make_fixture_includes_guide_feature_pairs_and_diagnostics(tmp_path):
    bench = _load_bench()

    fixture = bench.make_fixture(tmp_path)

    guides = fixture.left_targets.read_text(encoding="utf-8")
    features = fixture.right_targets.read_text(encoding="utf-8")
    reads = fixture.reads.read_text(encoding="utf-8")
    assert "GUIDE_A\tACGTAC\tGENE_A" in guides
    assert "GUIDE_D\tACGTAT\tGENE_D" in guides
    assert "HTO_A\tGGAACC\tcell_hash_A" in features
    assert "ADT_CD3\tCCTTAA\tCD3" in features
    assert "@exact_a_hto\nACGTACGGAACC" in reads
    assert "@ambiguous_guide_ad\nACGTAGGGAACC" in reads
    assert "@left_unmatched\nGGGGGGGGAACC" in reads
    assert "@invalid_short\nACGT" in reads
    assert fixture.expected == {
        "total_reads": 7,
        "assigned_pairs": 3,
        "pair_ambiguous": 1,
        "left_unmatched": 1,
        "right_unmatched": 1,
        "invalid": 1,
    }


def test_validation_detects_perturb_seq_summary_mismatches(tmp_path):
    bench = _load_bench()
    summary = tmp_path / "summary.json"
    summary.write_text(
        json.dumps(
            {
                "total_reads": 7,
                "assigned_pairs": 3,
                "pair_ambiguous": 1,
                "left_unmatched": 1,
                "right_unmatched": 1,
                "invalid": 1,
            }
        ),
        encoding="utf-8",
    )

    assert bench.validation_mismatches(summary, {"assigned_pairs": 3, "invalid": 1}) == []
    assert bench.validation_mismatches(summary, {"assigned_pairs": 4, "invalid": 1}) == [
        "assigned_pairs expected 4 observed 3"
    ]


def test_public_exact_slice_baseline_counts_crispr_guides(tmp_path):
    bench = _load_bench()
    targets = tmp_path / "guides.tsv"
    reads = tmp_path / "reads.fastq"
    targets.write_text(
        "target_id\ttarget_seq\tgene\n"
        "RAB1A-2\tTATTTCCTGGTTCGCCGGC\tRAB1A\n",
        encoding="utf-8",
    )
    reads.write_text(
        "@r0\n"
        "CAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCTTAGCTCTTAAACTATTTCCTGGTTCGCCGGCC\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        "@r1\n"
        "CAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCTTAGCTCTTAAACGCCCGCATCGTCAGCACGTT\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        encoding="utf-8",
    )

    stats = bench.exact_slice_hash_stats(targets, reads, target_start=63, target_length=19)

    assert stats["n_reads"] == "2"
    assert stats["n_targets"] == "1"
    assert stats["assigned_unique"] == "1"
    assert stats["assigned_exact"] == "1"
    assert stats["unmatched_reads"] == "1"
