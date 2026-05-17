import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
COMPARE = ROOT / "scripts" / "compare_count_tables.py"


def _load_compare():
    spec = importlib.util.spec_from_file_location("compare_count_tables", COMPARE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_counts(path: Path, guide: str, count: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"sgrna\tGene\tsample\n{guide}\tGENE\t{count}\n", encoding="utf-8")


def test_count_agreement_artifacts_use_repo_relative_paths(tmp_path):
    compare = _load_compare()
    left = ROOT / "examples" / "crispr_guides" / "output" / "left.test.tsv"
    right = ROOT / "examples" / "crispr_guides" / "output" / "right.test.tsv"
    _write_counts(left, "guide_a", 3)
    _write_counts(right, "guide_a", 3)
    try:
        summary, _detail = compare.compare("example", left, right, "left", "right")
    finally:
        left.unlink(missing_ok=True)
        right.unlink(missing_ok=True)

    assert summary["left_path"] == "examples/crispr_guides/output/left.test.tsv"
    assert summary["right_path"] == "examples/crispr_guides/output/right.test.tsv"
    assert str(ROOT) not in summary["left_path"]
    assert str(ROOT) not in summary["right_path"]


def test_makefile_has_real_competitor_benchmark_gate_target():
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")

    assert "bench-real-competitors:" in makefile
    assert "DOTMATCH_PUBLIC_SUBSAMPLE:-10000" in makefile
    assert "scripts/run_public_crispr_benchmark.py --small --run-mageck --run-cutadapt --run-bowtie2 --run-guide-counter" in makefile
    assert "scripts/compare_count_tables.py" in makefile
    assert "scripts/check_public_crispr_claim_gate.py" in makefile
    assert "scripts/bench_barcode_demux.py --reads examples/barcode_demux/data/SRR391079.subsample100000.fastq.gz" in makefile
    assert "scripts/check_barcode_comparison_gate.py" in makefile
