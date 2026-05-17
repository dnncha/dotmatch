import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_benchmark_report_references_committed_artifacts() -> None:
    report = (ROOT / "docs" / "benchmarks" / "README.md").read_text(encoding="utf-8")

    linked_svgs = set(re.findall(r"\]\(([^)]+\.svg)\)", report))
    assert {
        "exact_speedup_vs_edlib.svg",
        "threshold_speedup_heatmap.svg",
        "batch_assignment_throughput.svg",
    }.issubset(linked_svgs)
    for asset in linked_svgs:
        assert (ROOT / "docs" / "benchmarks" / asset).is_file()

    assert (ROOT / "docs" / "benchmarks" / "native" / "README.md").is_file()
    assert (ROOT / "docs" / "native-comparator-scope.md").is_file()
