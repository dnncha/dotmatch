from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_benchmark_report_documents_native_comparator_scope_boundary() -> None:
    report = (ROOT / "docs" / "benchmarks" / "README.md").read_text(encoding="utf-8")
    generator = (ROOT / "scripts" / "generate_benchmark_report.py").read_text(encoding="utf-8")

    required = (
        "SeqAn/Parasail comparisons are not claimed until "
        "docs/native-comparator-scope.md requirements are met."
    )
    assert required in report
    assert required in generator
    assert "Comparative performance wording should use native C Edlib, SeqAn, and Parasail comparisons" not in report
    assert "Comparative performance wording should use native C Edlib, SeqAn, and Parasail comparisons" not in generator
