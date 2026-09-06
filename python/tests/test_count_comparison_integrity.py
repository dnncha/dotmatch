from scripts.compare_count_tables import parse_counts, compare
import pytest


def test_large_integer_is_not_rounded_through_float(tmp_path):
    path = tmp_path / "counts.tsv"
    path.write_text("sgRNA\tGene\tA\na\tG\t9007199254740993\n")
    assert parse_counts(path)["a"] == 9007199254740993


@pytest.mark.parametrize("value", ["", "NaN", "inf", "2.5", "1e3", "-1", "oops"])
def test_invalid_raw_counts_are_never_silently_ignored(tmp_path, value):
    path = tmp_path / "counts.tsv"
    path.write_text(f"sgRNA\tGene\tA\na\tG\t{value}\n")
    with pytest.raises(ValueError):
        parse_counts(path)


def test_duplicate_guide_ids_are_rejected(tmp_path):
    path = tmp_path / "counts.tsv"
    path.write_text("sgRNA\tGene\tA\na\tG\t1\na\tG\t2\n")
    with pytest.raises(ValueError):
        parse_counts(path)


def test_swapped_sample_counts_do_not_establish_matrix_identity(tmp_path):
    left = tmp_path / "left.tsv"
    right = tmp_path / "right.tsv"
    left.write_text("sgRNA\tGene\tA\tB\na\tG\t1\t9\n")
    right.write_text("sgRNA\tGene\tA\tB\na\tG\t9\t1\n")
    summary, _ = compare("test", left, right, "left", "right")
    assert summary["total_delta"] == "0"
    assert summary["aggregate_counts_identical"] == "true"
    assert summary["counts_identical"] == "false"


def test_named_sample_reordering_is_supported(tmp_path):
    left = tmp_path / "left.tsv"
    right = tmp_path / "right.tsv"
    left.write_text("sgRNA\tGene\tA\tB\na\tG\t1\t9\n")
    right.write_text("sgRNA\tGene\tB\tA\na\tG\t9\t1\n")
    summary, _ = compare("test", left, right, "left", "right")
    assert summary["counts_identical"] == "true"


def test_different_sample_names_require_review(tmp_path):
    left = tmp_path / "left.tsv"
    right = tmp_path / "right.tsv"
    left.write_text("sgRNA\tGene\tA\na\tG\t1\n")
    right.write_text("sgRNA\tGene\tB\na\tG\t1\n")
    summary, _ = compare("test", left, right, "left", "right")
    assert summary["counts_identical"] == ""
    assert summary["matrix_comparability"] == "review_axes_or_annotations"
