import gzip
import pytest
from dotmatch.core import load_targets
from dotmatch.cli import _read_targets
from dotmatch.target_io import read_target_table


@pytest.mark.parametrize("suffix", [".csv", ".csv.gz", ".CSV.GZ", ".tsv", ".tsv.gz"])
def test_named_columns_gzip_bom_and_order(tmp_path, suffix):
    delimiter = "," if ".csv" in suffix.lower() else "\t"
    content = (
        "\ufeff"
        + delimiter.join(["gene", "sequence", "id"])
        + "\r\n"
        + delimiter.join(["G1", "acgtnry", "a"])
        + "\r\n"
    )
    path = tmp_path / ("library" + suffix)
    path.write_bytes(
        gzip.compress(content.encode())
        if suffix.lower().endswith(".gz")
        else content.encode()
    )
    assert load_targets(path) == [("a", "ACGTNRY")]
    assert _read_targets(path)[0].gene == "G1"


@pytest.mark.parametrize("text", ["ACGT\nTGCA\n", "sequence\nACGT\nTGCA\n"])
def test_sequence_only_keeps_generated_ids(tmp_path, text):
    path = tmp_path / "targets.tsv"
    path.write_text(text)
    assert load_targets(path) == [("target_0", "ACGT"), ("target_1", "TGCA")]


@pytest.mark.parametrize(
    "text",
    [
        "id\tsequence\na\tACGT\na\tTGCA\n",
        "id\tsequence\n\tACGT\n",
        "sequence\tid\nggg\n",
        "id\tsequence\n",
        "id\tsequence\na\t\n",
        "id\tsequence\na\tA C G T\n",
    ],
)
def test_malformed_targets_fail_clearly(tmp_path, text):
    path = tmp_path / "targets.tsv"
    path.write_text(text)
    for reader in (load_targets, _read_targets):
        with pytest.raises(ValueError):
            reader(path)


def test_duplicate_sequences_are_retained(tmp_path):
    path = tmp_path / "targets.tsv"
    path.write_text("id\tsequence\na\tACGT\nb\tACGT\n")
    assert load_targets(path) == [("a", "ACGT"), ("b", "ACGT")]


def test_quoted_csv_and_reordered_headers(tmp_path):
    path = tmp_path / "library.csv"
    path.write_text('gene,id,sequence\n"Gene, annotation",a,ACGT\n')
    assert read_target_table(path)[0].gene == "Gene, annotation"
