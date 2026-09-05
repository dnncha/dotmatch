import json
from pathlib import Path

import pytest

from editwitness.fasta import init_from_fasta, read_fasta
from editwitness.io import InputError
from editwitness.sequence import find_all, reverse_complement

ROOT = Path(__file__).resolve().parents[1]


def test_init_from_local_fasta_maps_reverse_primer(demo):
    seq = demo.reference.sequence
    alt = next(x for x in "ACGT" if x != seq[450])
    manifest = init_from_fasta(ROOT / "examples/synthetic.fasta", seq[200:220], reverse_complement(seq[680:700]), 450, alt)
    assert manifest.assays[0].left_primer.start == 200
    assert manifest.assays[0].right_primer.start == 680
    assert len(manifest.hypotheses) == 3


@pytest.mark.parametrize("text", ["", "ACGT", ">\nACGT", ">one\nACGT\n>two\nTT", ">one\nACNT", ">one\n"])
def test_invalid_fasta(text, tmp_path):
    path = tmp_path / "bad.fasta"
    path.write_text(text)
    with pytest.raises(InputError):
        read_fasta(path)


def test_lowercase_normalized_and_overlapping_matches(tmp_path):
    path = tmp_path / "good.fasta"
    path.write_text(">example\nacgt\nacgt\n")
    assert read_fasta(path) == ("example", "ACGTACGT")
    assert find_all("AAAA", "AA") == (0, 1, 2)
    with pytest.raises(ValueError):
        find_all("ABC", "")


def test_init_rejects_ambiguous_missing_and_invalid_designs(demo):
    path = ROOT / "examples/synthetic.fasta"
    seq = demo.reference.sequence
    l, r = seq[200:220], reverse_complement(seq[680:700])
    for left, right, pos, alt in [("A", r, 450, "T"), ("N", r, 450, "T"),
                                  (l, r, 1000, "T"), (l, r, 450, "TT"), (l, r, 100, "T")]:
        with pytest.raises(InputError):
            init_from_fasta(path, left, right, pos, alt)


def test_fasta_encoding_size_and_reference_budget(tmp_path, monkeypatch):
    from editwitness import fasta
    path = tmp_path / "input.fasta"
    path.write_bytes(b">name\n\xff")
    with pytest.raises(InputError, match="UTF-8"):
        fasta.read_fasta(path)
    path.write_text(">name\n" + "A" * 20001)
    with pytest.raises(InputError, match="20,000"):
        fasta.read_fasta(path)
    monkeypatch.setattr(fasta, "MAX_INPUT_BYTES", 10)
    with pytest.raises(InputError, match="byte budget"):
        fasta.read_fasta(path)
