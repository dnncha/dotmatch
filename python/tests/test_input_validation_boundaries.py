"""Input throughput changes must preserve lossless counts and literal symbols."""
import io
import sys

import pytest

from dotmatch.count_io import parse_count_value, read_count_table
from dotmatch.fastq_io import iter_fastq_records


@pytest.mark.parametrize("field", ["sequence", "quality"])
@pytest.mark.parametrize("character", [chr(i) for i in range(256) if i not in (10, 13)]
                         + ["\u0661", "\u2003", "\uff11", "\U0001f9ec"])
def test_fastq_accepts_exactly_printable_ascii(field, character):
    sequence = f"A{character}T" if field == "sequence" else "ACT"
    quality = f"I{character}I" if field == "quality" else "III"
    text = f"@read description\n{sequence}\n+read\n{quality}\n"
    records = iter_fastq_records(io.StringIO(text), "reads.fastq")
    if 33 <= ord(character) <= 126:
        record, = records
        assert (record.read_id, record.seq, record.qual) == ("read", sequence.upper(), quality)
    else:
        with pytest.raises(ValueError, match="sequence|Phred\\+33"):
            list(records)


@pytest.mark.parametrize("text, expected", [
    ("0007", 7), (" 0007\t", 7), ("0" * 4096, 0),
    ("0" * 3096 + "9" * 1000, 10**1000 - 1), ("9" * 1000, 10**1000 - 1),
    ("9007199254740993", 9007199254740993), ("+0007", 7), ("-0", 0),
    ("007.0", 7), ("7e2", 700),
])
def test_integer_fast_path_preserves_limits_and_decimal_fallback(text, expected):
    assert parse_count_value(text, "guide", "sample") == expected


@pytest.mark.parametrize("text, message", [
    ("1" + "0" * 1000, "count too large"), ("0" * 4097, "count too large"),
    ("\u0661", "non-numeric count"), ("\uff11", "non-numeric count"),
    ("1_000", "non-numeric count"), ("-1", "negative count"),
    ("1.5", "non-integer count"), ("NaN", "non-finite count"),
    ("0e999999", "count too large"),
])
def test_invalid_counts_retain_location_and_error_class(text, message):
    with pytest.raises(ValueError, match=f"{message} for guide/sample"):
        parse_count_value(text, "guide", "sample")


def test_long_counts_ignore_configurable_integer_string_limit():
    if not hasattr(sys, "set_int_max_str_digits"):
        pytest.skip("integer-string conversion limits unavailable on this Python")
    previous = sys.get_int_max_str_digits()
    try:
        sys.set_int_max_str_digits(640)
        assert parse_count_value("0" * 3096 + "9" * 1000, "g", "s") == 10**1000 - 1
    finally:
        sys.set_int_max_str_digits(previous)


def test_unselected_sample_is_still_validated(tmp_path):
    path = tmp_path / "counts.tsv"
    path.write_text("sgRNA\tGene\tA\tB\ng\tG\t7\t\u0661\n", encoding="utf-8")
    with pytest.raises(ValueError, match="non-numeric count for g/B"):
        read_count_table(path, sample_cols=["A"])
