from __future__ import annotations

import io

import pytest

from dotmatch.fastq_io import iter_fastq_records


@pytest.mark.parametrize("codepoint", [*range(256), 0x100, 0x2028, 0xD800, 0xFFFF, 0x10FFFF])
@pytest.mark.parametrize("field", ["sequence", "quality"])
def test_fastq_printable_ascii_contract(codepoint, field):
    character = chr(codepoint)
    seq = "A" + character + "T" if field == "sequence" else "ACT"
    qual = "I" + character + "I" if field == "quality" else "III"
    records = iter_fastq_records(io.StringIO(f"@read\n{seq}\n+\n{qual}\n"), "fixture.fastq")
    if 33 <= codepoint <= 126:
        record, = records
        assert record.seq == seq.upper()
        assert record.qual == qual
    else:
        with pytest.raises(ValueError, match="record 1"):
            list(records)
