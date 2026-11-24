"""Strict four-line FASTQ input shared by Python workflows.

Sequence symbols compare literally. This parser checks the file structure and
printable ASCII, not whether a read originated from a particular DNA target.
"""
from __future__ import annotations

import gzip
import zlib
from pathlib import Path
from typing import Any, Iterator, NamedTuple, TextIO

NATIVE_MAX_READ_BASES = 8191
_NATIVE_SPACE = " \t\n\r\v\f"


def _valid_printable_ascii(text: str) -> bool:
    return text.isascii() and text.isprintable() and " " not in text


def _native_fields(text: str) -> list[str]:
    fields: list[str] = []
    current: list[str] = []
    for char in text:
        if char in _NATIVE_SPACE:
            if current:
                fields.append("".join(current))
                current = []
        else:
            current.append(char)
    if current:
        fields.append("".join(current))
    return fields


class FastqRecord(NamedTuple):
    read_id: str
    seq: str
    qual: str


def iter_fastq_records(handle: TextIO, source: str | Path, *, content_digest: Any = None) -> Iterator[FastqRecord]:
    """Consume a text handle without closing it; hash before normalization."""
    ordinal = 0
    bare_cr_line = False

    def next_line() -> str:
        nonlocal bare_cr_line
        line = handle.readline()
        if line and bare_cr_line:
            raise ValueError(
                f"FASTQ {source} has a bare carriage-return line ending; lines end with \\n or \\r\\n "
                "(open text handles with newline='\\n')"
            )
        bare_cr_line = line.endswith("\r")
        return line

    while True:
        header = next_line()
        if not header:
            return
        ordinal += 1
        seq, plus, qual = next_line(), next_line(), next_line()
        if content_digest is not None:
            for line in (header, seq, plus, qual):
                content_digest.update(line.encode("utf-8"))
        location = f"in {source} at record {ordinal} (line {4 * ordinal - 3})"
        if not seq or not plus or not qual:
            raise ValueError(f"truncated FASTQ record {location}")
        header, seq, plus, qual = (text.rstrip("\r\n") for text in (header, seq, plus, qual))
        if not header.startswith("@") or not plus.startswith("+"):
            raise ValueError(f"invalid FASTQ record {location}: expected @ header and + separator")
        identifiers = _native_fields(header[1:])
        if not identifiers:
            raise ValueError(f"invalid FASTQ record {location}: missing read identifier")
        if not _valid_printable_ascii(identifiers[0]):
            raise ValueError(f"invalid FASTQ read identifier {location}")
        if not seq or len(seq) != len(qual):
            raise ValueError(f"invalid FASTQ record {location}: sequence and quality lengths differ or are empty")
        if not _valid_printable_ascii(seq):
            raise ValueError(f"invalid FASTQ sequence {location}: expected printable ASCII without whitespace")
        if len(seq) > NATIVE_MAX_READ_BASES:
            raise ValueError(f"invalid FASTQ sequence {location}: sequence is longer than {NATIVE_MAX_READ_BASES} bases")
        if not _valid_printable_ascii(qual):
            raise ValueError(f"invalid Phred+33 quality {location}: expected ASCII 33–126")
        repeated = _native_fields(plus[1:])
        if repeated and repeated[0] != identifiers[0]:
            raise ValueError(f"invalid FASTQ separator {location}: repeated read identifier does not match")
        yield FastqRecord(identifiers[0], seq.upper(), qual)


def iter_fastq(path: str | Path, *, content_digest: Any = None) -> Iterator[FastqRecord]:
    """Read plain/gzipped FASTQ, optionally hashing original decompressed bytes.

    The entire iterator must be consumed before its digest describes the whole
    input. Uppercasing and newline removal happen after digest updates.
    """
    source = Path(path)
    opener = gzip.open if source.name.lower().endswith(".gz") else open
    try:
        with opener(source, "rt", encoding="utf-8", newline="\n") as handle:
            yield from iter_fastq_records(handle, source, content_digest=content_digest)
    except (EOFError, zlib.error, gzip.BadGzipFile) as exc:
        raise ValueError(f"corrupt or truncated gzip FASTQ {source}: {exc}") from exc
