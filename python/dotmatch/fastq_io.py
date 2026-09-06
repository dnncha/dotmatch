"""Strict four-line FASTQ input shared by Python workflows.

Sequence symbols compare literally. This parser checks the file structure and
printable ASCII, not whether a read originated from a particular DNA target.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, TextIO


@dataclass(frozen=True)
class FastqRecord:
    read_id: str
    seq: str
    qual: str


def iter_fastq_records(handle: TextIO, source: str | Path, *, content_digest: Any = None) -> Iterator[FastqRecord]:
    """Consume a text handle without closing it; hash before normalization."""
    ordinal = 0
    while True:
        header = handle.readline()
        if not header:
            return
        ordinal += 1
        seq, plus, qual = handle.readline(), handle.readline(), handle.readline()
        if content_digest is not None:
            for line in (header, seq, plus, qual):
                content_digest.update(line.encode("utf-8"))
        location = f"in {source} at record {ordinal} (line {4 * ordinal - 3})"
        if not seq or not plus or not qual:
            raise ValueError(f"truncated FASTQ record {location}")
        header, seq, plus, qual = (text.rstrip("\r\n") for text in (header, seq, plus, qual))
        if not header.startswith("@") or not plus.startswith("+"):
            raise ValueError(f"invalid FASTQ record {location}: expected @ header and + separator")
        identifiers = header[1:].split()
        if not identifiers:
            raise ValueError(f"invalid FASTQ record {location}: missing read identifier")
        if any(ord(ch) < 32 or ord(ch) == 127 for ch in identifiers[0]):
            raise ValueError(f"invalid FASTQ read identifier {location}")
        if not seq or len(seq) != len(qual):
            raise ValueError(f"invalid FASTQ record {location}: sequence and quality lengths differ or are empty")
        if not seq.isascii() or any(ord(ch) <= 32 or ord(ch) >= 127 for ch in seq):
            raise ValueError(f"invalid FASTQ sequence {location}: expected printable ASCII without whitespace")
        if any(ord(ch) < 33 or ord(ch) > 126 for ch in qual):
            raise ValueError(f"invalid Phred+33 quality {location}: expected ASCII 33–126")
        if plus[1:].strip() and plus[1:].split()[0] != identifiers[0]:
            raise ValueError(f"invalid FASTQ separator {location}: repeated read identifier does not match")
        yield FastqRecord(identifiers[0], seq.upper(), qual)


def iter_fastq(path: str | Path, *, content_digest: Any = None) -> Iterator[FastqRecord]:
    """Read plain/gzipped FASTQ, optionally hashing original decompressed bytes.

    The entire iterator must be consumed before its digest describes the whole
    input. Uppercasing and newline removal happen after digest updates.
    """
    source = Path(path)
    opener = gzip.open if source.name.lower().endswith(".gz") else open
    with opener(source, "rt", encoding="utf-8", newline="") as handle:
        yield from iter_fastq_records(handle, source, content_digest=content_digest)
