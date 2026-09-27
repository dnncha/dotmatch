#!/usr/bin/env python3
"""Cross-path FASTQ structure regressions for the native CLI."""

from __future__ import annotations

import gzip
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("DOTMATCH_BIN", ROOT / "dotmatch")).resolve()
FASTQ_IO_SPEC = importlib.util.spec_from_file_location("dotmatch_fastq_io", ROOT / "python/dotmatch/fastq_io.py")
if FASTQ_IO_SPEC is None or FASTQ_IO_SPEC.loader is None:
    raise RuntimeError("could not load Python FASTQ oracle")
FASTQ_IO = importlib.util.module_from_spec(FASTQ_IO_SPEC)
sys.modules[FASTQ_IO_SPEC.name] = FASTQ_IO
FASTQ_IO_SPEC.loader.exec_module(FASTQ_IO)
iter_fastq = FASTQ_IO.iter_fastq
TARGETS = b"g\tACGTACGT\n"
VALID = b"@ok\nACGTACGT\n+ok\nIIIIIIII\n"


INVALID = {
    "missing-id": (b"@\nACGTACGT\n+\nIIIIIIII\n", 1),
    "control-id": (b"@bad\x01\nACGTACGT\n+\nIIIIIIII\n", 1),
    "nonascii-id": (b"@bad\xff\nACGTACGT\n+\nIIIIIIII\n", 1),
    "bad-header": (b"read\nACGTACGT\n+\nIIIIIIII\n", 1),
    "empty-sequence": (b"@bad\n\n+bad\n\n", 2),
    "sequence-space": (b"@bad\nACGT CGT\n+bad\nIIIIIIII\n", 2),
    "sequence-nonascii": (b"@bad\nACGT\xffCGT\n+bad\nIIIIIIII\n", 2),
    "bad-plus": (b"@bad\nACGTACGT\n-bad\nIIIIIIII\n", 3),
    "mismatched-plus-id": (b"@bad\nACGTACGT\n+other\nIIIIIIII\n", 3),
    "quality-space": (b"@bad\nACGTACGT\n+bad\nIIIIIII \n", 4),
    "truncated": (b"@bad\nACGTACGT\n+bad\n", 4),
}


def write_fastq(path: Path, payload: bytes) -> None:
    if path.name.lower().endswith(".gz"):
        with path.open("wb") as raw_handle:
            with gzip.GzipFile(fileobj=raw_handle, mode="wb", mtime=0) as handle:
                handle.write(payload)
    else:
        path.write_bytes(payload)


def count_command(targets: Path, reads: Path, output: Path, *, full_record: bool) -> list[str]:
    command = [
        str(BIN),
        "count",
        "--targets",
        str(targets),
        "--reads",
        str(reads),
        "--sample-label",
        "sample",
        "--target-start",
        "0",
        "--target-length",
        "8",
        "--k",
        "0",
        "--metric",
        "hamming",
        "--format",
        "mageck",
        "--out",
        str(output),
    ]
    if full_record:
        command.extend(["--assignments", str(output.with_suffix(".assignments.tsv"))])
    return command


def main() -> None:
    if not BIN.exists():
        raise AssertionError(f"native CLI does not exist: {BIN}")
    with tempfile.TemporaryDirectory(prefix="dotmatch-fastq-validation-") as raw_tmp:
        tmp = Path(raw_tmp)
        targets = tmp / "targets.tsv"
        targets.write_bytes(TARGETS)

        for name, (bad_record, line_offset) in INVALID.items():
            for suffix in (".fastq", ".FASTQ.GZ"):
                reads = tmp / f"{name}{suffix}"
                write_fastq(reads, VALID + bad_record)
                try:
                    list(iter_fastq(reads))
                except (UnicodeDecodeError, ValueError, EOFError, gzip.BadGzipFile):
                    pass
                else:
                    raise AssertionError(f"Python parser accepted {reads.name}")

                for full_record in (False, True):
                    mode = "full" if full_record else "sequence"
                    output = tmp / f"{name}-{suffix.replace('.', '_')}-{mode}.tsv"
                    assignments = output.with_suffix(".assignments.tsv")
                    completed = subprocess.run(
                        count_command(targets, reads, output, full_record=full_record),
                        text=True,
                        capture_output=True,
                        check=False,
                    )
                    if completed.returncode == 0:
                        raise AssertionError(f"native {mode} path accepted {reads.name}")
                    expected = f"record 2 (line {4 + line_offset})"
                    if expected not in completed.stderr:
                        raise AssertionError(
                            f"native {mode} path lost context for {reads.name}:\n{completed.stderr}"
                        )
                    if output.exists() or assignments.exists():
                        raise AssertionError(f"failed {mode} parse left a completed-looking output for {reads.name}")

        valid_reads = tmp / "valid.FASTQ.GZ"
        long_description = b"x" * 20000
        valid_payload = (
            b"@   ok description\r\nACGTACGT\r\n+ok note\r\nIIIIIIII\r\n"
            + b"@long "
            + long_description
            + b"\nACGTACGT\n+long\nIIIIIIII\n"
        )
        write_fastq(valid_reads, valid_payload)
        assert [record.read_id for record in iter_fastq(valid_reads)] == ["ok", "long"]
        output = tmp / "valid.tsv"
        completed = subprocess.run(
            count_command(targets, valid_reads, output, full_record=False),
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode != 0:
            raise AssertionError(f"native path rejected valid CRLF/long-header gzip:\n{completed.stderr}")
        if "g\t\t2" not in output.read_text(encoding="utf-8"):
            raise AssertionError("valid uppercase-extension gzip did not produce two exact counts")

        empty_reads = tmp / "empty.fastq"
        empty_reads.write_bytes(b"")
        assert list(iter_fastq(empty_reads)) == []
        empty_output = tmp / "empty.tsv"
        completed = subprocess.run(
            count_command(targets, empty_reads, empty_output, full_record=False),
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode != 0 or "g\t\t0" not in empty_output.read_text(encoding="utf-8"):
            raise AssertionError(f"empty FASTQ behavior diverged:\n{completed.stderr}")

        for name, payload in (
            ("corrupt.fastq.gz", b"not a gzip stream"),
            ("truncated.fastq.gz", (tmp / "valid.FASTQ.GZ").read_bytes()[:-8]),
        ):
            reads = tmp / name
            reads.write_bytes(payload)
            try:
                list(iter_fastq(reads))
            except (UnicodeDecodeError, ValueError, EOFError, OSError, gzip.BadGzipFile):
                pass
            else:
                raise AssertionError(f"Python parser accepted {name}")
            bad_output = tmp / f"{name}.tsv"
            completed = subprocess.run(
                count_command(targets, reads, bad_output, full_record=False),
                text=True,
                capture_output=True,
                check=False,
            )
            if completed.returncode == 0 or bad_output.exists():
                raise AssertionError(f"native path accepted corrupt gzip {name}")

    print("native FASTQ validation parity: PASS")


if __name__ == "__main__":
    main()
