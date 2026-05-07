#!/usr/bin/env python3
"""Run external demultiplexing competitors where their semantics fit.

This is intentionally conservative: Cutadapt is benchmarked only for anchored
5' barcode matching with substitutions-only semantics (`--no-indels`), which is
not identical to DotMatch edit-distance assignment when indels are allowed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import subprocess
import tempfile
import time
from pathlib import Path


def read_barcodes(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    with path.open() as fh:
        for i, line in enumerate(fh):
            line = line.rstrip("\n\r")
            if not line:
                continue
            if "\t" in line:
                ident, seq = line.split("\t", 1)
            else:
                ident, seq = str(i), line
            rows.append((ident, seq))
    return rows


def write_cutadapt_fasta(path: Path, barcodes: list[tuple[str, str]]) -> None:
    with path.open("w") as fh:
        for ident, seq in barcodes:
            fh.write(f">{ident}\n^{seq}\n")


def write_target_fasta(path: Path, barcodes: list[tuple[str, str]]) -> None:
    with path.open("w") as fh:
        for ident, seq in barcodes:
            fh.write(f">{ident}\n{seq}\n")


def extract_fastq(reads: Path, out: Path, start: int, length: int) -> None:
    opener = gzip.open if reads.suffix == ".gz" else open
    with opener(reads, "rt") as inp, out.open("w") as fh:
        while True:
            header = inp.readline()
            if not header:
                break
            seq = inp.readline()
            plus = inp.readline()
            qual = inp.readline()
            if not seq or not plus or not qual:
                raise RuntimeError("truncated FASTQ")
            seq = seq.rstrip("\n\r")
            qual = qual.rstrip("\n\r")
            end = start + length
            if end <= len(seq):
                obs = seq[start:end]
                qobs = qual[start:end]
            else:
                obs = ""
                qobs = ""
            fh.write(header)
            fh.write(obs + "\n")
            fh.write(plus)
            fh.write(qobs + "\n")


def count_fastq(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return sum(1 for _ in fh) // 4


def run_command(cmd: list[str]) -> tuple[float, int]:
    t0 = time.perf_counter()
    rc = subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode
    return time.perf_counter() - t0, rc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--barcodes", required=True, type=Path)
    parser.add_argument("--reads", required=True, type=Path)
    parser.add_argument("--barcode-start", required=True, type=int)
    parser.add_argument("--barcode-length", required=True, type=int)
    parser.add_argument("--k", required=True, type=int, choices=[0, 1])
    parser.add_argument("--dotmatch", default="./dotmatch")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--run-cutadapt", action="store_true")
    parser.add_argument("--run-bowtie2", action="store_true")
    parser.add_argument("--run-guide-counter", action="store_true")
    args = parser.parse_args()

    barcodes = read_barcodes(args.barcodes)
    n_reads = count_fastq(args.reads)
    rows: list[dict[str, str]] = []

    with tempfile.TemporaryDirectory(prefix="dotmatch-competitors-") as tmp:
        tmpdir = Path(tmp)
        extracted = tmpdir / "extracted.fastq"
        extract_fastq(args.reads, extracted, args.barcode_start, args.barcode_length)
        dotmatch_out = tmpdir / "dotmatch.tsv"
        seconds, rc = run_command(
            [
                args.dotmatch,
                "fastq-assign",
                "--barcodes",
                str(args.barcodes),
                "--reads",
                str(args.reads),
                "--barcode-start",
                str(args.barcode_start),
                "--barcode-length",
                str(args.barcode_length),
                "--k",
                str(args.k),
                "--out",
                str(dotmatch_out),
            ]
        )
        rows.append(
            {
                "tool": "dotmatch_fastq_assign",
                "semantics": "edit_distance",
                "n_reads": str(n_reads),
                "n_targets": str(len(barcodes)),
                "barcode_length": str(args.barcode_length),
                "k": str(args.k),
                "seconds": f"{seconds:.6f}",
                "reads_per_sec": f"{n_reads / seconds:.1f}" if rc == 0 else "0.0",
                "exit_code": str(rc),
            }
        )

        if args.run_cutadapt:
            cutadapt = shutil.which("cutadapt")
            if cutadapt is None:
                raise SystemExit("cutadapt not found on PATH")
            fasta = tmpdir / "barcodes.fasta"
            write_cutadapt_fasta(fasta, barcodes)
            error_rate = args.k / args.barcode_length
            seconds, rc = run_command(
                [
                    cutadapt,
                    "-e",
                    f"{error_rate:.8f}",
                    "--no-indels",
                    "-g",
                    f"file:{fasta}",
                    "--action=none",
                    "-o",
                    "/dev/null",
                    str(extracted),
                ]
            )
            rows.append(
                {
                    "tool": "cutadapt_anchored_no_indels",
                    "semantics": "hamming_distance",
                    "n_reads": str(n_reads),
                    "n_targets": str(len(barcodes)),
                    "barcode_length": str(args.barcode_length),
                    "k": str(args.k),
                    "seconds": f"{seconds:.6f}",
                    "reads_per_sec": f"{n_reads / seconds:.1f}" if rc == 0 else "0.0",
                    "exit_code": str(rc),
                }
            )

        if args.run_bowtie2:
            bowtie2 = shutil.which("bowtie2")
            bowtie2_build = shutil.which("bowtie2-build")
            if bowtie2 is None or bowtie2_build is None:
                raise SystemExit("bowtie2 and bowtie2-build must be on PATH")
            fasta = tmpdir / "targets.fasta"
            index_prefix = tmpdir / "targets"
            sam = tmpdir / "bowtie2.sam"
            write_target_fasta(fasta, barcodes)
            _, build_rc = run_command([bowtie2_build, str(fasta), str(index_prefix)])
            seed_len = str(min(args.barcode_length, 20))
            seconds, rc = run_command(
                [
                    bowtie2,
                    "-x",
                    str(index_prefix),
                    "-U",
                    str(extracted),
                    "-S",
                    str(sam),
                    "-N",
                    str(args.k),
                    "-L",
                    seed_len,
                    "--end-to-end",
                    "--very-sensitive",
                ]
            )
            if build_rc != 0:
                rc = build_rc
            rows.append(
                {
                    "tool": "bowtie2_extracted_end_to_end",
                    "semantics": "aligner_workflow",
                    "n_reads": str(n_reads),
                    "n_targets": str(len(barcodes)),
                    "barcode_length": str(args.barcode_length),
                    "k": str(args.k),
                    "seconds": f"{seconds:.6f}",
                    "reads_per_sec": f"{n_reads / seconds:.1f}" if rc == 0 else "0.0",
                    "exit_code": str(rc),
                }
            )

        if args.run_guide_counter:
            guide_counter = shutil.which("guide-counter")
            if guide_counter is None:
                raise SystemExit("guide-counter not found on PATH")
            raise SystemExit(
                "guide-counter comparison requires CRISPR-specific sample/library configuration; "
                "use scripts/run_public_crispr_benchmark.py --run-guide-counter"
            )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "tool",
                "semantics",
                "n_reads",
                "n_targets",
                "barcode_length",
                "k",
                "seconds",
                "reads_per_sec",
                "exit_code",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
