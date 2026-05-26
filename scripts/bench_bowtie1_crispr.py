#!/usr/bin/env python3
"""Bounded Bowtie 1 Hamming comparator for CRISPR guide windows."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path


SEQUENCE_COLUMNS = ("gRNA.sequence", "Seq", "sequence", "Sequence", "guide", "sgRNA", "sgrna", "target")
ID_COLUMNS = ("id", "sgRNAID", "sgrna", "name", "guide", "target_id")
BASE_FIELDS = [
    "tool",
    "version",
    "workflow",
    "semantics",
    "n_reads",
    "n_targets",
    "target_start",
    "target_length",
    "k",
    "seconds",
    "reads_per_sec",
    "exit_code",
    "command",
    "assigned_reads",
    "ambiguous_reads",
    "rejected_reads",
]
METADATA_FIELDS = [
    "dataset_id",
    "sample_id",
    "manifest",
    "repeat",
    "requested_records_per_sample",
    "run_level",
]


def public_text(value: str | Path) -> str:
    text = str(value)
    private_tmp = "/" + "private/tmp/"
    tmp_root = "/" + "tmp/"
    var_folders = "/" + "var/folders/"
    dotmatch_tmp = "/" + "tmp/dotmatch"
    text = text.replace(private_tmp, tmp_root)
    text = re.sub(re.escape(var_folders) + r'[^,\s"]*/([^/,\s"]+)', r"<tmp>/\1", text)
    text = re.sub(re.escape(dotmatch_tmp) + r'[^,\s"]*', r"<tmp>/dotmatch", text)
    return text


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def count_fastq(path: Path) -> int:
    with _open_text(path) as fh:
        return sum(1 for _ in fh) // 4


def extract_guide_fastq(reads: list[Path], out: Path, start: int, length: int) -> int:
    if start < 0:
        raise ValueError("start must be non-negative")
    if length <= 0:
        raise ValueError("length must be positive")
    n_reads = 0
    prefix_headers = len(reads) > 1
    with out.open("w", encoding="utf-8") as out_fh:
        for read_index, read_path in enumerate(reads):
            with _open_text(read_path) as inp:
                while True:
                    header = inp.readline()
                    if not header:
                        break
                    seq = inp.readline()
                    plus = inp.readline()
                    qual = inp.readline()
                    if not seq or not plus or not qual:
                        raise RuntimeError(f"truncated FASTQ: {read_path}")
                    seq = seq.rstrip("\n\r")
                    qual = qual.rstrip("\n\r")
                    end = start + length
                    if end > len(seq) or end > len(qual):
                        guide = ""
                        guide_qual = ""
                    else:
                        guide = seq[start:end]
                        guide_qual = qual[start:end]
                    if prefix_headers and header.startswith("@"):
                        header = f"@{read_index}:{header[1:]}"
                    out_fh.write(header)
                    out_fh.write(guide + "\n")
                    out_fh.write(plus)
                    out_fh.write(guide_qual + "\n")
                    n_reads += 1
    return n_reads


def _dialect_for(path: Path) -> csv.Dialect:
    sample = path.read_text(encoding="utf-8", errors="replace")[:4096]
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t")
    except csv.Error:
        return csv.excel_tab if "\t" in sample.splitlines()[0] else csv.excel


def _first_present(row: dict[str, str], names: tuple[str, ...]) -> str | None:
    for name in names:
        if name in row and row[name].strip():
            return row[name].strip()
    lower = {key.lower(): key for key in row}
    for name in names:
        key = lower.get(name.lower())
        if key is not None and row[key].strip():
            return row[key].strip()
    return None


def read_guides(path: Path) -> list[tuple[str, str]]:
    dialect = _dialect_for(path)
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, dialect=dialect)
        guides: list[tuple[str, str]] = []
        for i, row in enumerate(reader, start=1):
            seq = _first_present(row, SEQUENCE_COLUMNS)
            if not seq:
                raise ValueError(f"could not find guide sequence column in {path}")
            ident = _first_present(row, ID_COLUMNS) or f"guide_{i}"
            guides.append((ident, seq.upper()))
    if not guides:
        raise ValueError(f"no guides found in {path}")
    return guides


def write_fasta(path: Path, guides: list[tuple[str, str]]) -> None:
    seen: set[str] = set()
    with path.open("w", encoding="utf-8") as fh:
        for ident, seq in guides:
            fasta_id = "_".join(ident.split())
            if not fasta_id:
                raise ValueError("empty guide identifier")
            if fasta_id in seen:
                raise ValueError(f"duplicate guide identifier: {fasta_id}")
            seen.add(fasta_id)
            fh.write(f">{fasta_id}\n{seq}\n")


def bowtie_command(bowtie: str, index_prefix: Path, fastq: Path, out: Path, k: int) -> list[str]:
    return [
        bowtie,
        "-q",
        "-v",
        str(k),
        "--best",
        "--strata",
        "--norc",
        "-a",
        str(index_prefix),
        str(fastq),
        str(out),
    ]


def parse_bowtie_assignments(path: Path, n_reads: int) -> dict[str, str]:
    targets_by_read: dict[str, set[str]] = {}
    if path.exists():
        with path.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip():
                    continue
                fields = line.rstrip("\n\r").split("\t")
                if len(fields) < 3:
                    continue
                read_id, target_id = fields[0], fields[2]
                targets_by_read.setdefault(read_id, set()).add(target_id)
    assigned = sum(1 for targets in targets_by_read.values() if len(targets) == 1)
    ambiguous = sum(1 for targets in targets_by_read.values() if len(targets) > 1)
    rejected = max(0, n_reads - assigned - ambiguous)
    return {
        "assigned_reads": str(assigned),
        "ambiguous_reads": str(ambiguous),
        "rejected_reads": str(rejected),
    }


def tool_version(exe: str) -> str:
    resolved = shutil.which(exe)
    if resolved is None:
        return "not_installed"
    try:
        proc = subprocess.run(
            [resolved, "--version"],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
    except Exception:
        return "unknown"
    for line in proc.stdout.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped.replace(",", ";")
    return "unknown"


def command_text(cmd: list[str]) -> str:
    return " ".join(public_text(arg) for arg in cmd)


def make_row(
    args: argparse.Namespace,
    n_reads: int,
    n_targets: int,
    seconds: float,
    exit_code: int,
    command: list[str],
    version: str,
    stats: dict[str, str] | None = None,
) -> dict[str, str]:
    row = {
        "tool": f"bowtie1_crispr_hamming_k{args.k}",
        "version": version,
        "workflow": args.workflow,
        "semantics": f"hamming_k{args.k}_no_indels_bowtie1_v",
        "n_reads": str(n_reads),
        "n_targets": str(n_targets),
        "target_start": str(args.target_start),
        "target_length": str(args.target_length),
        "k": str(args.k),
        "seconds": f"{seconds:.6f}",
        "reads_per_sec": f"{n_reads / seconds:.1f}" if seconds > 0 and exit_code == 0 else "0.0",
        "exit_code": str(exit_code),
        "command": command_text(command),
        "assigned_reads": "",
        "ambiguous_reads": "",
        "rejected_reads": "",
    }
    if stats:
        row.update(stats)
    for field in METADATA_FIELDS:
        value = getattr(args, field, "")
        if value not in (None, ""):
            row[field] = str(value)
    for item in args.metadata:
        if "=" not in item:
            raise ValueError(f"--metadata must be KEY=VALUE, got {item!r}")
        key, value = item.split("=", 1)
        if not key:
            raise ValueError("--metadata key must not be empty")
        row[key] = value
    return row


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = BASE_FIELDS[:]
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_command(cmd: list[str]) -> int:
    return subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--guides", required=True, type=Path, help="DotMatch guide CSV/TSV")
    parser.add_argument("--reads", required=True, action="append", type=Path, help="FASTQ or FASTQ.gz; repeatable")
    parser.add_argument("--target-start", required=True, type=int, help="0-based guide window start in each read")
    parser.add_argument("--target-length", required=True, type=int, help="fixed guide window length")
    parser.add_argument("--k", required=True, type=int, choices=range(0, 4), help="Bowtie -v mismatch bound")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--workflow", default="")
    parser.add_argument("--dataset-id", default="")
    parser.add_argument("--sample-id", default="")
    parser.add_argument("--manifest", default="")
    parser.add_argument("--repeat", default="")
    parser.add_argument("--requested-records-per-sample", default="")
    parser.add_argument("--run-level", default="")
    parser.add_argument("--metadata", action="append", default=[], help="extra CSV metadata as KEY=VALUE")
    parser.add_argument("--bowtie", default="bowtie")
    parser.add_argument("--bowtie-build", default="bowtie-build")
    parser.add_argument("--allow-missing", action="store_true", help="write a not_installed row if Bowtie 1 is absent")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    guides = read_guides(args.guides)
    n_reads = sum(count_fastq(path) for path in args.reads)
    n_targets = len(guides)
    bowtie = shutil.which(args.bowtie)
    bowtie_build = shutil.which(args.bowtie_build)

    with tempfile.TemporaryDirectory(prefix="dotmatch-bowtie1-") as tmp:
        tmpdir = Path(tmp)
        index_prefix = tmpdir / "guides"
        extracted = tmpdir / "guide_windows.fastq"
        bowtie_out = tmpdir / "bowtie.out"
        align_cmd = bowtie_command(args.bowtie, index_prefix, extracted, bowtie_out, args.k)
        if bowtie is None or bowtie_build is None:
            if not args.allow_missing:
                missing = args.bowtie if bowtie is None else args.bowtie_build
                raise SystemExit(f"{missing} not found on PATH; pass --allow-missing to emit a not_installed row")
            row = make_row(args, n_reads, n_targets, 0.0, 127, align_cmd, "not_installed")
            write_rows(args.out, [row])
            print(args.out)
            return

        fasta = tmpdir / "guides.fa"
        build_cmd = [bowtie_build, str(fasta), str(index_prefix)]
        start = time.perf_counter()
        extract_guide_fastq(args.reads, extracted, args.target_start, args.target_length)
        write_fasta(fasta, guides)
        build_rc = run_command(build_cmd)
        align_rc = 0
        stats: dict[str, str] | None = None
        if build_rc == 0:
            align_cmd = bowtie_command(bowtie, index_prefix, extracted, bowtie_out, args.k)
            align_rc = run_command(align_cmd)
            if align_rc == 0:
                stats = parse_bowtie_assignments(bowtie_out, n_reads)
        seconds = time.perf_counter() - start
        exit_code = build_rc if build_rc != 0 else align_rc
        command = [*build_cmd, "&&", *align_cmd]
        row = make_row(args, n_reads, n_targets, seconds, exit_code, command, tool_version(args.bowtie), stats)
        write_rows(args.out, [row])
        print(args.out)


if __name__ == "__main__":
    main()
