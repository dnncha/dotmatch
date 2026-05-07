#!/usr/bin/env python3
"""Benchmark inline-barcode FASTQ demultiplexing workflows.

The script can run on a user-supplied real barcode dataset:

    python3 scripts/bench_barcode_demux.py \
      --reads SRR391079.fastq.gz --barcodes barcodes.tsv \
      --barcode-start 1 --barcode-length auto --k 0 --run-cutadapt

If no reads/barcodes are supplied, it generates a small deterministic fixture.
Those fixture rows are useful for CI and graph plumbing, not for state-of-the-art
claims.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "barcode_demux.csv"
WORK = ROOT / "benchmarks" / "work" / "barcode_demux"


def parse_time_rss(path: Path) -> int:
    if not path.exists():
        return 0
    for line in path.read_text(errors="replace").splitlines():
        stripped = line.strip()
        if sys.platform == "darwin" and "maximum resident set size" in stripped:
            try:
                return int(stripped.split()[0]) // 1024
            except (IndexError, ValueError):
                return 0
        if "Maximum resident set size" in stripped:
            try:
                return int(stripped.rsplit(":", 1)[1].strip())
            except (IndexError, ValueError):
                return 0
    return 0


def timed_command(cmd: list[str], time_out: Path) -> list[str]:
    time_bin = "/usr/bin/time"
    if not Path(time_bin).exists():
        return cmd
    if sys.platform == "darwin":
        return [time_bin, "-l", "-o", str(time_out), *cmd]
    return [time_bin, "-v", "-o", str(time_out), *cmd]


def run(cmd: list[str], allow_missing: bool = False) -> tuple[float, int, int]:
    with tempfile.NamedTemporaryFile(prefix="dotmatch-demux-time-", delete=False) as tmp:
        time_path = Path(tmp.name)
    actual_cmd = timed_command(cmd, time_path)
    start = time.perf_counter()
    try:
        rc = subprocess.run(actual_cmd, cwd=ROOT, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode
    except FileNotFoundError:
        if allow_missing:
            return 0.0, 127, 0
        raise
    seconds = time.perf_counter() - start
    peak_rss_kb = parse_time_rss(time_path)
    time_path.unlink(missing_ok=True)
    return seconds, rc, peak_rss_kb


def count_fastq(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return sum(1 for _ in fh) // 4


def read_barcodes(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            delim = "," if "," in line and "\t" not in line else "\t"
            fields = line.split(delim)
            if fields[0].lower() in {"id", "sample", "sample_id", "sgrna"}:
                continue
            if len(fields) == 1:
                rows.append((str(len(rows)), fields[0].upper()))
            else:
                rows.append((fields[0], fields[1].upper()))
    return rows


def write_cutadapt_fasta(barcodes: Path, out: Path) -> None:
    with out.open("w") as fh:
        for name, seq in read_barcodes(barcodes):
            fh.write(f">{name}\n^{seq}\n")


def safe_filename(name: str) -> str:
    safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in name)
    return safe or "barcode"


def make_fixture(work: Path, records: int) -> tuple[Path, Path, str]:
    work.mkdir(parents=True, exist_ok=True)
    barcodes = work / "fixture_barcodes.tsv"
    reads = work / "fixture_reads.fastq.gz"
    barcode_rows = [
        ("sample_A", "ACGTACGT"),
        ("sample_B", "TTGGAACC"),
        ("sample_C", "GATCTAGC"),
        ("sample_D", "CCGTAATG"),
    ]
    with barcodes.open("w") as fh:
        for name, seq in barcode_rows:
            fh.write(f"{name}\t{seq}\n")
    payload = "GATTACAGATTACA"
    with gzip.open(reads, "wt") as fh:
        for i in range(records):
            name, seq = barcode_rows[i % len(barcode_rows)]
            observed = seq
            if i % 17 == 0:
                observed = seq[:-1] + ("A" if seq[-1] != "A" else "C")
            if i % 43 == 0:
                observed = "NNNNNNNN"
            full = observed + payload
            fh.write(f"@fixture_{i}_{name}\n{full}\n+\n{'I' * len(full)}\n")
    return reads, barcodes, "synthetic_inline_barcode_fixture"


def write_fastq_record(out, header: str, seq: str, plus: str, qual: str) -> None:
    out.write(header)
    out.write(seq)
    out.write(plus)
    out.write(qual)


def hash_splitter_exact(
    reads: Path,
    barcodes: list[tuple[str, str]],
    out_dir: Path,
    barcode_start: int = 0,
) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    by_length: dict[int, dict[str, list[str]]] = {}
    for name, seq in barcodes:
        by_length.setdefault(len(seq), {}).setdefault(seq, []).append(name)
    lengths = sorted(by_length, reverse=True)
    handles = {}
    assigned = 0
    unmatched = 0
    opener = gzip.open if reads.suffix == ".gz" else open
    try:
        unknown = (out_dir / "unknown.fastq").open("w", encoding="utf-8")
        handles["__unknown__"] = unknown
        with opener(reads, "rt", encoding="utf-8") as fh:
            while True:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline()
                plus = fh.readline()
                qual = fh.readline()
                if not seq or not plus or not qual:
                    raise RuntimeError("FASTQ ended mid-record")
                sample = None
                stripped = seq.rstrip("\n")
                observed = stripped[barcode_start:]
                for length in lengths:
                    names = by_length[length].get(observed[:length])
                    if names and len(names) == 1:
                        sample = names[0]
                        break
                if sample is None:
                    unmatched += 1
                    write_fastq_record(unknown, header, seq, plus, qual)
                    continue
                assigned += 1
                if sample not in handles:
                    handles[sample] = (out_dir / f"{safe_filename(sample)}.fastq").open("w", encoding="utf-8")
                write_fastq_record(handles[sample], header, seq, plus, qual)
    finally:
        for handle in handles.values():
            handle.close()
    return {"assigned_reads": str(assigned), "unmatched_reads": str(unmatched)}


def command_text(cmd: list[str]) -> str:
    return " ".join(cmd)


def tool_version(exe: str, args: list[str]) -> str:
    found = shutil.which(exe)
    if found is None:
        return "not_installed"
    try:
        p = subprocess.run([found, *args], check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except Exception:
        return "unknown"
    for line in p.stdout.splitlines():
        if line.strip():
            return line.strip().replace(",", ";")
    return "unknown"


def dotmatch_stats(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    data = json.loads(path.read_text())
    total = int(data.get("total_reads", 0))
    verified = int(data.get("candidates_verified", 0))
    return {
        "assigned_reads": str(data.get("assigned_unique", "")),
        "exact_reads": str(data.get("assigned_exact", "")),
        "corrected_reads": str(data.get("assigned_corrected", "")),
        "ambiguous_reads": str(data.get("ambiguous", "")),
        "unmatched_reads": str(int(data.get("unmatched", 0)) + int(data.get("invalid", 0))),
        "verified_per_read": f"{verified / total:.4f}" if total else "",
    }


def count_demux_outputs(path: Path) -> int:
    total = 0
    for item in path.glob("*.fastq"):
        if "unknown" in item.name:
            continue
        total += count_fastq(item)
    return total


def make_row(tool: str, version: str, workflow: str, semantics: str, reads: int, barcodes: int,
             barcode_len: int | str, k: int, metric: str, seconds: float, rc: int, peak_rss_kb: int,
             command: list[str], repeat: int, stats: dict[str, str] | None = None) -> dict[str, str]:
    row = {
        "tool": tool,
        "version": version,
        "workflow": workflow,
        "semantics": semantics,
        "repeat": str(repeat),
        "n_reads": str(reads),
        "n_barcodes": str(barcodes),
        "barcode_length": str(barcode_len),
        "k": str(k),
        "metric": metric,
        "seconds": f"{seconds:.6f}",
        "reads_per_sec": f"{reads / seconds:.1f}" if seconds > 0 and rc == 0 else "0.0",
        "peak_rss_kb": str(peak_rss_kb),
        "assigned_reads": "",
        "exact_reads": "",
        "corrected_reads": "",
        "ambiguous_reads": "",
        "unmatched_reads": "",
        "verified_per_read": "",
        "exit_code": str(rc),
        "command": command_text(command),
    }
    if stats:
        row.update(stats)
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads")
    parser.add_argument("--barcodes")
    parser.add_argument("--barcode-start", type=int, default=0)
    parser.add_argument("--barcode-length", default="8")
    parser.add_argument("--k", type=int, default=1, choices=[0, 1])
    parser.add_argument("--metric", choices=["hamming", "levenshtein"], default="hamming")
    parser.add_argument("--records", type=int, default=int(os.environ.get("DOTMATCH_BARCODE_RECORDS", "20000")))
    parser.add_argument("--workflow-name", default="")
    parser.add_argument("--run-cutadapt", action="store_true")
    parser.add_argument("--run-hash-splitter", action="store_true")
    parser.add_argument("--repeats", type=int, default=int(os.environ.get("DOTMATCH_BARCODE_REPEATS", "1")))
    parser.add_argument("--out", default=str(RAW))
    args = parser.parse_args()
    if args.repeats < 1:
        raise SystemExit("--repeats must be >= 1")
    barcode_length_arg = str(args.barcode_length)
    auto_barcode_length = barcode_length_arg == "auto"
    if not auto_barcode_length:
        try:
            barcode_length_value = int(barcode_length_arg)
        except ValueError as exc:
            raise SystemExit("--barcode-length must be a positive integer or auto") from exc
        if barcode_length_value <= 0:
            raise SystemExit("--barcode-length must be a positive integer or auto")
    else:
        barcode_length_value = 0
    if auto_barcode_length and args.run_cutadapt and args.k != 0:
        raise SystemExit("--barcode-length auto with Cutadapt is only supported for --k 0 exact-prefix comparisons")

    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    WORK.mkdir(parents=True, exist_ok=True)
    if args.reads and args.barcodes:
        reads = Path(args.reads).resolve()
        barcodes = Path(args.barcodes).resolve()
        workflow = args.workflow_name or "real_inline_barcode_user_supplied"
    else:
        reads, barcodes, workflow = make_fixture(WORK, args.records)

    n_reads = count_fastq(reads)
    barcode_rows = read_barcodes(barcodes)
    n_barcodes = len(barcode_rows)
    out_csv = Path(args.out)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []

    for repeat in range(1, args.repeats + 1):
        dotmatch_out = WORK / f"dotmatch_out_r{repeat}"
        shutil.rmtree(dotmatch_out, ignore_errors=True)
        dotmatch_summary = WORK / f"dotmatch_summary_r{repeat}.json"
        dotmatch_cmd = [
            str(ROOT / "dotmatch"),
            "demux",
            "--barcodes", str(barcodes),
            "--reads", str(reads),
            "--barcode-start", str(args.barcode_start),
            "--barcode-length", barcode_length_arg,
            "--k", str(args.k),
            "--metric", args.metric,
            "--out-dir", str(dotmatch_out),
            "--summary", str(dotmatch_summary),
        ]
        seconds, rc, peak_rss_kb = run(dotmatch_cmd)
        rows.append(make_row(
            "dotmatch_demux",
            "local",
            workflow,
            "fixed_position_unique_ambiguous_nomatch",
            n_reads,
            n_barcodes,
            barcode_length_arg,
            args.k,
            args.metric,
            seconds,
            rc,
            peak_rss_kb,
            dotmatch_cmd,
            repeat,
            dotmatch_stats(dotmatch_summary),
        ))

        if args.run_cutadapt:
            cutadapt_out = WORK / f"cutadapt_out_r{repeat}"
            shutil.rmtree(cutadapt_out, ignore_errors=True)
            cutadapt_out.mkdir(parents=True, exist_ok=True)
            fasta = WORK / "cutadapt_barcodes.fasta"
            write_cutadapt_fasta(barcodes, fasta)
            error_rate = "0" if args.k == 0 else f"{(args.k / barcode_length_value):.8f}"
            cutadapt_cmd = [
                "cutadapt",
                "-e", error_rate,
                "--no-indels",
                *(
                    ["-u", str(args.barcode_start)]
                    if args.barcode_start > 0 else []
                ),
                "-g", f"file:{fasta}",
                "--untrimmed-output", str(cutadapt_out / "unknown.fastq"),
                "-o", str(cutadapt_out / "{name}.fastq"),
                str(reads),
            ]
            seconds, rc, peak_rss_kb = run(cutadapt_cmd, allow_missing=True)
            assigned = count_demux_outputs(cutadapt_out) if rc == 0 else 0
            rows.append(make_row(
                "cutadapt_demux",
                tool_version("cutadapt", ["--version"]),
                workflow,
                "anchored_cutadapt_demux_no_indels",
                n_reads,
                n_barcodes,
                barcode_length_arg,
                args.k,
                "hamming",
                seconds,
                rc,
                peak_rss_kb,
                cutadapt_cmd,
                repeat,
                {"assigned_reads": str(assigned), "unmatched_reads": str(max(0, n_reads - assigned))},
            ))

        if args.run_hash_splitter:
            hash_out = WORK / f"hash_splitter_out_r{repeat}"
            shutil.rmtree(hash_out, ignore_errors=True)
            hash_cmd = [
                "python3", "scripts/bench_barcode_demux.py",
                "--reads", str(reads),
                "--barcodes", str(barcodes),
                "--barcode-start", str(args.barcode_start),
                "--barcode-length", barcode_length_arg,
                "--k", "0",
                "--metric", "hamming",
                "--run-hash-splitter",
                "--repeats", "1",
            ]
            start = time.perf_counter()
            hash_stats = hash_splitter_exact(reads, barcode_rows, hash_out, args.barcode_start)
            seconds = time.perf_counter() - start
            row = make_row(
                "hash_splitter_exact",
                "python_local",
                workflow,
                "longest_unique_exact_prefix_no_mismatch",
                n_reads,
                n_barcodes,
                barcode_length_arg,
                0,
                "exact",
                seconds,
                0,
                0,
                hash_cmd,
                repeat,
                hash_stats,
            )
            row["peak_rss_kb"] = ""
            rows.append(row)

    fields = [
        "tool", "version", "workflow", "semantics", "repeat", "n_reads", "n_barcodes", "barcode_length", "k", "metric",
        "seconds", "reads_per_sec", "peak_rss_kb", "assigned_reads", "exact_reads", "corrected_reads",
        "ambiguous_reads", "unmatched_reads", "verified_per_read", "exit_code", "command",
    ]
    with out_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(out_csv)


if __name__ == "__main__":
    main()
