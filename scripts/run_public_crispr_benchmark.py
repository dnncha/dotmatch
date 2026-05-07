#!/usr/bin/env python3
"""Run the MAGeCK/Yusa public CRISPR guide-counting benchmark workflow."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "crispr_guides"


def public_text(value: str | Path) -> str:
    text = str(value)
    root = str(ROOT)
    text = text.replace(root + os.sep, "")
    if text == root:
        text = "."
    private_tmp = "/" + "private/tmp/"
    tmp_root = "/" + "tmp/"
    var_folders = "/" + "var/folders/"
    dotmatch_tmp = "/" + "tmp/dotmatch"
    text = text.replace(private_tmp, tmp_root)
    text = re.sub(re.escape(var_folders) + r'[^,\s"]*/([^/,\s"]+)', r"<tmp>/\1", text)
    text = re.sub(re.escape(dotmatch_tmp) + r'[^,\s"]*/([^/,\s"]+)', r"<tmp>/\1", text)
    return text


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


def run(cmd: list[str], cwd: Path = ROOT, allow_missing: bool = False) -> tuple[float, int, int]:
    time_path: Path | None = None
    with tempfile.NamedTemporaryFile(prefix="dotmatch-time-", delete=False) as tmp:
        time_path = Path(tmp.name)
    actual_cmd = timed_command(cmd, time_path)
    t0 = time.perf_counter()
    try:
        rc = subprocess.run(actual_cmd, cwd=cwd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode
    except FileNotFoundError:
        if allow_missing:
            return 0.0, 127, 0
        raise
    seconds = time.perf_counter() - t0
    peak_rss_kb = parse_time_rss(time_path)
    time_path.unlink(missing_ok=True)
    return seconds, rc, peak_rss_kb


def command_text(cmd: list[str]) -> str:
    return " ".join(public_text(arg) for arg in cmd)


def tool_version(cmd: str, args: list[str]) -> str:
    exe = shutil.which(cmd)
    if exe is None:
        return "not_installed"
    try:
        p = subprocess.run([exe, *args], check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except Exception:
        return "unknown"
    for line in p.stdout.splitlines():
        line = line.strip()
        if line:
            return line.replace(",", ";")
    return "unknown"


def n_targets(path: Path) -> int:
    with path.open() as fh:
        return max(0, sum(1 for line in fh if line.strip()) - 1)


def dotmatch_stats(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    data = json.loads(path.read_text())
    assigned = exact = corrected = ambiguous = rejected = 0
    total_reads = candidates_considered = candidates_verified = 0
    for sample in data.get("samples", []):
        total_reads += int(sample.get("total_reads", 0))
        assigned += int(sample.get("assigned_unique", 0))
        exact += int(sample.get("assigned_exact", 0))
        corrected += int(sample.get("assigned_corrected", 0))
        ambiguous += int(sample.get("ambiguous", 0))
        rejected += int(sample.get("unmatched", 0)) + int(sample.get("invalid", 0))
        candidates_considered += int(sample.get("candidates_considered", 0))
        candidates_verified += int(sample.get("candidates_verified", 0))
    return {
        "assigned_reads": str(assigned),
        "exact_reads": str(exact),
        "corrected_reads": str(corrected),
        "ambiguous_reads": str(ambiguous),
        "rejected_reads": str(rejected),
        "candidates_per_read": f"{candidates_considered / total_reads:.4f}" if total_reads else "",
        "verified_per_read": f"{candidates_verified / total_reads:.4f}" if total_reads else "",
        "offset_mode": str(data.get("offset_mode", "")),
        "hamming_index": str(data.get("hamming_index", "")),
    }


def guide_counter_stats(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    mapped = total = 0
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            total += int(float(row.get("total_reads", 0) or 0))
            mapped += int(float(row.get("mapped_reads", 0) or 0))
    return {
        "assigned_reads": str(mapped),
        "exact_reads": "",
        "corrected_reads": "",
        "ambiguous_reads": "",
        "rejected_reads": str(max(0, total - mapped)),
        "overcount_reads": str(max(0, mapped - total)),
        "candidates_per_read": "",
        "verified_per_read": "",
    }


def mageck_stats(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    mapped = total = 0
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            total += int(float(row.get("Reads", 0) or 0))
            mapped += int(float(row.get("Mapped", 0) or 0))
    return {
        "assigned_reads": str(mapped),
        "exact_reads": str(mapped),
        "corrected_reads": "0",
        "ambiguous_reads": "",
        "rejected_reads": str(max(0, total - mapped)),
        "overcount_reads": "0",
        "candidates_per_read": "",
        "verified_per_read": "",
    }


def make_row(tool: str, version: str, workflow: str, semantics: str, n_reads: int, n_targets_value: int,
        seconds: float, rc: int, peak_rss_kb: int, command: list[str],
        stats: dict[str, str] | None = None) -> dict[str, str]:
    out = {
        "tool": tool,
        "version": version,
        "workflow": workflow,
        "semantics": semantics,
        "n_reads": str(n_reads),
        "n_targets": str(n_targets_value),
        "seconds": f"{seconds:.6f}",
        "reads_per_sec": f"{n_reads / seconds:.1f}" if seconds > 0 and rc == 0 else "0.0",
        "peak_rss_kb": str(peak_rss_kb),
        "assigned_reads": "",
        "exact_reads": "",
        "corrected_reads": "",
        "ambiguous_reads": "",
        "rejected_reads": "",
        "overcount_reads": "0",
        "candidates_per_read": "",
        "verified_per_read": "",
        "exit_code": str(rc),
        "command": command_text(command),
    }
    if stats:
        out.update(stats)
    return out


def count_fastq_gz(path: Path) -> int:
    import gzip

    with gzip.open(path, "rt") as fh:
        return sum(1 for _ in fh) // 4


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--small", action="store_true", help="download a small real FASTQ subsample instead of the full public FASTQ files")
    parser.add_argument("--out", default=str(ROOT / "benchmarks" / "raw" / "public_crispr_workflow.csv"))
    parser.add_argument("--run-mageck", action="store_true")
    parser.add_argument("--run-cutadapt", action="store_true")
    parser.add_argument("--run-bowtie2", action="store_true")
    parser.add_argument("--run-guide-counter", action="store_true")
    parser.add_argument("--skip-levenshtein", action="store_true",
                        help="skip DotMatch Levenshtein timing for Hamming/exact speed-lane runs")
    parser.add_argument("--dotmatch-threads", type=int, default=int(os.environ.get("DOTMATCH_COUNT_THREADS", "1")))
    args = parser.parse_args()

    data = EXAMPLE / "data"
    output = EXAMPLE / "output"
    output.mkdir(parents=True, exist_ok=True)
    fetch = [str(ROOT / "scripts" / "fetch_mageck_demo.py"), "--out", str(data)]
    if args.small:
        fetch.extend(["--subsample", os.environ.get("DOTMATCH_PUBLIC_SUBSAMPLE", "25")])
    subprocess.run(["python3", *fetch], cwd=ROOT, check=True)
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)

    reads = [data / "ERR376998.fastq.gz", data / "ERR376999.fastq.gz"]
    n_reads = sum(count_fastq_gz(p) for p in reads)
    n_target_rows = n_targets(data / "yusa_library.csv")
    workflow = "public_crispr_yusa_small" if args.small else "public_crispr_yusa_full"
    rows: list[dict[str, str]] = []

    dotmatch_exact_summary = output / "summary.exact.json"
    dotmatch_exact_cmd = [
            str(ROOT / "dotmatch"),
            "count",
            "--targets",
            str(data / "yusa_library.csv"),
            "--reads",
            str(reads[0]),
            "--reads",
            str(reads[1]),
            "--sample-label",
            "plasmid,ESC1",
            "--target-start",
            "23",
            "--target-length",
            "19",
            "--k",
            "0",
            "--metric",
            "hamming",
            "--format",
            "mageck",
            "--out",
            str(output / "counts.exact.mageck.tsv"),
            "--summary",
            str(dotmatch_exact_summary),
    ]
    if args.dotmatch_threads > 1:
        dotmatch_exact_cmd.extend(["--threads", str(args.dotmatch_threads)])
    seconds, rc, peak_rss_kb = run(dotmatch_exact_cmd, cwd=ROOT)
    rows.append(make_row("dotmatch_exact_k0", "local",
                    workflow, "exact_k0_no_errors", n_reads, n_target_rows, seconds, rc, peak_rss_kb,
                    dotmatch_exact_cmd, dotmatch_stats(dotmatch_exact_summary)))

    if not args.skip_levenshtein:
        dotmatch_lev_summary = output / "summary.levenshtein.json"
        dotmatch_lev_cmd = [
            str(ROOT / "dotmatch"),
            "count",
            "--targets",
            str(data / "yusa_library.csv"),
            "--reads",
            str(reads[0]),
            "--reads",
            str(reads[1]),
            "--sample-label",
            "plasmid,ESC1",
            "--target-start",
            "23",
            "--target-length",
            "19",
            "--k",
            "1",
            "--metric",
            "levenshtein",
            "--indel-window",
            "1",
            "--auto-offset",
            "5",
            "--auto-offset-sample",
            str(min(100000, max(1, n_reads // 2))),
            "--format",
            "mageck",
            "--out",
            str(output / "counts.levenshtein.mageck.tsv"),
            "--summary",
            str(dotmatch_lev_summary),
        ]
        if args.dotmatch_threads > 1:
            dotmatch_lev_cmd.extend(["--threads", str(args.dotmatch_threads)])
        seconds, rc, peak_rss_kb = run(dotmatch_lev_cmd, cwd=ROOT)
        rows.append(make_row("dotmatch_levenshtein_k1", "local",
                        workflow, "levenshtein_k1_substitution_insertion_deletion", n_reads, n_target_rows, seconds, rc,
                        peak_rss_kb, dotmatch_lev_cmd, dotmatch_stats(dotmatch_lev_summary)))

    dotmatch_ham_summary = output / "summary.hamming.json"
    dotmatch_ham_cmd = [
            str(ROOT / "dotmatch"),
            "count",
            "--targets",
            str(data / "yusa_library.csv"),
            "--reads",
            str(reads[0]),
            "--reads",
            str(reads[1]),
            "--sample-label",
            "plasmid,ESC1",
            "--target-start",
            "23",
            "--target-length",
            "19",
            "--k",
            "1",
            "--metric",
            "hamming",
            "--auto-offset",
            "5",
            "--auto-offset-sample",
            str(min(100000, max(1, n_reads // 2))),
            "--format",
            "mageck",
            "--out",
            str(output / "counts.hamming.mageck.tsv"),
            "--summary",
            str(dotmatch_ham_summary),
    ]
    if args.dotmatch_threads > 1:
        dotmatch_ham_cmd.extend(["--threads", str(args.dotmatch_threads)])
    seconds, rc, peak_rss_kb = run(dotmatch_ham_cmd, cwd=ROOT)
    rows.append(make_row("dotmatch_hamming_k1", "local",
                    workflow, "hamming_k1_no_indels", n_reads, n_target_rows, seconds, rc, peak_rss_kb, dotmatch_ham_cmd,
                    dotmatch_stats(dotmatch_ham_summary)))

    if args.run_mageck:
        mageck = shutil.which("mageck")
        if mageck is not None:
            cmd = [
                    mageck,
                    "count",
                    "-l",
                    str(data / "yusa_library.csv"),
                    "-n",
                    "mageck_exact_benchmark",
                    "--sample-label",
                    "plasmid,ESC1",
                    "--trim-5",
                    "23",
                    "--fastq",
                    str(reads[0]),
                    str(reads[1]),
            ]
            seconds, rc, peak_rss_kb = run(cmd, cwd=output)
            rows.append(make_row("mageck_count_exact", tool_version("mageck", ["--version"]), workflow,
                            "exact_fastq_count_trim5_23", n_reads, n_target_rows, seconds, rc, peak_rss_kb, cmd,
                            mageck_stats(output / "mageck_exact_benchmark.countsummary.txt")))
        else:
            rows.append(make_row("mageck_count_exact", "not_installed", workflow, "exact_fastq_count_trim5_23", n_reads,
                            n_target_rows, 0.0, 127, 0, ["mageck", "count"]))

    if args.run_guide_counter:
        guide_counter = shutil.which("guide-counter")
        if guide_counter is None:
            rows.append(make_row("guide_counter_one_mismatch", "not_installed", workflow, "hamming_k1_no_indels_auto_offset",
                            n_reads, n_target_rows, 0.0, 127, 0, ["guide-counter", "count"]))
        else:
            prefix = output / "guide_counter"
            cmd = [
                guide_counter,
                "count",
                "--input",
                str(reads[0]),
                str(reads[1]),
                "--samples",
                "plasmid",
                "ESC1",
                "--library",
                str(data / "yusa_library.csv"),
                "--output",
                str(prefix),
                "--offset-sample-size",
                str(min(100000, max(1, n_reads // 2))),
            ]
            seconds, rc, peak_rss_kb = run(cmd, cwd=output)
            rows.append(make_row("guide_counter_one_mismatch", "0.1.3", workflow,
                            "hamming_k1_no_indels_auto_offset", n_reads, n_target_rows, seconds, rc, peak_rss_kb, cmd,
                            guide_counter_stats(output / "guide_counter.stats.txt")))

    if args.run_cutadapt or args.run_bowtie2:
        tmp_targets = output / "targets.tsv"
        with (data / "yusa_library.csv").open() as inp, tmp_targets.open("w") as out_fh:
            reader = csv.DictReader(inp)
            for row in reader:
                out_fh.write(f"{row['id']}\t{row['gRNA.sequence']}\n")
        for read in reads:
            comp_out = output / f"competitors_{read.stem}.csv"
            cmd = [
                "python3",
                str(ROOT / "scripts" / "bench_competitors.py"),
                "--barcodes",
                str(tmp_targets),
                "--reads",
                str(read),
                "--barcode-start",
                "23",
                "--barcode-length",
                "19",
                "--k",
                "1",
                "--dotmatch",
                str(ROOT / "dotmatch"),
                "--out",
                str(comp_out),
            ]
            if args.run_cutadapt:
                cmd.append("--run-cutadapt")
            if args.run_bowtie2:
                cmd.append("--run-bowtie2")
            seconds, rc, peak_rss_kb = run(cmd, cwd=ROOT, allow_missing=True)
            rows.append(make_row(f"external_competitors_{read.name}", "see_competitor_csv", workflow,
                            "cutadapt_bowtie2_extracted_workflow", count_fastq_gz(read), n_target_rows, seconds, rc,
                            peak_rss_kb,
                            cmd))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "tool",
                "version",
                "workflow",
                "semantics",
                "n_reads",
                "n_targets",
                "seconds",
                "reads_per_sec",
                "peak_rss_kb",
                "assigned_reads",
                "exact_reads",
                "corrected_reads",
                "ambiguous_reads",
                "rejected_reads",
                "overcount_reads",
                "candidates_per_read",
                "verified_per_read",
                "offset_mode",
                "hamming_index",
                "exit_code",
                "command",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(out_path)


if __name__ == "__main__":
    main()
