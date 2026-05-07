#!/usr/bin/env python3
"""Benchmark multi-sample CRISPR counting on real public FASTQ subsamples."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from run_public_crispr_benchmark import count_fastq_gz, dotmatch_stats, guide_counter_stats, run, tool_version


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "examples" / "crispr_guides" / "data"
RAW = ROOT / "benchmarks" / "raw" / "public_crispr_sample_scaling.csv"


def write_samples(path: Path, reads: list[Path], n_samples: int) -> tuple[list[str], list[Path]]:
    labels: list[str] = []
    paths: list[Path] = []
    with path.open("w") as fh:
        fh.write("sample_id\tfastq\n")
        for i in range(n_samples):
            label = f"sample_{i + 1}"
            read = reads[i % len(reads)]
            labels.append(label)
            paths.append(read)
            fh.write(f"{label}\t{read}\n")
    return labels, paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records-per-sample", type=int, default=int(os.environ.get("DOTMATCH_SCALING_RECORDS", "100000")))
    parser.add_argument("--sample-counts", default=os.environ.get("DOTMATCH_SCALING_SAMPLE_COUNTS", "2,4,8"))
    parser.add_argument("--threads", type=int, default=int(os.environ.get("DOTMATCH_COUNT_THREADS", "8")))
    parser.add_argument("--out", default=str(RAW))
    parser.add_argument("--run-guide-counter", action="store_true")
    args = parser.parse_args()

    subprocess.run(["python3", str(ROOT / "scripts" / "fetch_mageck_demo.py"), "--out", str(DATA),
                    "--subsample", str(args.records_per_sample)], cwd=ROOT, check=True)
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)

    reads = [DATA / "ERR376998.fastq.gz", DATA / "ERR376999.fastq.gz"]
    sample_counts = [int(x) for x in args.sample_counts.split(",") if x.strip()]
    rows: list[dict[str, str]] = []

    with tempfile.TemporaryDirectory(prefix="dotmatch-sample-scaling-") as tmp_dir:
        tmp = Path(tmp_dir)
        for n_samples in sample_counts:
            samples_path = tmp / f"samples_{n_samples}.tsv"
            labels, paths = write_samples(samples_path, reads, n_samples)
            total_reads = sum(count_fastq_gz(p) for p in paths)
            dot_summary = tmp / f"dotmatch_{n_samples}.summary.json"
            dot_cmd = [
                str(ROOT / "dotmatch"),
                "crispr-count",
                "--library", str(DATA / "yusa_library.csv"),
                "--samples", str(samples_path),
                "--guide-start", "23",
                "--guide-length", "19",
                "--k", "1",
                "--metric", "hamming",
                "--auto-offset", "5",
                "--auto-offset-sample", str(args.records_per_sample),
                "--threads", str(min(args.threads, n_samples)),
                "--out", str(tmp / f"dotmatch_{n_samples}.counts.tsv"),
                "--summary", str(dot_summary),
            ]
            seconds, rc, peak_rss_kb = run(dot_cmd, cwd=ROOT)
            stats = dotmatch_stats(dot_summary)
            rows.append({
                "tool": "dotmatch_hamming_k1_threaded",
                "version": "local",
                "n_samples": str(n_samples),
                "records_per_sample": str(args.records_per_sample),
                "total_reads": str(total_reads),
                "threads": str(min(args.threads, n_samples)),
                "seconds": f"{seconds:.6f}",
                "reads_per_sec": f"{total_reads / seconds:.1f}" if seconds > 0 and rc == 0 else "0.0",
                "peak_rss_kb": str(peak_rss_kb),
                "assigned_reads": stats.get("assigned_reads", ""),
                "overcount_reads": "0",
                "exit_code": str(rc),
                "command": " ".join(dot_cmd),
            })

            if args.run_guide_counter:
                guide_counter = shutil.which("guide-counter")
                if guide_counter is None:
                    rows.append({
                        "tool": "guide_counter_one_mismatch",
                        "version": "not_installed",
                        "n_samples": str(n_samples),
                        "records_per_sample": str(args.records_per_sample),
                        "total_reads": str(total_reads),
                        "threads": "1",
                        "seconds": "0.000000",
                        "reads_per_sec": "0.0",
                        "peak_rss_kb": "0",
                        "assigned_reads": "",
                        "overcount_reads": "",
                        "exit_code": "127",
                        "command": "guide-counter count",
                    })
                else:
                    prefix = tmp / f"guide_counter_{n_samples}"
                    cmd = [
                        guide_counter,
                        "count",
                        "--input", *[str(p) for p in paths],
                        "--samples", *labels,
                        "--library", str(DATA / "yusa_library.csv"),
                        "--output", str(prefix),
                        "--offset-sample-size", str(args.records_per_sample),
                    ]
                    seconds, rc, peak_rss_kb = run(cmd, cwd=ROOT)
                    stats = guide_counter_stats(tmp / f"guide_counter_{n_samples}.stats.txt")
                    rows.append({
                        "tool": "guide_counter_one_mismatch",
                        "version": tool_version("guide-counter", ["--help"]).replace(",", ";"),
                        "n_samples": str(n_samples),
                        "records_per_sample": str(args.records_per_sample),
                        "total_reads": str(total_reads),
                        "threads": "1",
                        "seconds": f"{seconds:.6f}",
                        "reads_per_sec": f"{total_reads / seconds:.1f}" if seconds > 0 and rc == 0 else "0.0",
                        "peak_rss_kb": str(peak_rss_kb),
                        "assigned_reads": stats.get("assigned_reads", ""),
                        "overcount_reads": stats.get("overcount_reads", ""),
                        "exit_code": str(rc),
                        "command": " ".join(cmd),
                    })

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(out_path)


if __name__ == "__main__":
    main()
