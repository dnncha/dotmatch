#!/usr/bin/env python3
"""Benchmark production crispr-count CPU vs experimental Metal on public datasets."""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "crispr_cpu_metal.csv"
FIELDS = [
    "dataset",
    "backend",
    "metal_validate",
    "status",
    "records_per_sample",
    "n_samples",
    "total_reads",
    "n_targets",
    "guide_start",
    "guide_length",
    "k",
    "offset_mode",
    "auto_offset",
    "wall_seconds",
    "phase_total_seconds",
    "reads_per_sec",
    "backend_effective",
    "count_engine",
    "metal_validation",
    "exit_code",
    "count_match_cpu",
    "device",
    "notes",
]


def fnum(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def write_samples(path: Path, rows: list[tuple[str, Path]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("sample_id\tfastq\n")
        for sample_id, fastq in rows:
            fh.write(f"{sample_id}\t{fastq}\n")


def ensure_yusa(records_per_sample: int) -> dict[str, Any]:
    data = ROOT / "examples" / "crispr_guides" / "data"
    cmd = ["python3", str(ROOT / "scripts" / "fetch_mageck_demo.py"), "--out", str(data)]
    if records_per_sample > 0:
        cmd.extend(["--subsample", str(records_per_sample)])
    subprocess.run(cmd, cwd=ROOT, check=True)
    library = data / "yusa_library.csv"
    samples = [
        ("plasmid", data / "ERR376998.fastq.gz"),
        ("ESC1", data / "ERR376999.fastq.gz"),
    ]
    return {
        "dataset": "mageck_yusa",
        "library": library,
        "samples": samples,
        "guide_start": 23,
        "guide_length": 19,
        "metal_eligible": True,
        "offset_mode": "fixed",
        "auto_offset": 0,
        "extra": ["--ambiguity-policy", "best"],
    }


def ensure_sanson(records_per_sample: int) -> dict[str, Any]:
    data = ROOT / "examples" / "crispr_sanson_brunello" / "data"
    cmd = ["python3", str(ROOT / "scripts" / "fetch_sanson_brunello_demo.py"), "--out", str(data)]
    if records_per_sample > 0:
        cmd.extend(["--subsample", str(records_per_sample)])
    subprocess.run(cmd, cwd=ROOT, check=True)
    suffix = ".fastq.gz" if records_per_sample == 0 else f".subsample{records_per_sample}.fastq.gz"
    library = data / "broadgpp-brunello-library-corrected.txt"
    samples = [
        ("plasmid", data / f"plasmid{suffix}"),
        ("RepA", data / f"RepA{suffix}"),
        ("RepB", data / f"RepB{suffix}"),
        ("RepC", data / f"RepC{suffix}"),
    ]
    return {
        "dataset": "sanson_brunello",
        "library": library,
        "samples": samples,
        "guide_start": 20,
        "guide_length": 20,
        "metal_eligible": False,
        "offset_mode": "multi",
        "auto_offset": 20,
        "extra": [
            "--ambiguity-policy", "best",
            "--auto-offset", "20",
            "--auto-offset-sample", str(min(records_per_sample or 100000, 100000)),
            "--offset-mode", "multi",
            "--offset-min-fraction", "0.0025",
        ],
    }


def run_case(
    *,
    dotmatch: Path,
    dataset: dict[str, Any],
    records_per_sample: int,
    backend: str,
    metal_validate: bool,
    threads: int,
    out_dir: Path,
) -> dict[str, str]:
    row = {field: "" for field in FIELDS}
    row.update({
        "dataset": dataset["dataset"],
        "backend": backend,
        "metal_validate": "1" if metal_validate else "0",
        "records_per_sample": str(records_per_sample),
        "n_samples": str(len(dataset["samples"])),
        "guide_start": str(dataset["guide_start"]),
        "guide_length": str(dataset["guide_length"]),
        "k": "1",
        "offset_mode": str(dataset["offset_mode"]),
        "auto_offset": str(dataset["auto_offset"]),
        "device": platform.platform(),
    })
    if backend == "gpu-metal-experimental" and platform.system() != "Darwin":
        row.update({"status": "unavailable", "notes": "Metal is only available on Darwin"})
        return row
    if backend == "gpu-metal-experimental" and not dataset["metal_eligible"]:
        row.update({
            "status": "ineligible",
            "notes": "Metal requires fixed single-offset Hamming k=1; offset-mode multi blocks the production Metal path",
        })
        return row

    out_dir.mkdir(parents=True, exist_ok=True)
    samples_tsv = out_dir / "samples.tsv"
    write_samples(samples_tsv, dataset["samples"])
    counts = out_dir / "counts.mageck.tsv"
    summary = out_dir / "summary.json"
    cmd = [
        str(dotmatch),
        "crispr-count",
        "--library", str(dataset["library"]),
        "--samples", str(samples_tsv),
        "--guide-start", str(dataset["guide_start"]),
        "--guide-length", str(dataset["guide_length"]),
        "--k", "1",
        "--metric", "hamming",
        "--threads", str(threads),
        "--backend", backend,
        "--out", str(counts),
        "--summary", str(summary),
        *dataset["extra"],
    ]
    if metal_validate and backend == "gpu-metal-experimental":
        cmd.append("--metal-validate")
    start = time.perf_counter()
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True)
    wall = time.perf_counter() - start
    row["wall_seconds"] = f"{wall:.6f}"
    row["exit_code"] = str(result.returncode)

    summary_data: dict[str, Any] = {}
    if summary.exists():
        try:
            summary_data = json.loads(summary.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            summary_data = {}
    total_reads = sum(int(sample.get("total_reads", 0) or 0) for sample in summary_data.get("samples", []) or [])
    phase_total = fnum((summary_data.get("phase_seconds") or {}).get("total_before_summary"))
    row.update({
        "total_reads": str(total_reads),
        "n_targets": str(summary_data.get("n_targets", "")),
        "phase_total_seconds": f"{phase_total:.6f}" if phase_total else "",
        "reads_per_sec": f"{(total_reads / wall):.1f}" if wall > 0 and total_reads else "",
        "backend_effective": str(summary_data.get("backend_effective", "")),
        "count_engine": str(summary_data.get("count_engine", "")),
        "metal_validation": str(summary_data.get("metal_validation", "")),
    })
    if result.returncode == 0:
        row["status"] = "ok"
    elif backend == "gpu-metal-experimental" and "Metal backend unavailable" in result.stderr:
        row["status"] = "ineligible"
        row["notes"] = result.stderr.strip().splitlines()[-1] if result.stderr.strip() else "Metal unavailable"
    elif backend == "gpu-metal-experimental" and (
        summary_data.get("metal_validation") == "failed" or "Metal validation failed" in result.stderr
    ):
        row["status"] = "validation_failed"
        row["notes"] = "Metal validation failed against CPU authority checksum"
    else:
        row["status"] = "error"
        row["notes"] = (result.stderr.strip() or result.stdout.strip() or f"exit {result.returncode}")[:500]
    return row


def compare_counts(cpu_counts: Path, other_counts: Path) -> str:
    if not cpu_counts.exists() or not other_counts.exists():
        return ""
    cpu_lines = cpu_counts.read_text(encoding="utf-8").splitlines()
    other_lines = other_counts.read_text(encoding="utf-8").splitlines()
    if cpu_lines == other_lines:
        return "1"
    if len(cpu_lines) != len(other_lines):
        return "0"
    for left, right in zip(cpu_lines[1:], other_lines[1:]):
        if left.split("\t")[2:] != right.split("\t")[2:]:
            return "0"
    return "1"


def write_rows(rows: list[dict[str, str]], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in FIELDS})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dotmatch", default=os.environ.get("DOTMATCH_NATIVE_CLI", str(ROOT / "dotmatch")))
    parser.add_argument("--out", default=str(RAW))
    parser.add_argument(
        "--datasets",
        default=os.environ.get("DOTMATCH_CPU_METAL_DATASETS", "mageck_yusa,sanson_brunello"),
    )
    parser.add_argument(
        "--records-per-sample",
        default=os.environ.get("DOTMATCH_CPU_METAL_READS", "10000,100000"),
        help="comma-separated subsample sizes; use 0 for full FASTQs",
    )
    parser.add_argument("--threads", type=int, default=int(os.environ.get("DOTMATCH_COUNT_THREADS", "1")))
    parser.add_argument("--include-full", action="store_true", default=os.environ.get("DOTMATCH_CPU_METAL_FULL", "") == "1")
    args = parser.parse_args()

    dotmatch = Path(args.dotmatch)
    if not dotmatch.exists():
        raise SystemExit(f"dotmatch binary not found: {dotmatch}")

    sizes = [int(part.strip()) for part in args.records_per_sample.split(",") if part.strip()]
    if args.include_full and 0 not in sizes:
        sizes.append(0)

    fetchers = {
        "mageck_yusa": ensure_yusa,
        "sanson_brunello": ensure_sanson,
    }
    datasets = [name.strip() for name in args.datasets.split(",") if name.strip()]

    rows: list[dict[str, str]] = []
    with tempfile.TemporaryDirectory(prefix="dotmatch_cpu_metal_") as tmp:
        base = Path(tmp)
        for dataset_name in datasets:
            fetch = fetchers.get(dataset_name)
            if fetch is None:
                raise SystemExit(f"unknown dataset: {dataset_name}")
            for records in sizes:
                dataset = fetch(records)
                cpu_dir = base / dataset_name / str(records) / "cpu"
                cpu_row = run_case(
                    dotmatch=dotmatch,
                    dataset=dataset,
                    records_per_sample=records,
                    backend="cpu",
                    metal_validate=False,
                    threads=args.threads,
                    out_dir=cpu_dir,
                )
                rows.append(cpu_row)
                for metal_validate in (False, True):
                    label = "metal_validate" if metal_validate else "metal"
                    metal_dir = base / dataset_name / str(records) / label
                    metal_row = run_case(
                        dotmatch=dotmatch,
                        dataset=dataset,
                        records_per_sample=records,
                        backend="gpu-metal-experimental",
                        metal_validate=metal_validate,
                        threads=args.threads,
                        out_dir=metal_dir,
                    )
                    metal_row["count_match_cpu"] = compare_counts(cpu_dir / "counts.mageck.tsv", metal_dir / "counts.mageck.tsv")
                    rows.append(metal_row)

    write_rows(rows, Path(args.out))
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())