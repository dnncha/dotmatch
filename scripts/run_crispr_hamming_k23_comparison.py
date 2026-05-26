#!/usr/bin/env python3
"""Run DotMatch-vs-Bowtie 1 Hamming k2/k3 CRISPR comparator rows."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
from run_public_crispr_benchmark import (
    ROOT,
    command_text,
    count_fastq_gz,
    dotmatch_stats,
    make_row,
    n_targets,
    run,
)


RAW = ROOT / "benchmarks" / "raw"
DEFAULT_OUT = RAW / "crispr_comparison_hamming_k23_comparators.csv"


def load_manifest(path: Path) -> dict[str, object]:
    if not path.exists():
        raise SystemExit(f"missing dataset manifest: {path}")
    return json.loads(path.read_text())


def sample_labels(manifest: dict[str, object]) -> list[str]:
    return [str(s["sample_id"]) for s in manifest.get("samples", [])]  # type: ignore[index]


def sample_fastqs(manifest: dict[str, object]) -> list[Path]:
    return [Path(str(s["fastq"])) for s in manifest.get("samples", [])]  # type: ignore[index]


def manifest_read_count(manifest: dict[str, object], reads: list[Path]) -> int:
    total = 0
    for sample in manifest.get("samples", []):
        if not isinstance(sample, dict):
            continue
        value = sample.get("written_records", sample.get("expected_full_records"))
        if value is None:
            continue
        try:
            total += int(value)
        except (TypeError, ValueError):
            pass
    return total if total > 0 else sum(count_fastq_gz(read) for read in reads)


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh))


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_ks(text: str) -> list[int]:
    out: list[int] = []
    for part in text.split(","):
        part = part.strip().removeprefix("k")
        if not part:
            continue
        value = int(part)
        if value not in {2, 3}:
            raise ValueError("--ks supports only 2 and 3 for this comparator lane")
        out.append(value)
    if not out:
        raise ValueError("--ks must contain at least one value")
    return out


def command_with_reads(base: list[str], reads: list[Path]) -> list[str]:
    out = list(base)
    for read in reads:
        out.extend(["--reads", str(read)])
    return out


def dotmatch_command(dotmatch: Path, library: Path, reads: list[Path], labels: list[str], target_start: int,
                     target_length: int, k: int, out_dir: Path, threads: int) -> tuple[list[str], Path]:
    summary = out_dir / f"summary.dotmatch.hamming.k{k}.json"
    cmd = command_with_reads([
        str(dotmatch),
        "count",
        "--targets", str(library),
        "--sample-label", ",".join(labels),
        "--target-start", str(target_start),
        "--target-length", str(target_length),
        "--k", str(k),
        "--metric", "hamming",
        "--ambiguity-policy", "best",
        "--format", "mageck",
        "--out", str(out_dir / f"counts.dotmatch.hamming.k{k}.mageck.tsv"),
        "--summary", str(summary),
    ], reads)
    if threads > 1:
        cmd.extend(["--threads", str(threads)])
    return cmd, summary


def bowtie_command(script: Path, library: Path, reads: list[Path], target_start: int, target_length: int,
                   k: int, out_csv: Path, workflow: str, dataset_id: str, requested: str,
                   allow_missing: bool) -> list[str]:
    cmd = [
        sys.executable,
        str(script),
        "--guides", str(library),
        "--target-start", str(target_start),
        "--target-length", str(target_length),
        "--k", str(k),
        "--out", str(out_csv),
        "--workflow", workflow,
        "--dataset-id", dataset_id,
        "--requested-records-per-sample", requested,
        "--run-level", "subsample" if requested != "full" else "full",
    ]
    for read in reads:
        cmd.extend(["--reads", str(read)])
    if allow_missing:
        cmd.append("--allow-missing")
    return cmd


def run_bowtie(cmd: list[str], out_csv: Path) -> dict[str, str]:
    subprocess.run(cmd, cwd=ROOT, check=True)
    rows = read_rows(out_csv)
    if len(rows) != 1:
        raise RuntimeError(f"expected exactly one Bowtie row in {out_csv}, found {len(rows)}")
    return rows[0]


def comparator_row(dataset_id: str, requested: str, k: int, n_reads: int, n_targets_value: int,
                   dotmatch_row: dict[str, str], bowtie_row: dict[str, str], artifact: Path,
                   manifest: Path) -> dict[str, str]:
    dm_rps = float(dotmatch_row.get("reads_per_sec") or 0.0)
    bt_rps = float(bowtie_row.get("reads_per_sec") or 0.0)
    speedup = dm_rps / bt_rps if dm_rps > 0 and bt_rps > 0 else 0.0
    status = "ok" if dotmatch_row.get("exit_code") == "0" and bowtie_row.get("exit_code") == "0" else "missing_comparator"
    return {
        "dataset": dataset_id,
        "dataset_id": dataset_id,
        "k": str(k),
        "records_per_sample": requested,
        "requested_records_per_sample": requested,
        "comparison": f"{dataset_id}:dotmatch_hamming_k{k}_vs_bowtie1",
        "dotmatch_tool": f"dotmatch_hamming_k{k}",
        "bowtie1_tool": f"bowtie1_crispr_hamming_k{k}",
        "comparator": "bowtie1",
        "dotmatch_reads_per_sec": f"{dm_rps:.1f}",
        "bowtie1_reads_per_sec": f"{bt_rps:.1f}",
        "speedup": f"{speedup:.2f}" if speedup else "",
        "status": status,
        "semantics": f"Hamming k={k}, no indels, same-strand fixed guide window",
        "artifact": str(artifact.relative_to(ROOT)) if artifact.is_relative_to(ROOT) else str(artifact),
        "manifest": str(manifest),
        "n_reads": str(n_reads),
        "n_targets": str(n_targets_value),
        "dotmatch_seconds": dotmatch_row.get("seconds", ""),
        "bowtie1_seconds": bowtie_row.get("seconds", ""),
        "dotmatch_assigned_reads": dotmatch_row.get("assigned_reads", ""),
        "bowtie1_assigned_reads": bowtie_row.get("assigned_reads", ""),
        "bowtie1_ambiguous_reads": bowtie_row.get("ambiguous_reads", ""),
        "bowtie1_rejected_reads": bowtie_row.get("rejected_reads", ""),
        "dotmatch_command": dotmatch_row.get("command", ""),
        "bowtie1_command": bowtie_row.get("command", ""),
        "exit_code": "0" if status == "ok" else bowtie_row.get("exit_code", dotmatch_row.get("exit_code", "1")),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out", default=str(DEFAULT_OUT), type=Path)
    parser.add_argument("--ks", default="2,3")
    parser.add_argument("--workflow-name")
    parser.add_argument("--requested-records-per-sample", default="")
    parser.add_argument("--dotmatch-threads", type=int, default=int(os.environ.get("DOTMATCH_COUNT_THREADS", "1")))
    parser.add_argument("--allow-missing-bowtie", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    manifest = load_manifest(args.manifest)
    dataset_id = str(manifest.get("dataset_id", "dataset"))
    workflow = args.workflow_name or f"{dataset_id}_hamming_k23"
    reads = sample_fastqs(manifest)
    labels = sample_labels(manifest)
    if not reads or len(reads) != len(labels):
        raise SystemExit("manifest must contain matching sample IDs and FASTQ paths")
    for read in reads:
        if not read.exists():
            raise SystemExit(f"missing FASTQ: {read}")
    library = Path(str(manifest["library"]))
    target_start = int(manifest.get("target_start", 0))
    target_length = int(manifest.get("guide_length", manifest.get("target_length", 20)))
    n_reads = manifest_read_count(manifest, reads)
    n_target_rows = n_targets(library)
    requested = args.requested_records_per_sample or str(manifest.get("subsample_records", ""))
    if not requested:
        requested = "full" if all("expected_full_records" in s for s in manifest.get("samples", []) if isinstance(s, dict)) else str(n_reads)

    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    out_path = args.out
    rows = read_rows(out_path) if args.resume and out_path.exists() else []
    detail_dir = out_path.parent / "crispr_hamming_k23_details"
    detail_dir.mkdir(parents=True, exist_ok=True)
    work_dir = ROOT / "benchmarks" / "work" / "crispr_hamming_k23"
    work_dir.mkdir(parents=True, exist_ok=True)

    for k in parse_ks(args.ks):
        if args.resume and any(
            row.get("dataset") == dataset_id and row.get("k") == str(k) and
            row.get("records_per_sample") == requested and row.get("status") == "ok"
            for row in rows
        ):
            continue
        rows = [
            row for row in rows
            if not (row.get("dataset") == dataset_id and row.get("k") == str(k) and row.get("records_per_sample") == requested)
        ]
        dm_cmd, dm_summary = dotmatch_command(
            ROOT / "dotmatch", library, reads, labels, target_start, target_length, k, work_dir, args.dotmatch_threads
        )
        dm_seconds, dm_rc, dm_rss = run(dm_cmd, cwd=ROOT)
        dm_row = make_row(
            f"dotmatch_hamming_k{k}", "local", workflow, f"hamming_k{k}_no_indels_fixed_window",
            n_reads, n_target_rows, dm_seconds, dm_rc, dm_rss, dm_cmd, dotmatch_stats(dm_summary)
        )

        bowtie_csv = detail_dir / f"{dataset_id}_hamming_k{k}_{requested}_bowtie1.csv"
        bt_cmd = bowtie_command(
            ROOT / "scripts" / "bench_bowtie1_crispr.py", library, reads, target_start, target_length,
            k, bowtie_csv, workflow, dataset_id, requested, args.allow_missing_bowtie
        )
        bt_row = run_bowtie(bt_cmd, bowtie_csv)
        rows.append(comparator_row(dataset_id, requested, k, n_reads, n_target_rows, dm_row, bt_row, bowtie_csv, args.manifest))
        write_rows(out_path, rows)

    write_rows(out_path, rows)
    print(out_path)


if __name__ == "__main__":
    main()
