#!/usr/bin/env python3
"""Run one real CRISPR dataset benchmark from a DotMatch dataset manifest."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import tempfile
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
from run_public_crispr_benchmark import (
    ROOT,
    command_text,
    count_fastq_gz,
    dotmatch_stats,
    guide_counter_stats,
    mageck_stats,
    make_row,
    n_targets,
    run,
    tool_version,
)


def load_manifest(path: Path) -> dict[str, object]:
    if not path.exists():
        raise SystemExit(f"missing dataset manifest: {path}")
    return json.loads(path.read_text())


def sample_labels(manifest: dict[str, object]) -> list[str]:
    return [str(s["sample_id"]) for s in manifest.get("samples", [])]  # type: ignore[index]


def sample_scope(labels: list[str]) -> str:
    return labels[0] if len(labels) == 1 else ""


def sample_fastqs(manifest: dict[str, object]) -> list[Path]:
    return [Path(str(s["fastq"])) for s in manifest.get("samples", [])]  # type: ignore[index]


def manifest_read_count(manifest: dict[str, object]) -> int | None:
    total = 0
    samples = manifest.get("samples", [])
    if not isinstance(samples, list) or not samples:
        return None
    for sample in samples:
        if not isinstance(sample, dict):
            return None
        value = sample.get("written_records", sample.get("expected_full_records"))
        if value is None:
            return None
        try:
            total += int(value)
        except (TypeError, ValueError):
            return None
    return total


def output_root(manifest: dict[str, object]) -> Path:
    dataset_id = str(manifest.get("dataset_id", "dataset"))
    return ROOT / "examples" / dataset_id / "output"


def command_with_reads(base: list[str], reads: list[Path]) -> list[str]:
    out = list(base)
    for read in reads:
        out.extend(["--reads", str(read)])
    return out


def guide_counter_command(exe: str, reads: list[Path], labels: list[str], library: Path, prefix: Path,
                          offset_sample: int) -> list[str]:
    return [
        exe,
        "count",
        "--input",
        *[str(r) for r in reads],
        "--samples",
        *labels,
        "--library",
        str(library),
        "--output",
        str(prefix),
        "--offset-sample-size",
        str(offset_sample),
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out", default=str(ROOT / "benchmarks" / "raw" / "crispr_dataset_workflow.csv"))
    parser.add_argument("--workflow-name")
    parser.add_argument("--run-mageck", action="store_true")
    parser.add_argument("--run-guide-counter", action="store_true")
    parser.add_argument("--skip-levenshtein", action="store_true",
                        help="skip DotMatch Levenshtein timing for Hamming/exact speed-lane runs")
    parser.add_argument("--dotmatch-threads", type=int, default=int(os.environ.get("DOTMATCH_COUNT_THREADS", "1")))
    args = parser.parse_args()

    manifest = load_manifest(Path(args.manifest))
    dataset_id = str(manifest.get("dataset_id", "dataset"))
    workflow = args.workflow_name or f"{dataset_id}_real"
    out_dir = output_root(manifest)
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)

    reads = sample_fastqs(manifest)
    labels = sample_labels(manifest)
    if not reads or len(reads) != len(labels):
        raise SystemExit("manifest must contain matching sample IDs and FASTQ paths")
    for read in reads:
        if not read.exists():
            raise SystemExit(f"missing FASTQ: {read}")
    library = Path(str(manifest["library"]))
    n_reads = manifest_read_count(manifest)
    if n_reads is None:
        n_reads = sum(count_fastq_gz(p) for p in reads)
    n_target_rows = n_targets(library)
    target_start = str(manifest.get("target_start", 0))
    guide_length = str(manifest.get("guide_length", 20))
    auto_offset = int(manifest.get("auto_offset", 0) or 0)
    auto_offset_sample = int(manifest.get("auto_offset_sample", min(100000, max(1, n_reads // max(1, len(reads))))))
    offset_mode = str(manifest.get("offset_mode", "best"))
    offset_min_fraction = float(manifest.get("offset_min_fraction", 0.005))
    mageck_trim5 = str(manifest.get("mageck_trim5", target_start))
    label_arg = ",".join(labels)

    rows: list[dict[str, str]] = []

    common = [
        str(ROOT / "dotmatch"),
        "count",
        "--targets", str(library),
        "--sample-label", label_arg,
        "--target-start", target_start,
        "--target-length", guide_length,
        "--format", "mageck",
    ]
    if args.dotmatch_threads > 1:
        common.extend(["--threads", str(args.dotmatch_threads)])
    if auto_offset:
        common.extend(["--auto-offset", str(auto_offset), "--auto-offset-sample", str(auto_offset_sample)])
    common.extend(["--offset-mode", offset_mode, "--offset-min-fraction", f"{offset_min_fraction:.8g}"])

    exact_summary = out_dir / "summary.exact.json"
    exact_cmd = command_with_reads([
        *common,
        "--k", "0",
        "--metric", "hamming",
        "--out", str(out_dir / "counts.exact.mageck.tsv"),
        "--summary", str(exact_summary),
    ], reads)
    seconds, rc, rss = run(exact_cmd, cwd=ROOT)
    rows.append(make_row("dotmatch_exact_k0", "local", workflow, "exact_k0_no_errors",
                         n_reads, n_target_rows, seconds, rc, rss, exact_cmd, dotmatch_stats(exact_summary)))

    ham_summary = out_dir / "summary.hamming.json"
    ham_cmd = command_with_reads([
        *common,
        "--k", "1",
        "--metric", "hamming",
        "--out", str(out_dir / "counts.hamming.mageck.tsv"),
        "--summary", str(ham_summary),
    ], reads)
    seconds, rc, rss = run(ham_cmd, cwd=ROOT)
    rows.append(make_row("dotmatch_hamming_k1", "local", workflow, "hamming_k1_no_indels",
                         n_reads, n_target_rows, seconds, rc, rss, ham_cmd, dotmatch_stats(ham_summary)))

    if not args.skip_levenshtein:
        lev_summary = out_dir / "summary.levenshtein.json"
        lev_cmd = command_with_reads([
            *common,
            "--k", "1",
            "--metric", "levenshtein",
            "--indel-window", "1",
            "--out", str(out_dir / "counts.levenshtein.mageck.tsv"),
            "--summary", str(lev_summary),
        ], reads)
        seconds, rc, rss = run(lev_cmd, cwd=ROOT)
        rows.append(make_row("dotmatch_levenshtein_k1", "local", workflow,
                             "levenshtein_k1_substitution_insertion_deletion",
                             n_reads, n_target_rows, seconds, rc, rss, lev_cmd, dotmatch_stats(lev_summary)))

    if args.run_mageck:
        mageck = shutil.which("mageck")
        if mageck is None:
            rows.append(make_row("mageck_count_exact", "not_installed", workflow,
                                 f"exact_fastq_count_trim5_{mageck_trim5}", n_reads, n_target_rows,
                                 0.0, 127, 0, ["mageck", "count"]))
        else:
            cmd = [
                mageck, "count",
                "-l", str(library),
                "-n", f"{dataset_id}_mageck_exact",
                "--sample-label", label_arg,
                "--trim-5", mageck_trim5,
                "--fastq",
                *[str(r) for r in reads],
            ]
            seconds, rc, rss = run(cmd, cwd=out_dir)
            rows.append(make_row("mageck_count_exact", tool_version("mageck", ["--version"]), workflow,
                                 f"exact_fastq_count_trim5_{mageck_trim5}",
                                 n_reads, n_target_rows, seconds, rc, rss, cmd,
                                 mageck_stats(out_dir / f"{dataset_id}_mageck_exact.countsummary.txt")))

    if args.run_guide_counter:
        gc = shutil.which("guide-counter")
        if gc is None:
            rows.append(make_row("guide_counter_one_mismatch", "not_installed", workflow,
                                 "hamming_k1_no_indels_auto_offset", n_reads, n_target_rows,
                                 0.0, 127, 0, ["guide-counter", "count"]))
        else:
            prefix = out_dir / "guide_counter"
            cmd = guide_counter_command(gc, reads, labels, library, prefix, auto_offset_sample)
            seconds, rc, rss = run(cmd, cwd=out_dir)
            rows.append(make_row("guide_counter_one_mismatch", "0.1.3", workflow,
                                 "hamming_k1_no_indels_auto_offset",
                                 n_reads, n_target_rows, seconds, rc, rss, cmd,
                                 guide_counter_stats(out_dir / "guide_counter.stats.txt")))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    sample_id = sample_scope(labels)
    for row in rows:
        row["dataset_id"] = dataset_id
        row["sample_id"] = sample_id
        row["manifest"] = str(Path(args.manifest))
        row["offset_mode"] = offset_mode
        row["offset_min_fraction"] = f"{offset_min_fraction:.8g}"
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(out_path)


if __name__ == "__main__":
    main()
