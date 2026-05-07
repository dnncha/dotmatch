#!/usr/bin/env python3
"""Validate CRISPR SOTA datasets against the native Edlib oracle."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import gzip
import json
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw"
WORK = ROOT / "benchmarks" / "work" / "crispr_sota_edlib_validation"
DEFAULT_OUT = RAW / "crispr_sota_edlib_validation.csv"


def run_json(cmd: list[str]) -> dict[str, object]:
    p = subprocess.run(cmd, cwd=ROOT, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode != 0:
        raise SystemExit(p.stderr or p.stdout or f"command failed: {' '.join(cmd)}")
    return json.loads(p.stdout)


def write_fastq_prefix(src: Path, dest: Path, records: int) -> None:
    with gzip.open(src, "rt") as inp, gzip.open(dest, "wt") as out:
        for _ in range(records * 4):
            line = inp.readline()
            if not line:
                break
            out.write(line)


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text())


def as_int(value: str | None, default: int = 0) -> int:
    if value is None or value == "":
        return default
    return int(float(value))


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh))


def write_csv_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames or ["dataset", "sample"])
        writer.writeheader()
        writer.writerows(rows)


def replace_sample_row(rows: list[dict[str, str]], row: dict[str, str]) -> None:
    key = (row.get("dataset", ""), row.get("sample", ""))
    for i, existing in enumerate(rows):
        if (existing.get("dataset", ""), existing.get("sample", "")) == key:
            rows[i] = row
            return
    rows.append(row)


def validation_row_complete(row: dict[str, str], dataset: str, sample: str,
                            min_checked: int, min_records: int) -> bool:
    return (
        row.get("dataset") == dataset
        and row.get("sample") == sample
        and row.get("mismatches") == "0"
        and as_int(row.get("checked_reads")) >= min_checked
        and as_int(row.get("records_available_for_validation")) >= min_records
    )


def find_complete_row(rows: list[dict[str, str]], dataset: str, sample: str,
                      min_checked: int, min_records: int) -> dict[str, str] | None:
    for row in rows:
        if validation_row_complete(row, dataset, sample, min_checked, min_records):
            return row
    return None


def strata(assignments: Path) -> dict[str, int]:
    counts = {
        "exact": 0,
        "corrected": 0,
        "ambiguous": 0,
        "unmatched": 0,
        "contains_n": 0,
        "offset_shift_candidate": 0,
        "indel_window_candidate": 0,
    }
    if not assignments.exists():
        return counts
    with assignments.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            observed = row.get("observed_seq", "")
            status = row.get("status", "")
            correction = row.get("correction", "")
            if "N" in observed:
                counts["contains_n"] += 1
            if status == "ambiguous":
                counts["ambiguous"] += 1
            elif status == "none":
                counts["unmatched"] += 1
            elif correction == "exact":
                counts["exact"] += 1
            elif correction in {"substitution", "insertion", "deletion", "other"}:
                counts["corrected"] += 1
            if correction in {"insertion", "deletion"}:
                counts["indel_window_candidate"] += 1
            if len(observed) > 0 and (observed.startswith("N") or observed.endswith("N")):
                counts["offset_shift_candidate"] += 1
    return counts


def dotmatch_assignment_strata(dataset_id: str, targets: Path, reads: Path, sample: str, target_start: int,
                               target_length: int, auto_offset: int, auto_offset_sample: int,
                               offset_mode: str, offset_min_fraction: float,
                               records: int, tmp: Path) -> tuple[Path, int, str, dict[str, int]]:
    assignments = tmp / f"{dataset_id}.{sample}.assignments.tsv"
    summary = tmp / f"{dataset_id}.{sample}.summary.json"
    counts = tmp / f"{dataset_id}.{sample}.counts.tsv"
    cmd = [
        str(ROOT / "dotmatch"),
        "count",
        "--targets", str(targets),
        "--reads", str(reads),
        "--sample-label", sample,
        "--target-start", str(target_start),
        "--target-length", str(target_length),
        "--k", "1",
        "--metric", "levenshtein",
        "--indel-window", "1",
        "--format", "mageck",
        "--out", str(counts),
        "--assignments", str(assignments),
        "--ambiguous", "report",
        "--summary", str(summary),
    ]
    if auto_offset:
        cmd.extend(["--auto-offset", str(auto_offset), "--auto-offset-sample", str(min(auto_offset_sample, records))])
    cmd.extend(["--offset-mode", offset_mode, "--offset-min-fraction", f"{offset_min_fraction:.8g}"])
    env = os.environ.copy()
    env["DOTMATCH_MAX_READS"] = str(records)
    subprocess.run(cmd, cwd=ROOT, env=env, check=True, stdout=subprocess.DEVNULL)
    selected = target_start
    selected_list = str([target_start])
    data = load_json(summary)
    samples = data.get("samples", [])
    if samples:
        selected = int(samples[0].get("selected_target_start", target_start))  # type: ignore[index]
        selected_list = json.dumps(samples[0].get("selected_target_starts", [selected]))  # type: ignore[index]
    return assignments, selected, selected_list, strata(assignments)


def validate_one(dataset_id: str, sample: str, targets: Path, source_reads: Path, target_start: int,
                 target_length: int, auto_offset: int, auto_offset_sample: int, records: int,
                 offset_mode: str, offset_min_fraction: float, sample_size: int, tmp: Path,
                 edlib_threads: int = 1) -> dict[str, str]:
    validation_reads = tmp / f"{dataset_id}.{sample}.fastq.gz"
    write_fastq_prefix(source_reads, validation_reads, records)
    assignments, selected_start, selected_starts, observed_strata = dotmatch_assignment_strata(
        dataset_id, targets, validation_reads, sample, target_start, target_length,
        auto_offset, auto_offset_sample, offset_mode, offset_min_fraction, records, tmp
    )
    validate_cmd = [
        str(ROOT / "dotmatch"),
        "validate",
        "--targets", str(targets),
        "--reads", str(validation_reads),
        "--target-start", str(target_start),
        "--target-length", str(target_length),
        "--k", "1",
        "--metric", "levenshtein",
        "--indel-window", "1",
        "--oracle", "edlib",
        "--sample", str(sample_size),
        "--threads", str(edlib_threads),
    ]
    if auto_offset:
        validate_cmd.extend(["--auto-offset", str(auto_offset), "--auto-offset-sample", str(min(auto_offset_sample, records))])
    validate_cmd.extend(["--offset-mode", offset_mode, "--offset-min-fraction", f"{offset_min_fraction:.8g}"])
    result = run_json(validate_cmd)
    row = {
        "dataset": dataset_id,
        "sample": sample,
        "oracle": str(result.get("oracle", "edlib_native")),
        "checked_reads": str(result.get("checked_reads", "")),
        "mismatches": str(result.get("mismatches", "")),
        "k": str(result.get("k", "1")),
        "target_start": str(target_start),
        "selected_target_start": str(selected_start),
        "selected_target_starts": selected_starts,
        "target_length": str(target_length),
        "indel_window": "1",
        "offset_mode": offset_mode,
        "offset_min_fraction": f"{offset_min_fraction:.8g}",
        "records_available_for_validation": str(records),
        "assignments_path": str(assignments),
        "edlib_threads": str(edlib_threads),
        "oracle_strategy": str(result.get("oracle_strategy", "")),
        "edlib_alignments": str(result.get("edlib_alignments", "")),
        "bounded_windows": str(result.get("bounded_windows", "")),
        "fallback_windows": str(result.get("fallback_windows", "")),
    }
    row.update({f"stratum_{k}": str(v) for k, v in observed_strata.items()})
    return row


def yusa_dataset(records: int) -> tuple[str, Path, int, int, int, int, str, float, list[tuple[str, Path]]]:
    data = ROOT / "examples" / "crispr_guides" / "data"
    subprocess.run([
        "python3", str(ROOT / "scripts" / "fetch_mageck_demo.py"),
        "--out", str(data), "--subsample", str(records),
    ], cwd=ROOT, check=True)
    return (
        "mageck_yusa",
        data / "yusa_library.csv",
        23,
        19,
        5,
        min(100000, records),
        "best",
        0.005,
        [("plasmid", data / "ERR376998.fastq.gz"), ("ESC1", data / "ERR376999.fastq.gz")],
    )


def sanson_dataset(records: int) -> tuple[str, Path, int, int, int, int, str, float, list[tuple[str, Path]]]:
    data = ROOT / "examples" / "crispr_sanson_brunello" / "data"
    subprocess.run([
        "python3", str(ROOT / "scripts" / "fetch_sanson_brunello_demo.py"),
        "--out", str(data), "--subsample", str(records),
    ], cwd=ROOT, check=True)
    manifest = load_json(data / "manifest.json")
    samples = [(str(s["sample_id"]), Path(str(s["fastq"]))) for s in manifest.get("samples", [])]  # type: ignore[index]
    return (
        str(manifest.get("dataset_id", "sanson_brunello")),
        Path(str(manifest["library"])),
        int(manifest.get("target_start", 20)),
        int(manifest.get("guide_length", 20)),
        int(manifest.get("auto_offset", 20)),
        int(manifest.get("auto_offset_sample", 100000)),
        str(manifest.get("offset_mode", "multi")),
        float(manifest.get("offset_min_fraction", 0.005)),
        samples,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records-per-sample", type=int, default=int(os.environ.get("DOTMATCH_SOTA_VALIDATION_RECORDS", "10000")))
    parser.add_argument("--sample", type=int, default=int(os.environ.get("DOTMATCH_SOTA_VALIDATION_SAMPLE", "10000")))
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--work-dir", default=str(WORK),
                        help="persistent directory for validation FASTQ prefixes and assignment artifacts")
    parser.add_argument("--datasets", default=os.environ.get("DOTMATCH_SOTA_VALIDATION_DATASETS", "mageck_yusa,sanson_brunello"))
    parser.add_argument("--jobs", type=int, default=int(os.environ.get("DOTMATCH_SOTA_VALIDATION_JOBS", "1")),
                        help="number of samples to validate concurrently")
    parser.add_argument("--edlib-threads", type=int,
                        default=int(os.environ.get("DOTMATCH_SOTA_VALIDATION_EDLIB_THREADS", "1")),
                        help="threads per native Edlib validation helper")
    parser.add_argument("--resume", action="store_true",
                        help="reuse completed sample rows already present in --out and checkpoint after each sample")
    args = parser.parse_args()
    if args.records_per_sample <= 0 or args.sample <= 0:
        raise SystemExit("validation record counts must be positive")
    if args.jobs <= 0:
        raise SystemExit("--jobs must be positive")
    if args.edlib_threads <= 0:
        raise SystemExit("--edlib-threads must be positive")

    subprocess.run(["make", "dotmatch", "edlib-tools"], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    out_path = Path(args.out)
    rows: list[dict[str, str]] = read_csv_rows(out_path) if args.resume else []
    dataset_factories = {
        "mageck_yusa": yusa_dataset,
        "sanson_brunello": sanson_dataset,
    }
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    for dataset_name in [d.strip() for d in args.datasets.split(",") if d.strip()]:
        if dataset_name not in dataset_factories:
            raise SystemExit(f"unknown validation dataset: {dataset_name}")
        dataset_id, targets, target_start, target_length, auto_offset, auto_offset_sample, offset_mode, offset_min_fraction, samples = (
            dataset_factories[dataset_name](args.records_per_sample)
        )
        def run_sample(item: tuple[str, Path]) -> dict[str, str]:
            sample, reads = item
            return validate_one(
                dataset_id, sample, targets, reads, target_start, target_length,
                auto_offset, auto_offset_sample, args.records_per_sample,
                offset_mode, offset_min_fraction,
                min(args.sample, args.records_per_sample), work_dir, args.edlib_threads
            )
        sample_size = min(args.sample, args.records_per_sample)
        pending = [
            item for item in samples
            if not args.resume or find_complete_row(rows, dataset_id, item[0], sample_size, args.records_per_sample) is None
        ]
        if args.jobs == 1 or len(pending) <= 1:
            for item in pending:
                row = run_sample(item)
                replace_sample_row(rows, row)
                if args.resume:
                    write_csv_rows(out_path, rows)
        else:
            with ThreadPoolExecutor(max_workers=min(args.jobs, len(samples))) as executor:
                if args.resume:
                    futures = [executor.submit(run_sample, item) for item in pending]
                    for future in as_completed(futures):
                        row = future.result()
                        replace_sample_row(rows, row)
                        write_csv_rows(out_path, rows)
                else:
                    for row in executor.map(run_sample, pending):
                        replace_sample_row(rows, row)

    write_csv_rows(out_path, rows)
    if any(row.get("mismatches") != "0" for row in rows):
        raise SystemExit("CRISPR SOTA Edlib validation mismatches detected")
    print(out_path)


if __name__ == "__main__":
    main()
