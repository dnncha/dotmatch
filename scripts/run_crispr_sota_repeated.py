#!/usr/bin/env python3
"""Run repeated CRISPR SOTA evidence rows across Yusa and Sanson/Brunello."""

from __future__ import annotations

import argparse
import csv
import os
import platform
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw"
DEFAULT_OUT = RAW / "crispr_sota_repeated.csv"
SANSON_DATA = ROOT / "examples" / "crispr_sanson_brunello" / "data"


def parse_sizes(text: str) -> list[int]:
    out: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        value = int(part)
        if value <= 0:
            raise ValueError("read sizes must be positive")
        out.append(value)
    if not out:
        raise ValueError("at least one read size is required")
    return out


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return "unknown"


def machine_metadata() -> dict[str, str]:
    return {
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python": platform.python_version(),
        "dotmatch_commit": git_commit(),
    }


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh))


def expected_tools(args: argparse.Namespace) -> set[str]:
    tools = {"dotmatch_exact_k0", "dotmatch_hamming_k1"}
    if not args.skip_levenshtein:
        tools.add("dotmatch_levenshtein_k1")
    if args.run_mageck:
        tools.add("mageck_count_exact")
    if args.run_guide_counter:
        tools.add("guide_counter_one_mismatch")
    return tools


def should_run_subsamples(args: argparse.Namespace) -> bool:
    return not getattr(args, "full_only", False)


def should_run_full(args: argparse.Namespace) -> bool:
    return bool(getattr(args, "full", False) or getattr(args, "full_only", False))


def selected_sanson_samples(args: argparse.Namespace) -> list[str]:
    text = getattr(args, "sanson_samples", "plasmid,RepA,RepB,RepC")
    return [part.strip() for part in text.split(",") if part.strip()]


def full_sample_scope(dataset: str, args: argparse.Namespace) -> str:
    if dataset != "sanson_brunello":
        return ""
    samples = selected_sanson_samples(args)
    return samples[0] if len(samples) == 1 else ""


def full_run_level(dataset: str, args: argparse.Namespace) -> str:
    return "full_sample" if full_sample_scope(dataset, args) else "full"


def scoped_manifest_path(manifest: Path, sample_id: str) -> Path:
    if not sample_id:
        return manifest
    return manifest.with_name(f"{manifest.stem}.{sample_id}{manifest.suffix}")


def same_run(row: dict[str, str], dataset: str, requested_records: str, repeat: int, run_level: str,
             sample_id: str = "") -> bool:
    return (
        row.get("dataset_id") == dataset
        and row.get("requested_records_per_sample") == requested_records
        and row.get("repeat") == str(repeat)
        and row.get("run_level") == run_level
        and (not sample_id or row.get("sample_id") == sample_id)
    )


def run_complete(rows: list[dict[str, str]], dataset: str, requested_records: str, repeat: int,
                 run_level: str, tools: set[str], sample_id: str = "") -> bool:
    successful = {
        row.get("tool", "")
        for row in rows
        if same_run(row, dataset, requested_records, repeat, run_level, sample_id)
        and row.get("exit_code") == "0"
    }
    return tools.issubset(successful)


def without_run(rows: list[dict[str, str]], dataset: str, requested_records: str, repeat: int,
                run_level: str, sample_id: str = "") -> list[dict[str, str]]:
    return [
        row for row in rows
        if not same_run(row, dataset, requested_records, repeat, run_level, sample_id)
    ]


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_cmd(cmd: list[str], env: dict[str, str] | None = None) -> None:
    subprocess.run(cmd, cwd=ROOT, env=env, check=True)


def fetch_sanson(records_per_sample: int | None, args: argparse.Namespace) -> Path:
    cmd = ["python3", str(ROOT / "scripts" / "fetch_sanson_brunello_demo.py"), "--out", str(SANSON_DATA)]
    if records_per_sample is None:
        cmd.extend(["--subsample", "0"])
    else:
        cmd.extend(["--subsample", str(records_per_sample)])
    if getattr(args, "sanson_samples", ""):
        cmd.extend(["--samples", args.sanson_samples])
    run_cmd(cmd)
    return SANSON_DATA / "manifest.json"


def run_yusa(records_per_sample: int | None, args: argparse.Namespace, tmp_out: Path) -> list[dict[str, str]]:
    cmd = ["python3", str(ROOT / "scripts" / "run_public_crispr_benchmark.py"), "--out", str(tmp_out)]
    env = os.environ.copy()
    if records_per_sample is not None:
        cmd.append("--small")
        env["DOTMATCH_PUBLIC_SUBSAMPLE"] = str(records_per_sample)
    if args.run_mageck:
        cmd.append("--run-mageck")
    if args.run_guide_counter:
        cmd.append("--run-guide-counter")
    if args.skip_levenshtein:
        cmd.append("--skip-levenshtein")
    if args.dotmatch_threads > 1:
        cmd.extend(["--dotmatch-threads", str(args.dotmatch_threads)])
    run_cmd(cmd, env=env)
    rows = read_csv(tmp_out)
    for row in rows:
        row["dataset_id"] = "mageck_yusa"
    return rows


def run_sanson(records_per_sample: int | None, args: argparse.Namespace, tmp_out: Path) -> list[dict[str, str]]:
    manifest = fetch_sanson(records_per_sample, args)
    sample_id = full_sample_scope("sanson_brunello", args) if records_per_sample is None else ""
    scoped_manifest = scoped_manifest_path(manifest, sample_id)
    if scoped_manifest != manifest:
        shutil.copyfile(manifest, scoped_manifest)
        manifest = scoped_manifest
    cmd = [
        "python3", str(ROOT / "scripts" / "run_crispr_dataset_benchmark.py"),
        "--manifest", str(manifest),
        "--out", str(tmp_out),
    ]
    if records_per_sample is None:
        cmd.extend(["--workflow-name", "sanson_brunello_full"])
    else:
        cmd.extend(["--workflow-name", f"sanson_brunello_{records_per_sample}"])
    if args.run_mageck:
        cmd.append("--run-mageck")
    if args.run_guide_counter:
        cmd.append("--run-guide-counter")
    if args.skip_levenshtein:
        cmd.append("--skip-levenshtein")
    if args.dotmatch_threads > 1:
        cmd.extend(["--dotmatch-threads", str(args.dotmatch_threads)])
    run_cmd(cmd)
    return read_csv(tmp_out)


def selected_dataset_runners(args: argparse.Namespace):
    available = {
        "mageck_yusa": run_yusa,
        "sanson_brunello": run_sanson,
    }
    selected = []
    for name in [part.strip() for part in args.datasets.split(",") if part.strip()]:
        if name not in available:
            raise SystemExit(f"unknown dataset: {name}")
        selected.append((name, available[name]))
    if not selected:
        raise SystemExit("at least one dataset is required")
    return selected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--read-sizes", default=os.environ.get("DOTMATCH_SOTA_READ_SIZES", "10000,100000"))
    parser.add_argument("--repeats", type=int, default=int(os.environ.get("DOTMATCH_SOTA_REPEATS", "5")))
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--run-mageck", action="store_true")
    parser.add_argument("--run-guide-counter", action="store_true")
    parser.add_argument("--skip-levenshtein", action="store_true",
                        help="skip DotMatch Levenshtein rows for exact/Hamming speed-lane runs")
    parser.add_argument("--dotmatch-threads", type=int, default=int(os.environ.get("DOTMATCH_COUNT_THREADS", "1")))
    parser.add_argument("--datasets", default=os.environ.get("DOTMATCH_SOTA_DATASETS", "mageck_yusa,sanson_brunello"),
                        help="comma-separated dataset ids to run: mageck_yusa,sanson_brunello")
    parser.add_argument("--sanson-samples", default=os.environ.get("DOTMATCH_SANSON_SAMPLES", "plasmid,RepA,RepB,RepC"),
                        help="comma-separated Sanson/Brunello samples to fetch/run; use one sample for full_sample rows")
    parser.add_argument("--full", action="store_true", help="also run one full FASTQ timing per dataset unless --full-repeats is set")
    parser.add_argument("--full-only", action="store_true",
                        help="run only full FASTQ timing rows; useful for resuming missing strict-gate evidence")
    parser.add_argument("--full-repeats", type=int, default=int(os.environ.get("DOTMATCH_SOTA_FULL_REPEATS", "1")))
    parser.add_argument("--resume", action="store_true",
                        help="reuse completed rows already present in --out and checkpoint after each dataset run")
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")
    if args.full_repeats <= 0:
        raise SystemExit("--full-repeats must be positive")

    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    sizes = parse_sizes(args.read_sizes)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, str]] = read_csv(out_path) if args.resume and out_path.exists() else []
    required_tools = expected_tools(args)
    meta = machine_metadata()
    dataset_runners = selected_dataset_runners(args)

    with tempfile.TemporaryDirectory(prefix="dotmatch-crispr-sota-") as tmp_dir:
        tmp = Path(tmp_dir)
        if should_run_subsamples(args):
            for records in sizes:
                for repeat in range(1, args.repeats + 1):
                    for dataset, runner in dataset_runners:
                        requested = str(records)
                        if args.resume and run_complete(all_rows, dataset, requested, repeat, "subsample", required_tools):
                            continue
                        all_rows = without_run(all_rows, dataset, requested, repeat, "subsample")
                        tmp_out = tmp / f"{dataset}_{records}_{repeat}.csv"
                        rows = runner(records, args, tmp_out)
                        for row in rows:
                            row["repeat"] = str(repeat)
                            row["requested_records_per_sample"] = requested
                            row["run_level"] = "subsample"
                            row.setdefault("offset_mode", "best" if dataset == "mageck_yusa" else "")
                            row.update(meta)
                        all_rows.extend(rows)
                        if args.resume:
                            write_rows(out_path, all_rows)
        if should_run_full(args):
            for repeat in range(1, args.full_repeats + 1):
                for dataset, runner in dataset_runners:
                    run_level = full_run_level(dataset, args)
                    sample_id = full_sample_scope(dataset, args)
                    if args.resume and run_complete(all_rows, dataset, "full", repeat, run_level,
                                                    required_tools, sample_id=sample_id):
                        continue
                    all_rows = without_run(all_rows, dataset, "full", repeat, run_level, sample_id=sample_id)
                    tmp_out = tmp / f"{dataset}_full_{repeat}.csv"
                    rows = runner(None, args, tmp_out)
                    for row in rows:
                        row["repeat"] = str(repeat)
                        row["requested_records_per_sample"] = "full"
                        row["run_level"] = run_level
                        row.setdefault("offset_mode", "best" if dataset == "mageck_yusa" else "")
                        row.update(meta)
                    all_rows.extend(rows)
                    if args.resume:
                        write_rows(out_path, all_rows)

    if not all_rows:
        raise SystemExit("no benchmark rows were produced")
    write_rows(out_path, all_rows)
    print(out_path)


if __name__ == "__main__":
    main()
