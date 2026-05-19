#!/usr/bin/env python3
"""Validate public CRISPR DotMatch assignments against Edlib."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "examples" / "crispr_guides" / "data"
OUT = ROOT / "examples" / "crispr_guides" / "output"
RAW = ROOT / "benchmarks" / "raw" / "public_crispr_edlib_validation.csv"


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


def run_assignments(reads: Path, sample: str, sample_limit: int, assignments: Path) -> None:
    cmd = [
        str(ROOT / "dotmatch"),
        "count",
        "--targets", str(DATA / "yusa_library.csv"),
        "--reads", str(reads),
        "--sample-label", sample,
        "--target-start", "23",
        "--target-length", "19",
        "--k", "1",
        "--metric", "levenshtein",
        "--indel-window", "1",
        "--format", "mageck",
        "--out", str(OUT / f"validation.{sample}.counts.tsv"),
        "--assignments", str(assignments),
        "--ambiguous", "report",
        "--summary", str(OUT / f"validation.{sample}.summary.json"),
    ]
    env = os.environ.copy()
    env["DOTMATCH_MAX_READS"] = str(sample_limit)
    subprocess.run(cmd, cwd=ROOT, env=env, check=True, stdout=subprocess.DEVNULL)


def strata(assignments: Path) -> dict[str, int]:
    counts = {"exact": 0, "corrected": 0, "ambiguous": 0, "unmatched": 0, "contains_n": 0}
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
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records-per-sample", type=int, default=int(os.environ.get("DOTMATCH_PUBLIC_VALIDATION_RECORDS", "1000")))
    parser.add_argument("--sample", type=int, default=int(os.environ.get("DOTMATCH_PUBLIC_VALIDATION_SAMPLE", "1000")))
    parser.add_argument("--out", default=str(RAW))
    args = parser.parse_args()

    subprocess.run(["python3", str(ROOT / "scripts" / "fetch_mageck_demo.py"), "--out", str(DATA),
                    "--subsample", str(args.records_per_sample)], cwd=ROOT, check=True)
    subprocess.run(["make", "dotmatch", "edlib-tools"], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)

    rows: list[dict[str, str]] = []
    samples = [("plasmid", DATA / "ERR376998.fastq.gz"), ("ESC1", DATA / "ERR376999.fastq.gz")]
    with tempfile.TemporaryDirectory(prefix="dotmatch-edlib-validation-") as tmp_dir:
        tmp = Path(tmp_dir)
        for label, reads in samples:
            validation_reads = tmp / f"{label}.fastq.gz"
            write_fastq_prefix(reads, validation_reads, args.records_per_sample)
            assignments = tmp / f"{label}.assignments.tsv"
            run_assignments(validation_reads, label, args.records_per_sample, assignments)
            observed_strata = strata(assignments)
            result = run_json([
                str(ROOT / "dotmatch"),
                "validate",
                "--targets", str(DATA / "yusa_library.csv"),
                "--reads", str(validation_reads),
                "--target-start", "23",
                "--target-length", "19",
                "--k", "1",
                "--indel-window", "1",
                "--oracle", "edlib",
                "--sample", str(args.sample),
            ])
            row = {
                "dataset": "mageck_yusa",
                "sample": label,
                "oracle": str(result.get("oracle", "edlib_native")),
                "checked_reads": str(result.get("checked_reads", "")),
                "mismatches": str(result.get("mismatches", "")),
                "k": str(result.get("k", "1")),
                "target_start": str(result.get("target_start", "23")),
                "target_length": str(result.get("target_length", "19")),
                "indel_window": str(result.get("indel_window", "1")),
            }
            row.update({f"stratum_{k}": str(v) for k, v in observed_strata.items()})
            rows.append(row)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    if any(row["mismatches"] != "0" for row in rows):
        raise SystemExit("Edlib validation mismatches detected")
    print(out_path)


if __name__ == "__main__":
    main()
