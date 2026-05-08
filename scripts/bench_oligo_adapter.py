#!/usr/bin/env python3
"""Generate oligo/adapter fixed-window assignment evidence.

This benchmark exercises DotMatch's known-target counting layer on short
adapter-like oligos and an optional public adapter-prefix fixture. It is not
adapter trimming, read merging, or UMI grouping evidence.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "oligo_adapter.csv"
WORK = ROOT / "benchmarks" / "work" / "oligo_adapter"
REPORT = ROOT / "docs" / "benchmarks" / "oligo_adapter" / "README.md"
WORKFLOW = "synthetic_oligo_adapter_fixture"
PUBLIC_METADATA = ROOT / "examples" / "oligo_adapter" / "data" / "metadata.json"
PUBLIC_WORKFLOW = "public_fast_adapter_truseq_r1"


class Fixture:
    def __init__(self, targets: Path, reads: Path, expected: dict[str, int]) -> None:
        self.targets = targets
        self.reads = reads
        self.expected = expected


def public_text(value) -> str:
    text = str(value)
    root = str(ROOT)
    text = text.replace(root + os.sep, "")
    if text == root:
        text = "."
    var_folders = "/" + "var/folders/"
    private_tmp = "/" + "private/tmp/"
    text = re.sub(re.escape(var_folders) + r"[^,\s\"]*/([^/,\s\"]+)", r"<tmp>/\1", text)
    text = re.sub(re.escape(private_tmp) + r"[^,\s\"]*/([^/,\s\"]+)", r"<tmp>/\1", text)
    return text


def command_text(cmd: list[str]) -> str:
    return " ".join(public_text(arg) for arg in cmd)


def make_fixture(work: Path) -> Fixture:
    work.mkdir(parents=True, exist_ok=True)
    targets = work / "adapter_oligos.tsv"
    reads = work / "adapter_reads.fastq"
    target_rows = [
        ("ADAPTER_A", "ACGTACGTACGT", "adapter_A"),
        ("ADAPTER_B", "TTTTCCCCAAAA", "adapter_B"),
        ("ADAPTER_D", "ACGTACGTACGA", "adapter_D"),
    ]
    read_rows = [
        ("exact_adapter_a", "ACGTACGTACGT"),
        ("exact_adapter_b", "TTTTCCCCAAAA"),
        ("exact_adapter_d", "ACGTACGTACGA"),
        ("corrected_adapter_b", "TTTTCCCCAAAT"),
        ("ambiguous_adapter_ad", "ACGTACGTACGC"),
        ("unmatched", "CCCCCCCCCCCC"),
    ]
    with targets.open("w", encoding="utf-8") as fh:
        fh.write("target_id\ttarget_seq\tgene\n")
        for row in target_rows:
            fh.write("\t".join(row) + "\n")
    with reads.open("w", encoding="utf-8") as fh:
        for name, seq in read_rows:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
    return Fixture(
        targets=targets,
        reads=reads,
        expected={
            "total_reads": 6,
            "assigned_unique": 4,
            "assigned_exact": 3,
            "assigned_corrected": 1,
            "ambiguous": 1,
            "unmatched": 1,
        },
    )


def validation_mismatches(summary_path: Path, expected: dict[str, int]) -> list[str]:
    observed = json.loads(summary_path.read_text(encoding="utf-8"))
    if isinstance(observed.get("samples"), list) and observed["samples"]:
        observed = observed["samples"][0]
    failures = []
    for key, expected_value in expected.items():
        value = int(observed.get(key, 0))
        if value != expected_value:
            failures.append(f"{key} expected {expected_value} observed {value}")
    return failures


def summary_stats(summary_path: Path) -> dict[str, str]:
    data = json.loads(summary_path.read_text(encoding="utf-8"))
    sample = data["samples"][0] if isinstance(data.get("samples"), list) and data["samples"] else data
    unmatched = int(sample.get("unmatched", 0)) + int(sample.get("invalid", 0))
    return {
        "n_reads": str(sample.get("total_reads", "")),
        "n_targets": str(data.get("n_targets", "")),
        "assigned_unique": str(sample.get("assigned_unique", "")),
        "assigned_exact": str(sample.get("assigned_exact", "")),
        "corrected_reads": str(sample.get("assigned_corrected", "")),
        "ambiguous_reads": str(sample.get("ambiguous", "")),
        "unmatched_reads": str(unmatched),
        "candidates_verified": str(sample.get("candidates_verified", "")),
        "alphabet_policy": str(data.get("alphabet_policy", "")),
    }


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def resolve_repo_metadata_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def read_target_sequences(path: Path) -> dict[str, str]:
    targets: dict[str, str] = {}
    with open_text(path) as fh:
        first_data = True
        for raw in fh:
            line = raw.rstrip("\n\r")
            if not line:
                continue
            cols = line.split("\t")
            if first_data and {"target_id", "target_seq"} & {col.strip().lower() for col in cols}:
                first_data = False
                continue
            first_data = False
            if len(cols) >= 2:
                targets[cols[1].strip().upper()] = cols[0].strip()
    if not targets:
        raise RuntimeError(f"no adapter targets found in {path}")
    return targets


def exact_slice_hash_stats(targets: Path, reads: Path, target_start: int, target_length: int) -> dict[str, str]:
    target_by_seq = read_target_sequences(targets)
    total = 0
    assigned = 0
    with open_text(reads) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip().upper()
            plus = fh.readline()
            qual = fh.readline()
            if not seq or not plus or not qual:
                raise RuntimeError(f"truncated FASTQ record in {reads}")
            if target_start + target_length <= len(seq):
                observed = seq[target_start:target_start + target_length]
                if observed in target_by_seq:
                    assigned += 1
            total += 1
    return {
        "n_reads": str(total),
        "n_targets": str(len(target_by_seq)),
        "assigned_unique": str(assigned),
        "assigned_exact": str(assigned),
        "corrected_reads": "0",
        "ambiguous_reads": "0",
        "unmatched_reads": str(total - assigned),
    }


def run_benchmark(dotmatch: Path, work: Path) -> list[dict[str, str]]:
    fixture = make_fixture(work)
    counts = work / "adapter_counts.tsv"
    summary = work / "adapter_summary.json"
    assignments = work / "adapter_assignments.tsv"
    sample_qc = work / "adapter_sample_qc.tsv"
    cmd = [
        str(dotmatch),
        "count",
        "--targets",
        str(fixture.targets),
        "--reads",
        str(fixture.reads),
        "--sample-label",
        "oligo_adapter_fixture",
        "--target-start",
        "0",
        "--target-length",
        "12",
        "--k",
        "1",
        "--metric",
        "hamming",
        "--format",
        "dotmatch",
        "--out",
        str(counts),
        "--summary",
        str(summary),
        "--assignments",
        str(assignments),
        "--ambiguous",
        "report",
        "--sample-qc",
        str(sample_qc),
    ]
    start = time.perf_counter()
    proc = subprocess.run(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    seconds = time.perf_counter() - start
    stats = summary_stats(summary) if summary.exists() else {
        "n_reads": "0",
        "n_targets": "0",
        "assigned_unique": "0",
        "assigned_exact": "0",
        "corrected_reads": "0",
        "ambiguous_reads": "0",
        "unmatched_reads": "0",
        "candidates_verified": "0",
        "alphabet_policy": "",
    }
    mismatches = validation_mismatches(summary, fixture.expected) if summary.exists() else ["summary missing"]
    return [
        {
            "tool": "dotmatch_count",
            "workflow": WORKFLOW,
            "status": "smoke",
            "seconds": f"{seconds:.6f}",
            "exit_code": str(proc.returncode),
            "target_start": "0",
            "target_length": "12",
            "k": "1",
            "metric": "hamming",
            "validation_mismatches": str(len(mismatches)),
            "validation_notes": "; ".join(mismatches),
            "command": command_text(cmd),
            **stats,
        }
    ]


def public_paths(metadata_path: Path) -> tuple[dict, Path, Path]:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if not metadata.get("evidence_ready"):
        raise RuntimeError(f"oligo/adapter metadata is not evidence-ready: {metadata_path}")
    targets = resolve_repo_metadata_path(str(metadata["targets"]))
    reads = resolve_repo_metadata_path(str(metadata["local_fastq"]))
    if not targets.is_file():
        raise RuntimeError(f"missing oligo/adapter targets: {targets}")
    if not reads.is_file():
        raise RuntimeError(f"missing oligo/adapter FASTQ: {reads}")
    return metadata, targets, reads


def run_public_dotmatch(
    dotmatch: Path,
    metadata_path: Path,
    metadata: dict,
    targets: Path,
    reads: Path,
    work: Path,
    k: int,
) -> dict[str, str]:
    work.mkdir(parents=True, exist_ok=True)
    target_start = int(metadata["target_start"])
    target_length = int(metadata["target_length"])
    prefix = f"public_fast_adapter_truseq_r1_k{k}"
    counts = work / f"{prefix}_counts.tsv"
    summary = work / f"{prefix}_summary.json"
    assignments = work / f"{prefix}_assignments.tsv"
    sample_qc = work / f"{prefix}_sample_qc.tsv"
    cmd = [
        str(dotmatch),
        "count",
        "--targets",
        str(targets),
        "--reads",
        str(reads),
        "--sample-label",
        "fast_adapter_truseq_r1",
        "--target-start",
        str(target_start),
        "--target-length",
        str(target_length),
        "--k",
        str(k),
        "--metric",
        "hamming",
        "--format",
        "dotmatch",
        "--out",
        str(counts),
        "--summary",
        str(summary),
        "--assignments",
        str(assignments),
        "--ambiguous",
        "report",
        "--sample-qc",
        str(sample_qc),
    ]
    start = time.perf_counter()
    proc = subprocess.run(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    seconds = time.perf_counter() - start
    ok = proc.returncode == 0 and summary.exists()
    row = {
        "tool": "dotmatch_count",
        "workflow": PUBLIC_WORKFLOW,
        "status": "supported",
        "target_start": str(target_start),
        "target_length": str(target_length),
        "k": str(k),
        "metric": "hamming",
        "seconds": f"{seconds:.6f}",
        "exit_code": str(proc.returncode),
        "validation_mismatches": "0" if ok else "1",
        "validation_notes": "" if ok else "summary missing or command failed",
        "command": command_text(cmd),
        "counts": public_text(counts),
        "summary": public_text(summary),
        "assignments": public_text(assignments),
        "sample_qc": public_text(sample_qc),
        "metadata": public_text(metadata_path),
    }
    if summary.exists():
        row.update(summary_stats(summary))
    return row


def run_public_exact_baseline(metadata_path: Path, metadata: dict, targets: Path, reads: Path) -> dict[str, str]:
    target_start = int(metadata["target_start"])
    target_length = int(metadata["target_length"])
    stats = exact_slice_hash_stats(targets, reads, target_start, target_length)
    cmd = [
        "python3",
        "scripts/bench_oligo_adapter.py",
        "--include-public",
        "--metadata",
        str(metadata_path),
    ]
    return {
        "tool": "exact_slice_hash",
        "workflow": PUBLIC_WORKFLOW,
        "status": "supported",
        "target_start": str(target_start),
        "target_length": str(target_length),
        "k": "0",
        "metric": "exact",
        "seconds": "0.000000",
        "exit_code": "0",
        "validation_mismatches": "0",
        "validation_notes": "transparent exact substring baseline over the fixed public adapter-prefix window",
        "command": command_text(cmd),
        "counts": "",
        "summary": "",
        "assignments": "",
        "sample_qc": "",
        "metadata": public_text(metadata_path),
        "alphabet_policy": "",
        "candidates_verified": "",
        **stats,
    }


def _set_public_validation(dotmatch_k0: dict[str, str], dotmatch_k1: dict[str, str], baseline: dict[str, str]) -> None:
    notes: list[str] = []
    if dotmatch_k0.get("assigned_unique") != baseline.get("assigned_unique"):
        notes.append("DotMatch k=0 assigned_unique differs from exact-slice baseline")
    if dotmatch_k0.get("assigned_exact") != baseline.get("assigned_exact"):
        notes.append("DotMatch k=0 assigned_exact differs from exact-slice baseline")
    try:
        if int(dotmatch_k1.get("assigned_unique") or 0) < int(dotmatch_k0.get("assigned_unique") or 0):
            notes.append("DotMatch k=1 assigned fewer reads than k=0")
    except ValueError:
        notes.append("could not parse public assigned_unique counts")
    if notes:
        dotmatch_k0["validation_mismatches"] = str(int(dotmatch_k0.get("validation_mismatches") or 0) + len(notes))
        dotmatch_k0["validation_notes"] = "; ".join(filter(None, [dotmatch_k0.get("validation_notes", ""), *notes]))


def run_public_benchmark(dotmatch: Path, metadata_path: Path, work: Path) -> list[dict[str, str]]:
    metadata, targets, reads = public_paths(metadata_path)
    k0 = run_public_dotmatch(dotmatch, metadata_path, metadata, targets, reads, work, k=0)
    k1 = run_public_dotmatch(dotmatch, metadata_path, metadata, targets, reads, work, k=1)
    baseline = run_public_exact_baseline(metadata_path, metadata, targets, reads)
    _set_public_validation(k0, k1, baseline)
    return [k0, k1, baseline]


def write_csv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "tool",
        "workflow",
        "status",
        "n_targets",
        "n_reads",
        "target_start",
        "target_length",
        "k",
        "metric",
        "assigned_unique",
        "assigned_exact",
        "corrected_reads",
        "ambiguous_reads",
        "unmatched_reads",
        "candidates_verified",
        "validation_mismatches",
        "validation_notes",
        "seconds",
        "exit_code",
        "alphabet_policy",
        "command",
        "counts",
        "summary",
        "assignments",
        "sample_qc",
        "metadata",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_report(rows: list[dict[str, str]], report: Path) -> None:
    report.parent.mkdir(parents=True, exist_ok=True)
    synthetic = [row for row in rows if row.get("workflow") == WORKFLOW]
    public = [row for row in rows if row.get("workflow") == PUBLIC_WORKFLOW]
    lines = [
        "# Oligo/Adapter Assignment Evidence",
        "",
        "This report covers fixed-window assignment of short adapter-like oligos with DotMatch's known-target counting layer.",
        "",
        "The synthetic lane checks exact, one-substitution, ambiguous, and unmatched adapter-like oligos. The public lane uses the fast-adapter-trimming TruSeq R1 fixture and validates DotMatch k=0 against a transparent exact-slice hash baseline over the documented fixed window.",
        "",
        "## Synthetic Command",
        "",
        "```bash",
        synthetic[0]["command"] if synthetic else "",
        "```",
        "",
        "## Raw Rows",
        "",
        "| tool | workflow | status | targets | reads | start | length | k | metric | assigned | exact | corrected | ambiguous | unmatched | validation mismatches |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| {tool} | {workflow} | {status} | {n_targets} | {n_reads} | {target_start} | {target_length} | {k} | {metric} | {assigned_unique} | {assigned_exact} | {corrected_reads} | {ambiguous_reads} | {unmatched_reads} | {validation_mismatches} |".format(**row)
        )
    if public:
        lines.extend(
            [
                "",
                "## Public Adapter-Prefix Lane",
                "",
                "- Dataset: fast-adapter-trimming `788707_20180313_S_R1.small.fastq.gz` with `adapters/truseq.fa` target prefixes.",
                "- Source repository: https://github.com/linsalrob/fast-adapter-trimming",
                "- Source license: MIT, as reported by the upstream GitHub repository metadata.",
                "- Fixed window: DotMatch uses `--target-start 229 --target-length 20` on R1.",
                "- Comparator semantics: the exact-slice hash baseline counts reads whose fixed R1 substring exactly matches a deduplicated TruSeq adapter-prefix target. It validates fixed-window assignment semantics, not trimming correctness.",
                "",
                "## Public Commands",
                "",
            ]
        )
        for row in public:
            lines.extend(["```bash", row["command"], "```", ""])
    lines.extend(
        [
            "",
            "## Evidence Boundary",
            "",
            "Use these lanes only to verify fixed-window known-oligo/adapter assignment, one-substitution rescue, and explicit ambiguous/unmatched diagnostics. Run `make oligo-adapter-smoke-gate` for smoke evidence and `make oligo-adapter-public-gate` for the public lane. The public lane supports adapter-prefix assignment wording for the checked R1 window only; it is not adapter trimming evidence. Primer removal, UMI grouping, read merging, or production workflow claims require separate comparator semantics, raw artifacts, validation, and a passing gate.",
            "",
        ]
    )
    report.write_text("\n".join(lines), encoding="utf-8")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dotmatch", default=str(ROOT / "dotmatch"))
    parser.add_argument("--out", default=str(RAW))
    parser.add_argument("--report", default=str(REPORT))
    parser.add_argument("--work", default=str(WORK))
    parser.add_argument("--include-public", action="store_true")
    parser.add_argument("--metadata", default=str(PUBLIC_METADATA))
    args = parser.parse_args(argv)

    rows = run_benchmark(Path(args.dotmatch), Path(args.work))
    if args.include_public:
        rows.extend(run_public_benchmark(Path(args.dotmatch), Path(args.metadata), Path(args.work)))
    write_csv(rows, Path(args.out))
    write_report(rows, Path(args.report))
    print(args.out)
    print(args.report)
    return 0 if all(row.get("exit_code") == "0" and row.get("validation_mismatches") == "0" for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
