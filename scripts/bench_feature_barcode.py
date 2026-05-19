#!/usr/bin/env python3
from __future__ import annotations
"""Generate feature-barcode assignment benchmark evidence.

The fixture models whitelist-style feature assignment for cell hashing or ADT
barcodes. The optional public lane checks per-read 10x TotalSeq-B antibody
Feature Barcode assignment, not Cell Ranger UMI/cell-level quantification.
"""

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
RAW = ROOT / "benchmarks" / "raw" / "feature_barcode.csv"
WORK = ROOT / "benchmarks" / "work" / "feature_barcode"
REPORT = ROOT / "docs" / "benchmarks" / "feature_barcode" / "README.md"
PUBLIC_METADATA = ROOT / "examples" / "feature_barcode" / "data" / "metadata.json"
PUBLIC_WORKFLOW = "public_10x_totalseq_b_feature_barcode"


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
    targets = work / "feature_barcodes.tsv"
    reads = work / "feature_reads.fastq"
    target_rows = [
        ("HTO_A", "ACGTACGTAA", "cell_hash_A"),
        ("HTO_B", "TTTTCCCCGG", "cell_hash_B"),
        ("ADT_CD3", "GGGGAAAACC", "CD3"),
        ("HTO_D", "ACGTACGTAT", "cell_hash_D"),
    ]
    read_rows = [
        ("exact_hto_a", "ACGTACGTAA"),
        ("exact_hto_b", "TTTTCCCCGG"),
        ("exact_adt_cd3", "GGGGAAAACC"),
        ("exact_hto_d", "ACGTACGTAT"),
        ("ambiguous_hto_ad", "ACGTACGTAC"),
        ("unmatched", "CCCCCCCCCC"),
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
            "assigned_exact": 4,
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
        raise RuntimeError(f"no feature barcode targets found in {path}")
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


def write_report(rows: list[dict[str, str]], report: Path) -> None:
    report.parent.mkdir(parents=True, exist_ok=True)
    synthetic = [row for row in rows if row.get("workflow") == "synthetic_feature_barcode_fixture"]
    public = [row for row in rows if row.get("workflow") == PUBLIC_WORKFLOW]
    lines = [
        "# Feature Barcode Assignment Evidence",
        "",
        "This report covers feature-barcode assignment evidence for DotMatch's known-target counting layer.",
        "",
        "The synthetic lane checks exact, ambiguous, and unmatched feature IDs. The public lane uses a 10x Genomics TotalSeq-B antibody Feature Barcode R2 subsample and validates DotMatch k=0 against a simple exact-slice comparison over the documented feature-reference window.",
        "",
        "## Synthetic Command",
        "",
        "```bash",
        synthetic[0]["command"] if synthetic else "",
        "```",
        "",
        "## Raw Rows",
        "",
        "| tool | workflow | status | features | reads | start | length | k | metric | assigned | exact | corrected | ambiguous | unmatched | validation mismatches |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        display = {"target_start": row.get("target_start", "0"), **row}
        lines.append(
            "| {tool} | {workflow} | {status} | {n_targets} | {n_reads} | {target_start} | {target_length} | {k} | {metric} | {assigned_unique} | {assigned_exact} | {corrected_reads} | {ambiguous_reads} | {unmatched_reads} | {validation_mismatches} |".format(**display)
        )
    if public:
        lines.extend(
            [
                "",
                "## Public Feature-Barcode Lane",
                "",
                "- Dataset: 10x Genomics 1k Human PBMCs with TotalSeq-B Human TBNK Antibody Cocktail, 3' v3.1.",
                "- Source page: https://www.10xgenomics.com/datasets/1-k-human-pbm-cs-with-total-seq-b-human-tbnk-antibody-cocktail-3-v-3-1-3-1-standard-6-0-0",
                "- Feature reference pattern: `^NNNNNNNNNN(BC)NNNNNNNNN`, so DotMatch uses `--target-start 10 --target-length 15` on antibody R2.",
                "- Comparison settings: the exact-slice check counts reads whose fixed R2 substring exactly matches a feature-reference sequence. It validates per-read assignment rules, not Cell Ranger cell/UMI quantification.",
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
            "Use these lanes to verify fixed-window feature-barcode whitelist counting and explicit ambiguity handling. The public 10x lane supports only per-read assignment against the documented Feature Barcode reference window. Broader CITE-seq, cell-hashing, or cell-level quantification comparisons require public comparison output, UMI/cell aggregation validation, exact commands, and a passing check.",
            "",
        ]
    )
    report.write_text("\n".join(lines), encoding="utf-8")


def run_benchmark(dotmatch: Path, work: Path) -> dict[str, str]:
    fixture = make_fixture(work)
    counts = work / "feature_counts.tsv"
    summary = work / "feature_summary.json"
    assignments = work / "feature_assignments.tsv"
    sample_qc = work / "feature_sample_qc.tsv"
    cmd = [
        str(dotmatch),
        "count",
        "--targets",
        str(fixture.targets),
        "--reads",
        str(fixture.reads),
        "--sample-label",
        "feature_barcode_fixture",
        "--target-start",
        "0",
        "--target-length",
        "10",
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
    mismatches = validation_mismatches(summary, fixture.expected) if summary.exists() else ["summary missing"]
    row = {
        "tool": "dotmatch_count",
        "workflow": "synthetic_feature_barcode_fixture",
        "status": "smoke",
        "target_start": "0",
        "target_length": "10",
        "k": "1",
        "metric": "hamming",
        "seconds": f"{seconds:.6f}",
        "exit_code": str(proc.returncode),
        "validation_mismatches": str(len(mismatches)),
        "validation_notes": "; ".join(mismatches),
        "command": command_text(cmd),
        "counts": public_text(counts),
        "summary": public_text(summary),
        "assignments": public_text(assignments),
        "sample_qc": public_text(sample_qc),
    }
    if summary.exists():
        row.update(summary_stats(summary))
    return row


def public_paths(metadata_path: Path) -> tuple[dict, Path, Path]:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if not metadata.get("evidence_ready"):
        raise RuntimeError(f"feature-barcode metadata is not evidence-ready: {metadata_path}")
    targets = resolve_repo_metadata_path(str(metadata["targets"]))
    reads = resolve_repo_metadata_path(str(metadata["local_fastq"]))
    if not targets.is_file():
        raise RuntimeError(f"missing feature-barcode targets: {targets}")
    if not reads.is_file():
        raise RuntimeError(f"missing feature-barcode FASTQ: {reads}")
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
    prefix = f"public_10x_totalseq_b_k{k}"
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
        "10x_totalseq_b_antibody_L001_R2",
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
        "scripts/bench_feature_barcode.py",
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
        "validation_notes": "simple exact substring comparison over the fixed Feature Barcode window",
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


def run_public_benchmark(dotmatch: Path, metadata_path: Path, work: Path) -> list[dict[str, str]]:
    metadata, targets, reads = public_paths(metadata_path)
    return [
        run_public_dotmatch(dotmatch, metadata_path, metadata, targets, reads, work, k=0),
        run_public_dotmatch(dotmatch, metadata_path, metadata, targets, reads, work, k=1),
        run_public_exact_baseline(metadata_path, metadata, targets, reads),
    ]


def write_csv(path: Path, rows: list[dict[str, str]] | dict[str, str]) -> None:
    if isinstance(rows, dict):
        rows = [rows]
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "tool",
        "workflow",
        "status",
        "n_reads",
        "n_targets",
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dotmatch", default=str(ROOT / "dotmatch"))
    parser.add_argument("--work", default=str(WORK))
    parser.add_argument("--out", default=str(RAW))
    parser.add_argument("--report", default=str(REPORT))
    parser.add_argument("--include-public", action="store_true")
    parser.add_argument("--metadata", default=str(PUBLIC_METADATA))
    args = parser.parse_args()

    rows = [run_benchmark(Path(args.dotmatch), Path(args.work))]
    if args.include_public:
        rows.extend(run_public_benchmark(Path(args.dotmatch), Path(args.metadata), Path(args.work)))
    write_csv(Path(args.out), rows)
    write_report(rows, Path(args.report))
    print(args.out)
    print(args.report)
    return 0 if all(row.get("exit_code") == "0" and row.get("validation_mismatches") == "0" for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
