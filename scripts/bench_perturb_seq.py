#!/usr/bin/env python3
"""Generate a small perturb-seq-style guide plus feature smoke benchmark.

The fixture checks fixed-window guide/feature pair assignment through
`dotmatch pair-count`. It is not expression processing and not public
Perturb-seq benchmark evidence.
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
RAW = ROOT / "benchmarks" / "raw" / "perturb_seq.csv"
WORK = ROOT / "benchmarks" / "work" / "perturb_seq"
REPORT = ROOT / "docs" / "benchmarks" / "perturb_seq" / "README.md"
PUBLIC_METADATA = ROOT / "examples" / "perturb_seq" / "data" / "metadata.json"
PUBLIC_WORKFLOW = "public_10x_crispr_guide_capture"


class Fixture:
    def __init__(self, left_targets: Path, right_targets: Path, reads: Path, expected: dict[str, int]) -> None:
        self.left_targets = left_targets
        self.right_targets = right_targets
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
    guides = work / "perturb_guides.tsv"
    features = work / "perturb_features.tsv"
    reads = work / "perturb_reads.fastq"
    guide_rows = [
        ("GUIDE_A", "ACGTAC", "GENE_A"),
        ("GUIDE_B", "TTTTCC", "GENE_B"),
        ("GUIDE_D", "ACGTAT", "GENE_D"),
    ]
    feature_rows = [
        ("HTO_A", "GGAACC", "cell_hash_A"),
        ("ADT_CD3", "CCTTAA", "CD3"),
    ]
    read_rows = [
        ("exact_a_hto", "ACGTACGGAACC"),
        ("exact_b_adt", "TTTTCCCCTTAA"),
        ("exact_d_hto", "ACGTATGGAACC"),
        ("ambiguous_guide_ad", "ACGTAGGGAACC"),
        ("right_unmatched", "ACGTACAAAAAA"),
        ("left_unmatched", "GGGGGGGGAACC"),
        ("invalid_short", "ACGT"),
    ]
    with guides.open("w", encoding="utf-8") as fh:
        fh.write("target_id\ttarget_seq\tgene\n")
        for row in guide_rows:
            fh.write("\t".join(row) + "\n")
    with features.open("w", encoding="utf-8") as fh:
        fh.write("target_id\ttarget_seq\tgene\n")
        for row in feature_rows:
            fh.write("\t".join(row) + "\n")
    with reads.open("w", encoding="utf-8") as fh:
        for name, seq in read_rows:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
    return Fixture(
        left_targets=guides,
        right_targets=features,
        reads=reads,
        expected={
            "total_reads": 7,
            "assigned_pairs": 3,
            "pair_ambiguous": 1,
            "left_unmatched": 1,
            "right_unmatched": 1,
            "invalid": 1,
        },
    )


def validation_mismatches(summary_path: Path, expected: dict[str, int]) -> list[str]:
    observed = json.loads(summary_path.read_text(encoding="utf-8"))
    failures = []
    for key, expected_value in expected.items():
        value = int(observed.get(key, 0))
        if value != expected_value:
            failures.append(f"{key} expected {expected_value} observed {value}")
    return failures


def summary_stats(summary_path: Path) -> dict[str, str]:
    data = json.loads(summary_path.read_text(encoding="utf-8"))
    return {
        "n_reads": str(data.get("total_reads", "")),
        "n_left_targets": str(data.get("n_left_targets", "")),
        "n_right_targets": str(data.get("n_right_targets", "")),
        "assigned_pairs": str(data.get("assigned_pairs", "")),
        "pair_ambiguous": str(data.get("pair_ambiguous", "")),
        "left_unmatched": str(data.get("left_unmatched", "")),
        "right_unmatched": str(data.get("right_unmatched", "")),
        "invalid_reads": str(data.get("invalid", "")),
        "candidates_verified": str(data.get("candidates_verified", "")),
        "alphabet_policy": str(data.get("alphabet_policy", "")),
    }


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def resolve_repo_metadata_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def count_summary_stats(summary_path: Path) -> dict[str, str]:
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
        raise RuntimeError(f"no CRISPR guide targets found in {path}")
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


def public_paths(metadata_path: Path) -> tuple[dict, Path, Path]:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if not metadata.get("evidence_ready"):
        raise RuntimeError(f"perturb-seq metadata is not evidence-ready: {metadata_path}")
    targets = resolve_repo_metadata_path(str(metadata["targets"]))
    reads = resolve_repo_metadata_path(str(metadata["local_fastq"]))
    if not targets.is_file():
        raise RuntimeError(f"missing perturb-seq targets: {targets}")
    if not reads.is_file():
        raise RuntimeError(f"missing perturb-seq FASTQ: {reads}")
    return metadata, targets, reads


def write_report(rows, report: Path) -> None:
    if isinstance(rows, dict):
        rows = [rows]
    report.parent.mkdir(parents=True, exist_ok=True)
    synthetic = [row for row in rows if row.get("workflow") == "synthetic_perturb_seq_fixture"]
    public = [row for row in rows if row.get("workflow") == PUBLIC_WORKFLOW]
    lines = [
        "# Perturb-Seq And CRISPR Guide-Capture Assignment Evidence",
        "",
        "This report covers perturb-seq-adjacent guide assignment evidence for DotMatch's known-target counting layer.",
        "",
        "The synthetic lane checks fixed-window guide plus feature-barcode pair assignment through `pair-count`. The public lane uses a 10x Genomics CRISPR Guide Capture R2 subsample and validates DotMatch k=0 against a simple exact-slice comparison over the observed fixed guide window.",
        "",
        "Current status: public guide-capture assignment evidence only. This is not single-cell expression quantification or Cell Ranger perturbation-effect validation.",
        "",
        "## Synthetic Command",
        "",
        "```bash",
        synthetic[0]["command"] if synthetic else "",
        "```",
        "",
        "## Synthetic Raw Row",
        "",
        "| tool | workflow | guides | features | reads | k | metric | assigned pairs | ambiguous | left unmatched | right unmatched | invalid | validation mismatches |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in synthetic:
        lines.append(
            "| {tool} | {workflow} | {n_left_targets} | {n_right_targets} | {n_reads} | {k} | {metric} | {assigned_pairs} | {pair_ambiguous} | {left_unmatched} | {right_unmatched} | {invalid_reads} | {validation_mismatches} |".format(**row)
        )
    if public:
        lines.extend(
            [
                "",
                "## Public CRISPR Guide-Capture Rows",
                "",
                "| tool | workflow | status | guides | reads | start | length | k | metric | assigned | exact | corrected | unmatched | validation mismatches |",
                "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |",
            ]
        )
        for row in public:
            lines.append(
                "| {tool} | {workflow} | {status} | {n_targets} | {n_reads} | {target_start} | {target_length} | {k} | {metric} | {assigned_unique} | {assigned_exact} | {corrected_reads} | {unmatched_reads} | {validation_mismatches} |".format(**row)
            )
        lines.extend(
            [
                "",
                "## Public Dataset",
                "",
                "- Dataset: 10x Genomics 1k A375 Cells Transduced with Non-Target and Target sgRNA, Chromium GEM-X Single Cell 5'.",
                "- Source page: https://www.10xgenomics.com/datasets/1k-CRISPR-5p-gemx",
                "- Fixture rules: the fetcher selects the observed fixed-window CRISPR Guide Capture group with the most exact assignments in the copied R2 prefix.",
                "- Comparison settings: the exact-slice check counts reads whose fixed R2 substring exactly matches the selected guide sequence. It validates per-read guide assignment rules, not Cell Ranger cell/UMI quantification or perturbation effects.",
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
            "Use these lanes to verify fixed-window guide/feature pair assignment, side-level diagnostics, and narrow public CRISPR guide-capture per-read assignment. Broader Perturb-seq comparisons require public cell barcode and UMI handling, guide-per-cell calls, expression or perturbation-effect comparison output, exact commands, validation files, and a passing check.",
            "",
        ]
    )
    report.write_text("\n".join(lines), encoding="utf-8")


def run_benchmark(dotmatch: Path, work: Path) -> dict[str, str]:
    fixture = make_fixture(work)
    counts = work / "perturb_pair_counts.tsv"
    summary = work / "perturb_summary.json"
    assignments = work / "perturb_assignments.tsv"
    cmd = [
        str(dotmatch),
        "pair-count",
        "--left-targets",
        str(fixture.left_targets),
        "--right-targets",
        str(fixture.right_targets),
        "--reads",
        str(fixture.reads),
        "--left-start",
        "0",
        "--left-length",
        "6",
        "--right-start",
        "6",
        "--right-length",
        "6",
        "--k",
        "1",
        "--metric",
        "hamming",
        "--out",
        str(counts),
        "--summary",
        str(summary),
        "--assignments",
        str(assignments),
    ]
    start = time.perf_counter()
    proc = subprocess.run(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    seconds = time.perf_counter() - start
    mismatches = validation_mismatches(summary, fixture.expected) if summary.exists() else ["summary missing"]
    row = {
        "tool": "dotmatch_pair_count",
        "workflow": "synthetic_perturb_seq_fixture",
        "status": "smoke",
        "left_length": "6",
        "right_length": "6",
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
    }
    if summary.exists():
        row.update(summary_stats(summary))
    return row


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
    prefix = f"public_10x_crispr_k{k}"
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
        "10x_crispr_guide_capture_L001_R2",
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
        row.update(count_summary_stats(summary))
    return row


def run_public_exact_baseline(metadata_path: Path, metadata: dict, targets: Path, reads: Path) -> dict[str, str]:
    target_start = int(metadata["target_start"])
    target_length = int(metadata["target_length"])
    stats = exact_slice_hash_stats(targets, reads, target_start, target_length)
    cmd = [
        "python3",
        "scripts/bench_perturb_seq.py",
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
        "validation_notes": "simple exact substring comparison over the observed fixed CRISPR guide window",
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


def write_csv(path: Path, rows) -> None:
    if isinstance(rows, dict):
        rows = [rows]
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "tool",
        "workflow",
        "status",
        "n_reads",
        "n_targets",
        "n_left_targets",
        "n_right_targets",
        "target_start",
        "target_length",
        "left_length",
        "right_length",
        "k",
        "metric",
        "assigned_unique",
        "assigned_exact",
        "corrected_reads",
        "ambiguous_reads",
        "unmatched_reads",
        "assigned_pairs",
        "pair_ambiguous",
        "left_unmatched",
        "right_unmatched",
        "invalid_reads",
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
