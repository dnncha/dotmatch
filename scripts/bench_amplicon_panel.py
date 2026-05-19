#!/usr/bin/env python3
"""Generate a small amplicon/panel target-assignment smoke benchmark.

The fixture is deterministic and intended for CI, report plumbing, and claim
discipline. It is not public clinical, diagnostic, or wet-lab benchmark
evidence.
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
RAW = ROOT / "benchmarks" / "raw" / "amplicon_panel.csv"
WORK = ROOT / "benchmarks" / "work" / "amplicon_panel"
REPORT = ROOT / "docs" / "benchmarks" / "amplicon_panel" / "README.md"
PUBLIC_METADATA = ROOT / "examples" / "amplicon_panel" / "data" / "metadata.json"
PUBLIC_WORKFLOW = "public_nfcore_artic_v3_amplicon_primer"


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
    targets = work / "panel_targets.tsv"
    reads = work / "panel_reads.fastq"
    target_rows = [
        ("AMP_A", "ACGTACGTACGT", "GENE_A"),
        ("AMP_B", "TTTTCCCCAAAA", "GENE_B"),
        ("AMP_C", "GGGGAAAACCCC", "GENE_C"),
        ("AMP_D", "ACGTACGTACGA", "GENE_D"),
    ]
    read_rows = [
        ("exact_a", "ACGTACGTACGT"),
        ("exact_b", "TTTTCCCCAAAA"),
        ("exact_c", "GGGGAAAACCCC"),
        ("exact_d", "ACGTACGTACGA"),
        ("ambiguous_ad", "ACGTACGTACGG"),
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
        raise RuntimeError(f"no amplicon/panel targets found in {path}")
    return targets


def exact_prefix_hash_stats(targets: Path, reads: Path, target_length: int) -> dict[str, str]:
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
            if len(seq) >= target_length and seq[:target_length] in target_by_seq:
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
        raise RuntimeError(f"amplicon/panel metadata is not evidence-ready: {metadata_path}")
    targets = resolve_repo_metadata_path(str(metadata["targets"]))
    reads = resolve_repo_metadata_path(str(metadata["local_fastq"]))
    if not targets.is_file():
        raise RuntimeError(f"missing amplicon/panel targets: {targets}")
    if not reads.is_file():
        raise RuntimeError(f"missing amplicon/panel FASTQ: {reads}")
    return metadata, targets, reads


def write_report(rows, report: Path) -> None:
    if isinstance(rows, dict):
        rows = [rows]
    report.parent.mkdir(parents=True, exist_ok=True)
    synthetic = [row for row in rows if row.get("workflow") == "synthetic_amplicon_panel_fixture"]
    public = [row for row in rows if row.get("workflow") == PUBLIC_WORKFLOW]
    lines = [
        "# Amplicon/Panel Assignment Evidence",
        "",
        "This report covers panel-style target assignment evidence for DotMatch's known-target counting layer.",
        "",
        "The synthetic lane checks exact, ambiguous, and unmatched fixed-window target assignment. The public lane uses nf-core viralrecon Illumina ARTIC V3 amplicon sample R1 and validates DotMatch k=0 against a simple exact-prefix comparison over the selected full-primer length group.",
        "",
        "Current status: public primer-start assignment evidence only. This is not amplicon consensus, variant calling, primer trimming, or clinical validation evidence.",
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
        display = {"target_start": row.get("target_start", "0"), **row}
        lines.append(
            "| {tool} | {workflow} | {status} | {n_targets} | {n_reads} | {target_start} | {target_length} | {k} | {metric} | {assigned_unique} | {assigned_exact} | {corrected_reads} | {ambiguous_reads} | {unmatched_reads} | {validation_mismatches} |".format(**display)
        )
    if public:
        lines.extend(
            [
                "",
                "## Public Amplicon/Panel Lane",
                "",
                "- Dataset: nf-core/test-datasets viralrecon Illumina amplicon `sample1_R1.fastq.gz`.",
                "- Primer panel: ARTIC V3 SARS-CoV-2 primer FASTA from the same nf-core/test-datasets branch.",
                "- Fixture rules: the fetcher selects the full-primer length group with the most exact R1 prefix assignments, then counts fixed-position primer starts.",
                "- Comparison settings: the exact-prefix check counts reads whose R1 prefix exactly matches the selected primer sequence. It validates per-read known-primer assignment rules, not consensus generation, primer trimming, variant calling, or clinical interpretation.",
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
            "Use these lanes to verify fixed-window known-target panel assignment, explicit ambiguity handling, and narrow public ARTIC primer-start per-read assignment. Broader amplicon/panel benchmarks require public full-assay comparison settings, consensus or variant-call validation where relevant, exact commands, validation files, and a passing check.",
            "",
        ]
    )
    report.write_text("\n".join(lines), encoding="utf-8")


def run_benchmark(dotmatch: Path, work: Path) -> dict[str, str]:
    fixture = make_fixture(work)
    counts = work / "panel_counts.tsv"
    summary = work / "panel_summary.json"
    assignments = work / "panel_assignments.tsv"
    sample_qc = work / "panel_sample_qc.tsv"
    cmd = [
        str(dotmatch),
        "count",
        "--targets",
        str(fixture.targets),
        "--reads",
        str(fixture.reads),
        "--sample-label",
        "amplicon_panel_fixture",
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
    mismatches = validation_mismatches(summary, fixture.expected) if summary.exists() else ["summary missing"]
    row = {
        "tool": "dotmatch_count",
        "workflow": "synthetic_amplicon_panel_fixture",
        "status": "smoke",
        "target_length": "12",
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
    prefix = f"public_nfcore_amplicon_k{k}"
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
        "nfcore_viralrecon_artic_v3_sample1_R1",
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
    target_length = int(metadata["target_length"])
    stats = exact_prefix_hash_stats(targets, reads, target_length)
    cmd = [
        "python3",
        "scripts/bench_amplicon_panel.py",
        "--include-public",
        "--metadata",
        str(metadata_path),
    ]
    return {
        "tool": "exact_prefix_hash",
        "workflow": PUBLIC_WORKFLOW,
        "status": "supported",
        "target_start": str(metadata["target_start"]),
        "target_length": str(target_length),
        "k": "0",
        "metric": "exact",
        "seconds": "0.000000",
        "exit_code": "0",
        "validation_mismatches": "0",
        "validation_notes": "simple exact prefix comparison over the selected ARTIC primer start window",
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
