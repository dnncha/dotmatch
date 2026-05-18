from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import math
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence, TextIO

from . import __version__
from .assayspec import command_assay
from .core import (
    MATCH_AMBIGUOUS,
    MATCH_INVALID,
    MATCH_NONE,
    MATCH_UNIQUE,
    Matcher,
    MatchResult,
    assign,
    distance,
)
from .native import find_native_cli, run_native_cli


DNA = "ACGT"
CRISPR_QC_THRESHOLDS = {
    "assignment_rate_min": 0.80,
    "ambiguous_rate_max": 0.05,
    "no_match_rate_max": 0.15,
    "invalid_rate_max": 0.02,
    "coverage_fraction_min": 0.90,
    "zero_count_fraction_max": 0.10,
    "gini_index_max": 0.50,
    "top_1pct_fraction_max": 0.30,
    "pairwise_sample_pearson_min": 0.80,
}


@dataclass(frozen=True)
class Target:
    target_id: str
    seq: str
    gene: str = ""


@dataclass(frozen=True)
class ReadRecord:
    read_id: str
    seq: str
    qual: str


@dataclass(frozen=True)
class BarcodeCandidate:
    start: int
    length: int
    score: float
    sampled_reads: int
    valid_reads: int
    unique: int
    exact: int
    ambiguous: int
    no_match: int
    invalid: int
    assignment_rate: float
    exact_rate: float
    ambiguous_rate: float
    no_match_rate: float
    invalid_rate: float


def _open_text(path: str | Path, mode: str = "rt") -> TextIO:
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode, encoding="utf-8")


def _read_targets(path: str | Path) -> list[Target]:
    targets: list[Target] = []
    with _open_text(path) as fh:
        first_data = True
        for raw in fh:
            line = raw.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if first_data and _looks_like_header(cols):
                first_data = False
                continue
            first_data = False
            if len(cols) == 1:
                seq = cols[0].strip().upper()
                target_id = f"target_{len(targets)}"
                gene = ""
            else:
                target_id = cols[0].strip()
                seq = cols[1].strip().upper()
                gene = cols[2].strip() if len(cols) > 2 else ""
            if not seq:
                raise ValueError(f"empty target sequence in {path}")
            targets.append(Target(target_id=target_id, seq=seq, gene=gene))
    if not targets:
        raise ValueError(f"no targets found in {path}")
    return targets


def _looks_like_header(cols: Sequence[str]) -> bool:
    normalized = {c.strip().lower() for c in cols[:3]}
    return bool(normalized & {"target_id", "guide_id", "barcode_id"}) and bool(
        normalized & {"target_seq", "guide_seq", "barcode_seq", "sequence", "seq"}
    )


def _iter_fastq(path: str | Path) -> Iterator[ReadRecord]:
    with _open_text(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not seq or not plus or not qual:
                raise ValueError("truncated FASTQ record")
            header = header.rstrip("\n\r")
            seq = seq.rstrip("\n\r").upper()
            plus = plus.rstrip("\n\r")
            qual = qual.rstrip("\n\r")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError("invalid FASTQ record")
            read_id = header[1:].split()[0]
            yield ReadRecord(read_id=read_id, seq=seq, qual=qual)


def _status_name(status: int) -> str:
    return {
        MATCH_INVALID: "invalid",
        MATCH_NONE: "none",
        MATCH_UNIQUE: "unique",
        MATCH_AMBIGUOUS: "ambiguous",
    }.get(status, f"unknown:{status}")


def _edit_kind(observed: str, target: str, dist: int) -> str:
    if dist == 0:
        return "exact"
    if dist != 1:
        return "other"
    if len(observed) == len(target):
        return "substitution"
    if len(observed) == len(target) + 1 and _one_delete_matches(observed, target):
        return "insertion"
    if len(observed) + 1 == len(target) and _one_delete_matches(target, observed):
        return "deletion"
    return "other"


def _one_delete_matches(longer: str, shorter: str) -> bool:
    i = j = edits = 0
    while i < len(longer) and j < len(shorter):
        if longer[i] == shorter[j]:
            i += 1
            j += 1
        else:
            edits += 1
            if edits > 1:
                return False
            i += 1
    return True


def _chunks(it: Iterable[ReadRecord], size: int) -> Iterator[list[ReadRecord]]:
    batch: list[ReadRecord] = []
    for item in it:
        batch.append(item)
        if len(batch) == size:
            yield batch
            batch = []
    if batch:
        yield batch


def _extract(seq: str, start: int, length: int) -> str | None:
    if start < 0 or length < 0:
        return None
    end = start + length
    if end > len(seq):
        return None
    return seq[start:end]


def _target_ambiguity_flags(targets: Sequence[Target], k: int) -> list[int]:
    flags = [0] * len(targets)
    if k < 1:
        seen: dict[str, int] = {}
        for i, target in enumerate(targets):
            prev = seen.get(target.seq)
            if prev is not None:
                flags[prev] = 1
                flags[i] = 1
            else:
                seen[target.seq] = i
        return flags

    for i, j, _dist in _near_target_pairs(targets, min(k, 1)):
        flags[i] = 1
        flags[j] = 1
    if k > 1:
        for i in range(len(targets)):
            if flags[i]:
                continue
            for j in range(i + 1, len(targets)):
                if distance(targets[i].seq, targets[j].seq) <= k:
                    flags[i] = flags[j] = 1
                    break
    return flags


def _neighbors_k1(seq: str) -> Iterator[str]:
    yield seq
    for i, base in enumerate(seq):
        for alt in DNA:
            if alt != base:
                yield seq[:i] + alt + seq[i + 1 :]
    for i in range(len(seq)):
        yield seq[:i] + seq[i + 1 :]
    for i in range(len(seq) + 1):
        for base in DNA:
            yield seq[:i] + base + seq[i:]


def _near_target_pairs(targets: Sequence[Target], k: int) -> Iterator[tuple[int, int, int]]:
    if k > 1:
        seen: set[tuple[int, int]] = set()
        for i in range(len(targets)):
            for j in range(i + 1, len(targets)):
                dist = distance(targets[i].seq, targets[j].seq)
                if dist <= k:
                    seen.add((i, j))
                    yield i, j, dist
        return

    by_seq: dict[str, list[int]] = defaultdict(list)
    for i, target in enumerate(targets):
        by_seq[target.seq].append(i)

    emitted: set[tuple[int, int]] = set()
    for i, target in enumerate(targets):
        for neighbor in _neighbors_k1(target.seq):
            for j in by_seq.get(neighbor, []):
                if i >= j:
                    continue
                pair = (i, j)
                if pair in emitted:
                    continue
                dist = distance(target.seq, targets[j].seq)
                if dist <= k:
                    emitted.add(pair)
                    yield i, j, dist


def _write_assignment_header(fh: TextIO) -> None:
    fh.write(
        "read_id\tobserved_seq\ttarget_id\ttarget_seq\tdistance\tstatus\t"
        "match_count\tsecond_best_distance\tcorrection\n"
    )


def _write_assignment_row(
    fh: TextIO,
    read_id: str,
    observed: str,
    targets: Sequence[Target],
    result: MatchResult,
    correction: str,
) -> None:
    if 0 <= result.target_index < len(targets):
        target = targets[result.target_index]
        target_id = target.target_id
        target_seq = target.seq
    else:
        target_id = ""
        target_seq = ""
    fh.write(
        f"{read_id}\t{observed}\t{target_id}\t{target_seq}\t{result.best_distance}\t"
        f"{_status_name(result.status)}\t{result.match_count}\t{result.second_best_distance}\t{correction}\n"
    )


def command_count(args: argparse.Namespace) -> int:
    targets = _read_targets(args.targets)
    matcher = Matcher([t.seq for t in targets])
    counts = {
        "exact": [0] * len(targets),
        "substitution": [0] * len(targets),
        "insertion": [0] * len(targets),
        "deletion": [0] * len(targets),
        "other": [0] * len(targets),
    }
    summary = {
        "total_reads": 0,
        "assigned_unique": 0,
        "assigned_exact": 0,
        "assigned_corrected": 0,
        "ambiguous": 0,
        "unmatched": 0,
        "invalid": 0,
        "k": args.k,
        "target_start": args.target_start,
        "target_length": args.target_length,
        "n_targets": len(targets),
        "candidates_considered": 0,
        "candidates_verified": 0,
    }

    assignment_fh = _open_text(args.assignments, "wt") if args.assignments else None
    try:
        if assignment_fh is not None:
            _write_assignment_header(assignment_fh)
        for batch in _chunks(_iter_fastq(args.reads), args.batch_size):
            observed: list[str] = []
            valid_positions: list[int] = []
            for pos, record in enumerate(batch):
                seq = _extract(record.seq, args.target_start, args.target_length)
                if seq is None:
                    summary["total_reads"] += 1
                    summary["invalid"] += 1
                    if assignment_fh is not None:
                        invalid = MatchResult(-1, -1, -1, 0, MATCH_INVALID)
                        _write_assignment_row(assignment_fh, record.read_id, "", targets, invalid, "invalid")
                    continue
                observed.append(seq)
                valid_positions.append(pos)

            results, stats = matcher.assign_with_stats(observed, k=args.k)
            summary["candidates_considered"] += stats.candidates_considered
            summary["candidates_verified"] += stats.candidates_verified
            for record_index, obs, result in zip(valid_positions, observed, results):
                record = batch[record_index]
                summary["total_reads"] += 1
                correction = "none"
                if result.status == MATCH_UNIQUE and 0 <= result.target_index < len(targets):
                    target = targets[result.target_index]
                    correction = _edit_kind(obs, target.seq, result.best_distance)
                    counts[correction][result.target_index] += 1
                    summary["assigned_unique"] += 1
                    if result.best_distance == 0:
                        summary["assigned_exact"] += 1
                    else:
                        summary["assigned_corrected"] += 1
                elif result.status == MATCH_AMBIGUOUS:
                    correction = "ambiguous"
                    summary["ambiguous"] += 1
                elif result.status == MATCH_NONE:
                    summary["unmatched"] += 1
                else:
                    correction = "invalid"
                    summary["invalid"] += 1

                if assignment_fh is not None and (
                    result.status != MATCH_AMBIGUOUS or args.ambiguous == "report"
                ):
                    _write_assignment_row(assignment_fh, record.read_id, obs, targets, result, correction)
    finally:
        matcher.close()
        if assignment_fh is not None:
            assignment_fh.close()

    ambiguity_flags = _target_ambiguity_flags(targets, args.k)
    with _open_text(args.out, "wt") as out:
        out.write(
            "target_id\ttarget_seq\tgene\tcount_exact\tcount_corrected_substitution\t"
            "count_corrected_insertion\tcount_corrected_deletion\tcount_corrected_other\t"
            "count_total\tambiguous_nearby\n"
        )
        for i, target in enumerate(targets):
            total = sum(bucket[i] for bucket in counts.values())
            out.write(
                f"{target.target_id}\t{target.seq}\t{target.gene}\t{counts['exact'][i]}\t"
                f"{counts['substitution'][i]}\t{counts['insertion'][i]}\t{counts['deletion'][i]}\t"
                f"{counts['other'][i]}\t{total}\t{ambiguity_flags[i]}\n"
            )

    if args.summary:
        with _open_text(args.summary, "wt") as fh:
            json.dump(summary, fh, indent=2, sort_keys=True)
            fh.write("\n")
    else:
        print(json.dumps(summary, sort_keys=True))
    return 0


def command_audit_targets(args: argparse.Namespace) -> int:
    targets = _read_targets(args.targets)
    pairs = list(_near_target_pairs(targets, args.k))
    duplicates = sum(1 for _i, _j, dist in pairs if dist == 0)
    min_distance = min((dist for _i, _j, dist in pairs), default=None)
    summary = {
        "n_targets": len(targets),
        "k": args.k,
        "duplicates": duplicates,
        "pairs_within_k": len(pairs),
        "unsafe_for_k": bool(pairs),
        "min_observed_pairwise_distance_within_k": min_distance,
    }

    if args.out:
        with _open_text(args.out, "wt") as out:
            out.write("target_id\ttarget_seq\tother_id\tother_seq\tdistance\n")
            for i, j, dist in pairs:
                out.write(
                    f"{targets[i].target_id}\t{targets[i].seq}\t"
                    f"{targets[j].target_id}\t{targets[j].seq}\t{dist}\n"
                )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


def command_validate(args: argparse.Namespace) -> int:
    targets = _read_targets(args.targets)
    matcher = Matcher([t.seq for t in targets])
    checked = 0
    mismatches = 0
    try:
        for batch in _chunks(_iter_fastq(args.reads), args.batch_size):
            observed = []
            for record in batch:
                seq = _extract(record.seq, args.target_start, args.target_length)
                if seq is not None:
                    observed.append(seq)
                if args.sample and len(observed) + checked >= args.sample:
                    break
            if not observed:
                continue
            indexed = matcher.assign(observed, k=args.k)
            oracle = assign(observed, [t.seq for t in targets], k=args.k)
            for obs, fast, slow in zip(observed, indexed, oracle):
                checked += 1
                if fast != slow:
                    mismatches += 1
                    if args.show_mismatches:
                        print(f"mismatch\t{obs}\tindexed={fast}\toracle={slow}", file=sys.stderr)
                if args.sample and checked >= args.sample:
                    break
            if args.sample and checked >= args.sample:
                break
    finally:
        matcher.close()

    summary = {
        "oracle": "native_scan",
        "checked_reads": checked,
        "mismatches": mismatches,
        "k": args.k,
        "target_start": args.target_start,
        "target_length": args.target_length,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if mismatches == 0 else 1


def command_crispr_qc(args: argparse.Namespace) -> int:
    counts = _read_crispr_count_matrix(args.counts)
    sample_qc = _read_sample_qc(args.sample_qc) if args.sample_qc else {}
    library = _read_crispr_library(args.library) if args.library else []
    report = _build_crispr_qc_report(counts, sample_qc, library, args.k, CRISPR_QC_THRESHOLDS)

    if args.out:
        with _open_text(args.out, "wt") as fh:
            json.dump(report, fh, indent=2, sort_keys=True)
            fh.write("\n")
    else:
        print(json.dumps(report, indent=2, sort_keys=True))
    if args.summary_tsv:
        _write_crispr_qc_summary_tsv(report, args.summary_tsv)
    if args.report:
        _write_crispr_qc_html(report, args.report)
    if args.fail_on_review and report["status"] != "pass":
        return 1
    return 0


def command_crispr_namespace(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(prog="dotmatch crispr", description="CRISPR guide-count setup, QC, and AssaySpec helpers.")
    sub = parser.add_subparsers(dest="command", required=True)

    init = sub.add_parser("init", help="write a CRISPR AssaySpec template")
    init.add_argument("--out", required=True)

    infer = sub.add_parser("infer", help="infer a CRISPR guide-count AssaySpec")
    infer.add_argument("--library", required=True)
    infer.add_argument("--reads", required=True)
    infer.add_argument("--sample-id", default="sample")
    infer.add_argument("--out", required=True)
    infer.add_argument("--report", required=True)
    infer.add_argument("--candidates")
    infer.add_argument("--max-reads", type=int, default=50000)
    infer.add_argument("--max-start", type=int, default=32)

    for name in ["check", "plan", "run"]:
        child = sub.add_parser(name, help=f"{name} a CRISPR AssaySpec")
        child.add_argument("spec")

    qc = sub.add_parser("qc", help="evaluate CRISPR guide-count QC and representation metrics")
    _add_crispr_qc_args(qc)

    args = parser.parse_args(list(argv))
    if args.command == "init":
        return command_assay(["init", "--template", "crispr", "--out", args.out])
    if args.command == "infer":
        assay_args = [
            "infer",
            "--mode",
            "count",
            "--assay-type",
            "crispr",
            "--targets",
            args.library,
            "--reads",
            args.reads,
            "--sample-id",
            args.sample_id,
            "--out",
            args.out,
            "--report",
            args.report,
            "--max-reads",
            str(args.max_reads),
            "--max-start",
            str(args.max_start),
        ]
        if args.candidates:
            assay_args.extend(["--candidates", args.candidates])
        return command_assay(assay_args)
    if args.command in {"check", "plan", "run"}:
        return command_assay([args.command, args.spec])
    if args.command == "qc":
        return command_crispr_qc(args)
    parser.error("unreachable")
    return 2


def command_barcode_namespace(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="dotmatch barcode",
        description="Fixed-window barcode inference, demultiplexing, audit, and autopsy.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    infer = sub.add_parser("infer", help="scan fixed windows and rank likely barcode offsets")
    infer.add_argument("--barcodes", required=True, help="TSV/CSV with barcode_id and barcode_seq")
    infer.add_argument("--reads", required=True, help="FASTQ or FASTQ.gz input")
    infer.add_argument("--scan-starts", default="0:30", help="start positions as A:B or comma-separated integers")
    infer.add_argument("--barcode-length", default="auto", help="integer length or auto from barcode sheet")
    infer.add_argument("--sample-reads", type=int, default=50000)
    infer.add_argument("--k", type=int, default=0)
    infer.add_argument("--metric", choices=["hamming", "levenshtein"], default="hamming")
    infer.add_argument("--out", required=True, help="ranked offset scan TSV")
    infer.add_argument("--summary", help="inference summary JSON")

    audit = sub.add_parser("audit", help="run native barcode-library collision audit")
    audit.add_argument("--barcodes", required=True)
    audit.add_argument("--k", type=int, default=1)
    audit.add_argument("--audit-mode", choices=["auto", "exact", "fast"], default="auto")
    audit.add_argument("--out-dir", required=True)

    demux = sub.add_parser("demux", help="run native fixed-window barcode demultiplexing")
    demux.add_argument("--barcodes", required=True)
    demux.add_argument("--reads", required=True)
    demux.add_argument("--barcode-start", type=int, required=True)
    demux.add_argument("--barcode-length", required=True)
    demux.add_argument("--k", type=int, default=0)
    demux.add_argument("--metric", choices=["hamming", "levenshtein"], default="hamming")
    demux.add_argument("--max-correction-qual", type=int)
    demux.add_argument("--out-dir", required=True)
    demux.add_argument("--summary")
    demux.add_argument("--assignments")
    demux.add_argument("--ambiguous-out")
    demux.add_argument("--unmatched-out")

    count = sub.add_parser("count", help="run native fixed-window barcode counting")
    count.add_argument("--barcodes", required=True)
    count.add_argument("--reads", required=True)
    count.add_argument("--barcode-start", type=int, required=True)
    count.add_argument("--barcode-length", type=int, required=True)
    count.add_argument("--k", type=int, default=0)
    count.add_argument("--metric", choices=["hamming", "levenshtein"], default="hamming")
    count.add_argument("--max-correction-qual", type=int)
    count.add_argument("--out", required=True)
    count.add_argument("--summary")
    count.add_argument("--assignments")
    count.add_argument("--format", choices=["dotmatch", "mageck"], default="dotmatch")

    autopsy = sub.add_parser("autopsy", help="infer, audit, demux, inspect unmatched reads, and write evidence report")
    autopsy.add_argument("--barcodes", required=True)
    autopsy.add_argument("--reads", required=True)
    autopsy.add_argument("--scan-starts", default="0:30")
    autopsy.add_argument("--barcode-length", default="auto")
    autopsy.add_argument("--sample-reads", type=int, default=50000)
    autopsy.add_argument("--k-values", default="0,1")
    autopsy.add_argument("--metric", choices=["hamming", "levenshtein"], default="hamming")
    autopsy.add_argument("--max-correction-qual", type=int)
    autopsy.add_argument("--top", type=int, default=100)
    autopsy.add_argument("--out-dir", required=True)

    report = sub.add_parser("report", help="rebuild an autopsy HTML/Markdown report from an output directory")
    report.add_argument("--out-dir", required=True)

    args = parser.parse_args(list(argv))
    try:
        if args.command == "infer":
            candidates = _barcode_infer_candidates(
                args.barcodes,
                args.reads,
                starts=_parse_scan_starts(args.scan_starts),
                length_arg=args.barcode_length,
                sample_reads=args.sample_reads,
                k=args.k,
                metric=args.metric,
            )
            if not candidates:
                raise ValueError("no barcode windows could be scored")
            _write_barcode_offset_scan(candidates, args.out)
            summary = _barcode_inference_summary(candidates, args.metric, args.k)
            if args.summary:
                _write_json(args.summary, summary)
            else:
                print(json.dumps(summary, indent=2, sort_keys=True))
            return 0
        if args.command == "audit":
            return _barcode_audit(args.barcodes, args.k, args.audit_mode, args.out_dir)
        if args.command == "demux":
            return _barcode_demux(args)
        if args.command == "count":
            return _barcode_count(args)
        if args.command == "autopsy":
            return _barcode_autopsy(args)
        if args.command == "report":
            out_dir = Path(args.out_dir)
            provenance = _read_json_if_exists(out_dir / "provenance.json")
            _write_barcode_autopsy_reports(out_dir, provenance)
            return 0
    except BrokenPipeError:
        return 1
    except Exception as exc:
        print(f"dotmatch barcode: {exc}", file=sys.stderr)
        return 2
    parser.error("unreachable")
    return 2


def _parse_scan_starts(value: str) -> list[int]:
    starts: list[int] = []
    for part in value.split(","):
        item = part.strip()
        if not item:
            continue
        if ":" in item:
            left, right = item.split(":", 1)
            begin = int(left)
            end = int(right)
            if begin > end:
                raise ValueError("--scan-starts range must be ascending")
            starts.extend(range(begin, end + 1))
        else:
            starts.append(int(item))
    if not starts:
        raise ValueError("--scan-starts did not contain any positions")
    if any(start < 0 for start in starts):
        raise ValueError("--scan-starts positions must be non-negative")
    return sorted(dict.fromkeys(starts))


def _barcode_lengths(targets: Sequence[Target], length_arg: str) -> list[int]:
    if length_arg == "auto":
        lengths = sorted({len(target.seq) for target in targets})
    else:
        lengths = [int(length_arg)]
    if not lengths or any(length <= 0 for length in lengths):
        raise ValueError("--barcode-length must be auto or a positive integer")
    return lengths


def _barcode_sample_reads(reads: str | Path, max_reads: int) -> list[ReadRecord]:
    if max_reads <= 0:
        raise ValueError("--sample-reads must be positive")
    records: list[ReadRecord] = []
    for record in _iter_fastq(reads):
        records.append(record)
        if len(records) >= max_reads:
            break
    if not records:
        raise ValueError(f"no reads found in {reads}")
    return records


def _barcode_infer_candidates(
    barcodes: str | Path,
    reads: str | Path,
    *,
    starts: Sequence[int],
    length_arg: str,
    sample_reads: int,
    k: int,
    metric: str,
) -> list[BarcodeCandidate]:
    if k < 0:
        raise ValueError("--k must be non-negative")
    if metric == "hamming" and k > 2:
        raise ValueError("hamming barcode inference supports k <= 2")
    targets = _read_targets(barcodes)
    records = _barcode_sample_reads(reads, sample_reads)
    candidates: list[BarcodeCandidate] = []
    for length in _barcode_lengths(targets, length_arg):
        compatible_targets = [target for target in targets if len(target.seq) == length]
        if not compatible_targets:
            continue
        matcher = Matcher([target.seq for target in compatible_targets])
        try:
            for start in starts:
                observed: list[str] = []
                invalid = 0
                for record in records:
                    seq = _extract(record.seq, start, length)
                    if seq is None:
                        invalid += 1
                    else:
                        observed.append(seq)
                results = matcher.assign(observed, k=k)
                unique = exact = ambiguous = no_match = 0
                for obs, result in zip(observed, results):
                    if result.status == MATCH_UNIQUE:
                        unique += 1
                        if result.best_distance == 0:
                            exact += 1
                    elif result.status == MATCH_AMBIGUOUS:
                        ambiguous += 1
                    elif result.status == MATCH_NONE:
                        no_match += 1
                    else:
                        invalid += 1
                sampled = len(records)
                valid = sampled - invalid
                assignment_rate = unique / sampled if sampled else 0.0
                exact_rate = exact / sampled if sampled else 0.0
                ambiguous_rate = ambiguous / sampled if sampled else 0.0
                no_match_rate = no_match / sampled if sampled else 0.0
                invalid_rate = invalid / sampled if sampled else 0.0
                score = assignment_rate - ambiguous_rate - invalid_rate
                candidates.append(
                    BarcodeCandidate(
                        start=start,
                        length=length,
                        score=score,
                        sampled_reads=sampled,
                        valid_reads=valid,
                        unique=unique,
                        exact=exact,
                        ambiguous=ambiguous,
                        no_match=no_match,
                        invalid=invalid,
                        assignment_rate=assignment_rate,
                        exact_rate=exact_rate,
                        ambiguous_rate=ambiguous_rate,
                        no_match_rate=no_match_rate,
                        invalid_rate=invalid_rate,
                    )
                )
        finally:
            matcher.close()
    return sorted(candidates, key=lambda c: (-c.score, -c.assignment_rate, c.ambiguous_rate, c.start, c.length))


def _write_barcode_offset_scan(candidates: Sequence[BarcodeCandidate], path: str | Path) -> None:
    columns = [
        "start",
        "length",
        "score",
        "sampled_reads",
        "valid_reads",
        "unique",
        "exact",
        "ambiguous",
        "no_match",
        "invalid",
        "assignment_rate",
        "exact_rate",
        "ambiguous_rate",
        "no_match_rate",
        "invalid_rate",
    ]
    with _open_text(path, "wt") as fh:
        fh.write("\t".join(columns) + "\n")
        for candidate in candidates:
            values = {
                "start": candidate.start,
                "length": candidate.length,
                "score": _fmt_rate(candidate.score),
                "sampled_reads": candidate.sampled_reads,
                "valid_reads": candidate.valid_reads,
                "unique": candidate.unique,
                "exact": candidate.exact,
                "ambiguous": candidate.ambiguous,
                "no_match": candidate.no_match,
                "invalid": candidate.invalid,
                "assignment_rate": _fmt_rate(candidate.assignment_rate),
                "exact_rate": _fmt_rate(candidate.exact_rate),
                "ambiguous_rate": _fmt_rate(candidate.ambiguous_rate),
                "no_match_rate": _fmt_rate(candidate.no_match_rate),
                "invalid_rate": _fmt_rate(candidate.invalid_rate),
            }
            fh.write("\t".join(str(values[column]) for column in columns) + "\n")


def _barcode_inference_summary(candidates: Sequence[BarcodeCandidate], metric: str, k: int) -> dict[str, object]:
    best = candidates[0]
    second = candidates[1] if len(candidates) > 1 else None
    separation = best.score - second.score if second is not None else best.score
    warnings: list[str] = []
    if best.assignment_rate == 0:
        warnings.append("no reads matched the barcode sheet in any scanned window")
    elif best.assignment_rate < 0.20:
        warnings.append("best scanned window has a low exact assignment rate; review the barcode sheet, read orientation, and window before treating the offset as ready")
    return {
        "schema_version": 1,
        "command": "barcode infer",
        "recommended_start": best.start,
        "recommended_length": best.length,
        "metric": metric,
        "k": k,
        "assignment_rate": best.assignment_rate,
        "ambiguous_rate": best.ambiguous_rate,
        "no_match_rate": best.no_match_rate,
        "invalid_rate": best.invalid_rate,
        "score": best.score,
        "separation_from_second_best": separation,
        "candidate_count": len(candidates),
        "warnings": warnings,
        "status": "review" if warnings else "ready",
    }


def _fmt_rate(value: float) -> str:
    return f"{value:.8f}"


def _write_json(path: str | Path, data: object) -> None:
    with _open_text(path, "wt") as fh:
        json.dump(data, fh, indent=2, sort_keys=True)
        fh.write("\n")


def _read_json_if_exists(path: str | Path) -> dict[str, object]:
    p = Path(path)
    if not p.exists():
        return {}
    with _open_text(p) as fh:
        data = json.load(fh)
    return data if isinstance(data, dict) else {}


def _native_completed(argv: Sequence[str]) -> dict[str, object]:
    native = find_native_cli()
    command = [str(native), *argv]
    completed = subprocess.run(command, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return {
        "argv": command,
        "exit_code": int(completed.returncode),
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def _require_native_success(result: dict[str, object]) -> None:
    if result["exit_code"] != 0:
        command = " ".join(str(part) for part in result["argv"])
        stderr = str(result.get("stderr", "")).strip()
        raise RuntimeError(f"native command failed ({command}): {stderr}")


def _barcode_audit(barcodes: str, k: int, audit_mode: str, out_dir: str) -> int:
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    result = _native_completed(["audit", "--targets", barcodes, "--k", str(k), "--audit-mode", audit_mode, "--out-dir", out_dir])
    _require_native_success(result)
    _copy_audit_aliases(Path(out_dir), Path(out_dir))
    return 0


def _copy_audit_aliases(source_dir: Path, dest_dir: Path) -> None:
    aliases = [
        ("collision_pairs.tsv", "collision_graph.tsv"),
        ("target_safety.tsv", "correction_safety.tsv"),
    ]
    dest_dir.mkdir(parents=True, exist_ok=True)
    for source_name, dest_name in aliases:
        source = source_dir / source_name
        if source.exists():
            shutil.copyfile(source, dest_dir / dest_name)
        elif not (dest_dir / dest_name).exists():
            (dest_dir / dest_name).write_text("", encoding="utf-8")


def _barcode_demux(args: argparse.Namespace) -> int:
    argv = [
        "demux",
        "--barcodes",
        args.barcodes,
        "--reads",
        args.reads,
        "--barcode-start",
        str(args.barcode_start),
        "--barcode-length",
        str(args.barcode_length),
        "--k",
        str(args.k),
        "--metric",
        args.metric,
        "--out-dir",
        args.out_dir,
    ]
    if args.max_correction_qual is not None:
        argv.extend(["--max-correction-qual", str(args.max_correction_qual)])
    if args.summary:
        argv.extend(["--summary", args.summary])
    if args.assignments:
        argv.extend(["--assignments", args.assignments])
    if args.ambiguous_out:
        argv.extend(["--ambiguous-out", args.ambiguous_out])
    if args.unmatched_out:
        argv.extend(["--unmatched-out", args.unmatched_out])
    result = _native_completed(argv)
    _require_native_success(result)
    return 0


def _barcode_count(args: argparse.Namespace) -> int:
    argv = [
        "count",
        "--targets",
        args.barcodes,
        "--reads",
        args.reads,
        "--target-start",
        str(args.barcode_start),
        "--target-length",
        str(args.barcode_length),
        "--k",
        str(args.k),
        "--metric",
        args.metric,
        "--out",
        args.out,
        "--format",
        args.format,
    ]
    if args.max_correction_qual is not None:
        argv.extend(["--max-correction-qual", str(args.max_correction_qual)])
    if args.summary:
        argv.extend(["--summary", args.summary])
    if args.assignments:
        argv.extend(["--assignments", args.assignments])
    result = _native_completed(argv)
    _require_native_success(result)
    return 0


def _parse_k_values(value: str) -> list[int]:
    values = [int(part.strip()) for part in value.split(",") if part.strip()]
    if not values:
        raise ValueError("--k-values must contain at least one integer")
    if any(k < 0 for k in values):
        raise ValueError("--k-values must be non-negative")
    return sorted(dict.fromkeys(values))


def _barcode_autopsy(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    commands: list[dict[str, object]] = []
    artifacts: list[str] = []

    candidates = _barcode_infer_candidates(
        args.barcodes,
        args.reads,
        starts=_parse_scan_starts(args.scan_starts),
        length_arg=args.barcode_length,
        sample_reads=args.sample_reads,
        k=0,
        metric=args.metric,
    )
    if not candidates:
        raise ValueError("no barcode windows could be scored")
    offset_scan = out_dir / "offset_scan.tsv"
    _write_barcode_offset_scan(candidates, offset_scan)
    artifacts.append(offset_scan.name)
    inference = _barcode_inference_summary(candidates, args.metric, 0)

    k_values = _parse_k_values(args.k_values)
    audit_k = max(k_values)
    audit_dir = out_dir / "audit"
    audit_result = _native_completed(
        ["audit", "--targets", args.barcodes, "--k", str(max(1, audit_k)), "--audit-mode", "auto", "--out-dir", str(audit_dir)]
    )
    commands.append(audit_result)
    _require_native_success(audit_result)
    _copy_audit_aliases(audit_dir, out_dir)
    artifacts.extend(["audit/", "collision_graph.tsv", "correction_safety.tsv"])

    best = candidates[0]
    demux_dir = out_dir / "demuxed"
    summary_path = out_dir / "summary.json"
    assignments_path = out_dir / "assignments.tsv"
    ambiguous_path = out_dir / "ambiguous.fastq"
    unmatched_path = out_dir / "unmatched.fastq"
    demux_argv = [
        "demux",
        "--barcodes",
        args.barcodes,
        "--reads",
        args.reads,
        "--barcode-start",
        str(best.start),
        "--barcode-length",
        str(args.barcode_length),
        "--k",
        str(audit_k),
        "--metric",
        args.metric,
        "--out-dir",
        str(demux_dir),
        "--summary",
        str(summary_path),
        "--assignments",
        str(assignments_path),
        "--ambiguous-out",
        str(ambiguous_path),
        "--unmatched-out",
        str(unmatched_path),
    ]
    if args.max_correction_qual is not None:
        demux_argv.extend(["--max-correction-qual", str(args.max_correction_qual)])
    demux_result = _native_completed(demux_argv)
    commands.append(demux_result)
    _require_native_success(demux_result)
    artifacts.extend(["demuxed/", "summary.json", "assignments.tsv", "ambiguous.fastq", "unmatched.fastq"])

    top_unmatched = out_dir / "top_unmatched.tsv"
    inspect_result = _native_completed(
        [
            "inspect-unmatched",
            "--targets",
            args.barcodes,
            "--reads",
            args.reads,
            "--target-start",
            str(best.start),
            "--target-length",
            str(best.length),
            "--k",
            str(min(audit_k, 1)),
            "--top",
            str(args.top),
            "--out",
            str(top_unmatched),
        ]
    )
    commands.append(inspect_result)
    _require_native_success(inspect_result)
    artifacts.append("top_unmatched.tsv")

    summary = _read_json_if_exists(summary_path)
    _write_barcode_sample_qc(out_dir / "sample_qc.tsv", summary)
    _write_barcode_counts(out_dir / "barcode_counts.tsv", demux_dir)
    audit_summary = _read_json_if_exists(audit_dir / "audit_summary.json")
    findings = _barcode_findings(inference, summary, audit_summary)
    _write_barcode_findings(out_dir / "findings.tsv", findings)
    artifacts.extend(["sample_qc.tsv", "barcode_counts.tsv", "findings.tsv"])
    artifacts.extend(["provenance.json", "multiqc_dotmatch_barcode_mqc.yaml"])

    provenance = {
        "schema_version": 1,
        "workflow": "barcode_autopsy",
        "barcodes": str(args.barcodes),
        "reads": str(args.reads),
        "recommended_start": best.start,
        "recommended_length": best.length,
        "k_values": k_values,
        "inference": inference,
        "commands": [_public_command_result(item) for item in commands],
        "artifacts": artifacts,
    }
    _write_json(out_dir / "provenance.json", provenance)
    _write_barcode_multiqc(out_dir / "multiqc_dotmatch_barcode_mqc.yaml", summary, inference)
    _write_barcode_autopsy_reports(out_dir, provenance)
    print(out_dir / "report.html")
    return 0


def _public_command_result(result: dict[str, object]) -> dict[str, object]:
    argv = list(result.get("argv", [])) if isinstance(result.get("argv", []), list) else []
    if argv:
        argv[0] = Path(str(argv[0])).name
    return {
        "argv": argv,
        "exit_code": result.get("exit_code", 0),
        "stderr": str(result.get("stderr", "")).strip(),
    }


def _write_barcode_sample_qc(path: Path, summary: dict[str, object]) -> None:
    total = _as_int(summary.get("total_reads"))
    assigned = _as_int(summary.get("assigned_unique"))
    ambiguous = _as_int(summary.get("ambiguous"))
    unmatched = _as_int(summary.get("unmatched"))
    invalid = _as_int(summary.get("invalid"))
    corrected = _as_int(summary.get("assigned_corrected"))
    exact = _as_int(summary.get("assigned_exact"))
    columns = [
        "sample_id",
        "total_reads",
        "assigned_unique",
        "assigned_exact",
        "assigned_corrected",
        "ambiguous",
        "unmatched",
        "invalid",
        "assignment_rate",
        "ambiguous_rate",
        "no_match_rate",
        "invalid_rate",
        "top_failure_reason",
    ]
    values = {
        "sample_id": "pooled",
        "total_reads": total,
        "assigned_unique": assigned,
        "assigned_exact": exact,
        "assigned_corrected": corrected,
        "ambiguous": ambiguous,
        "unmatched": unmatched,
        "invalid": invalid,
        "assignment_rate": _fmt_rate(assigned / total if total else 0.0),
        "ambiguous_rate": _fmt_rate(ambiguous / total if total else 0.0),
        "no_match_rate": _fmt_rate(unmatched / total if total else 0.0),
        "invalid_rate": _fmt_rate(invalid / total if total else 0.0),
        "top_failure_reason": _barcode_failure_reason(assigned, ambiguous, unmatched, invalid, total),
    }
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(columns) + "\n")
        fh.write("\t".join(str(values[column]) for column in columns) + "\n")


def _as_int(value: object) -> int:
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    if isinstance(value, str) and value.strip():
        return int(float(value))
    return 0


def _barcode_failure_reason(assigned: int, ambiguous: int, unmatched: int, invalid: int, total: int) -> str:
    if total == 0:
        return "no_reads"
    if invalid / total > 0.02:
        return "invalid_window"
    if ambiguous / total > 0.05:
        return "ambiguous_collision"
    if unmatched / total > 0.15:
        return "high_no_match"
    if assigned / total >= 0.80:
        return "clean"
    return "review"


def _barcode_findings(
    inference: dict[str, object],
    summary: dict[str, object],
    audit_summary: dict[str, object],
) -> list[dict[str, str]]:
    total = _as_int(summary.get("total_reads"))
    assigned = _as_int(summary.get("assigned_unique"))
    ambiguous = _as_int(summary.get("ambiguous"))
    unmatched = _as_int(summary.get("unmatched"))
    invalid = _as_int(summary.get("invalid"))
    findings: list[dict[str, str]] = []

    if str(inference.get("status") or "") == "review":
        findings.append(
            {
                "finding": "low_confidence_offset",
                "severity": "review",
                "evidence": f"offset_scan.tsv best assignment_rate={_fmt_float(inference.get('assignment_rate'))}",
                "meaning": "The sampled scan has a weak best window; the assay specification may be wrong or incomplete.",
                "next_action": "Review read side, barcode start, barcode length, target sheet, and orientation before production use.",
            }
        )
    if total and assigned / total < 0.80:
        findings.append(
            {
                "finding": "low_assignment_rate",
                "severity": "review",
                "evidence": f"summary.json assigned_unique={assigned} total_reads={total}",
                "meaning": "Fewer than 80% of reads were uniquely assigned under the selected fixed-window semantics.",
                "next_action": "Inspect offset_scan.tsv, top_unmatched.tsv, and the barcode sheet before trusting rescued reads.",
            }
        )
    if total and ambiguous / total > 0.05:
        findings.append(
            {
                "finding": "high_ambiguity_rate",
                "severity": "review",
                "evidence": f"summary.json ambiguous={ambiguous} total_reads={total}",
                "meaning": "A material fraction of reads is compatible with more than one barcode.",
                "next_action": "Do not rescue ambiguous reads into either sample without changing the barcode design or assignment policy.",
            }
        )
    if total and unmatched / total > 0.15:
        findings.append(
            {
                "finding": "high_no_match_rate",
                "severity": "review",
                "evidence": f"summary.json unmatched={unmatched} total_reads={total}",
                "meaning": "Many reads have no compatible barcode in the configured window.",
                "next_action": "Check the barcode sheet, read orientation, offset, and top unmatched sequences.",
            }
        )
    if total and invalid / total > 0.02:
        findings.append(
            {
                "finding": "high_invalid_window_rate",
                "severity": "review",
                "evidence": f"summary.json invalid={invalid} total_reads={total}",
                "meaning": "The configured window cannot be extracted from a notable fraction of reads.",
                "next_action": "Check read length, trimming, read side, and fixed-window coordinates.",
            }
        )
    if _audit_is_unsafe(audit_summary):
        findings.append(
            {
                "finding": "unsafe_one_edit_correction",
                "severity": "review",
                "evidence": "audit/audit_summary.json reports collision risk for k=1",
                "meaning": "At least one barcode is not safe for one-edit correction.",
                "next_action": "Use k=0 or redesign/fix colliding barcodes before enabling one-edit rescue.",
            }
        )
    if not findings:
        findings.append(
            {
                "finding": "no_review_findings",
                "severity": "pass",
                "evidence": "summary.json, offset_scan.tsv, and audit/audit_summary.json",
                "meaning": "No configured barcode-autopsy thresholds were crossed.",
                "next_action": "Review report.html and provenance.json, then use outputs in the downstream workflow.",
            }
        )
    return findings


def _fmt_float(value: object) -> str:
    if isinstance(value, (int, float)):
        return _fmt_rate(float(value))
    return str(value if value is not None else "")


def _audit_is_unsafe(audit_summary: dict[str, object]) -> bool:
    for key in ["safe_at_k1", "safe_for_k", "safe_for_audited_radius"]:
        if audit_summary.get(key) is False:
            return True
        if str(audit_summary.get(key) or "").strip().lower() in {"false", "no", "0"}:
            return True
    for key in ["risk_pairs_for_k1", "pairs_within_k", "one_edit_collision_pairs"]:
        if _as_int(audit_summary.get(key)) > 0:
            return True
    return False


def _write_barcode_findings(path: Path, findings: Sequence[dict[str, str]]) -> None:
    columns = ["finding", "severity", "evidence", "meaning", "next_action"]
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(columns) + "\n")
        for finding in findings:
            fh.write("\t".join(str(finding.get(column, "")) for column in columns) + "\n")


def _write_barcode_counts(path: Path, demux_dir: Path) -> None:
    with path.open("w", encoding="utf-8") as fh:
        fh.write("barcode_id\treads\n")
        if not demux_dir.exists():
            return
        for fastq in sorted(demux_dir.glob("*.fastq*")):
            fh.write(f"{fastq.name.split('.')[0]}\t{_count_fastq_records(fastq)}\n")


def _count_fastq_records(path: Path) -> int:
    lines = 0
    with _open_text(path) as fh:
        for _line in fh:
            lines += 1
    return lines // 4


def _write_barcode_multiqc(path: Path, summary: dict[str, object], inference: dict[str, object]) -> None:
    text = (
        "id: dotmatch_barcode_qc\n"
        "section_name: DotMatch Barcode QC\n"
        "plot_type: table\n"
        "data:\n"
        "  pooled:\n"
        f"    assignment_rate: {_fmt_rate(float(inference.get('assignment_rate', 0.0)))}\n"
        f"    ambiguous: {_as_int(summary.get('ambiguous'))}\n"
        f"    unmatched: {_as_int(summary.get('unmatched'))}\n"
        f"    invalid: {_as_int(summary.get('invalid'))}\n"
    )
    path.write_text(text, encoding="utf-8")


def _write_barcode_autopsy_reports(out_dir: Path, provenance: dict[str, object]) -> None:
    inference = provenance.get("inference") if isinstance(provenance.get("inference"), dict) else {}
    commands = provenance.get("commands") if isinstance(provenance.get("commands"), list) else []
    artifacts = provenance.get("artifacts") if isinstance(provenance.get("artifacts"), list) else []
    run_summary = _read_json_if_exists(out_dir / "summary.json")
    start = inference.get("recommended_start", provenance.get("recommended_start", ""))
    length = inference.get("recommended_length", provenance.get("recommended_length", ""))
    assignment_rate = inference.get("assignment_rate", "")
    inference_status = inference.get("status", "ready")
    inference_warnings = inference.get("warnings") if isinstance(inference.get("warnings"), list) else []
    ambiguous_rate = inference.get("ambiguous_rate", "")
    no_match_rate = inference.get("no_match_rate", "")
    total_reads = _as_int(run_summary.get("total_reads"))
    assigned_unique = _as_int(run_summary.get("assigned_unique"))
    assigned_exact = _as_int(run_summary.get("assigned_exact"))
    assigned_corrected = _as_int(run_summary.get("assigned_corrected"))
    ambiguous = _as_int(run_summary.get("ambiguous"))
    unmatched = _as_int(run_summary.get("unmatched"))
    invalid = _as_int(run_summary.get("invalid"))
    run_assignment_rate = _fmt_rate(assigned_unique / total_reads if total_reads else 0.0)
    run_ambiguous_rate = _fmt_rate(ambiguous / total_reads if total_reads else 0.0)
    run_no_match_rate = _fmt_rate(unmatched / total_reads if total_reads else 0.0)
    diagnosis = _barcode_failure_reason(assigned_unique, ambiguous, unmatched, invalid, total_reads)
    interpretation = _barcode_report_interpretation(diagnosis, str(inference_status))
    command_lines = [
        " ".join(str(part) for part in command.get("argv", [])) if isinstance(command, dict) else ""
        for command in commands
    ]
    artifact_lines = [str(item) for item in artifacts]
    findings = _read_barcode_findings(out_dir / "findings.tsv")
    finding_lines = [
        f"{row['severity']}: {row['finding']} - {row['meaning']} Next action: {row['next_action']}"
        for row in findings
    ]
    md = "\n".join(
        [
            "# Barcode Troubleshooting Report",
            "",
            "DotMatch reports how fixed-window barcode reads were assigned, rejected, or flagged for review.",
            "",
            "Speed is reported only after the comparator settings are documented.",
            "",
            "## Decision Summary",
            "",
            f"- Inference status: `{inference_status}`",
            f"- Assignment rate: `{run_assignment_rate}`",
            f"- Top finding: `{findings[0]['finding'] if findings else diagnosis}`",
            f"- Primary report: `report.html`",
            "",
            "## Comparator check",
            "",
            "Use this section with recorded comparator settings for the same barcode window, length, and correction policy.",
            "",
            "## Offset check",
            "",
            f"Highest-scoring sampled barcode window: start={start}, length={length}.",
            f"Exact assignment rate at that window: {assignment_rate}. Inference status: {inference_status}.",
            *(f"Warning: {warning}." for warning in inference_warnings),
            "",
            "## Barcode safety audit",
            "",
            "The audit output reports duplicate and nearby barcode pairs before one-edit correction is trusted.",
            "",
            "## Assignment summary",
            "",
            f"Run assignment rate: {run_assignment_rate}.",
            f"Exact assignments: {assigned_exact}. Corrected assignments: {assigned_corrected}.",
            f"Ambiguous reads: {ambiguous} ({run_ambiguous_rate}). Unmatched reads: {unmatched} ({run_no_match_rate}). Invalid windows: {invalid}.",
            f"Top failure reason: {diagnosis}. Top unmatched sequences are written to `top_unmatched.tsv`.",
            "",
            "### What this means",
            "",
            interpretation["meaning"],
            "",
            "### Next action",
            "",
            interpretation["next_action"],
            "",
            "Do not rescue ambiguous reads into either sample without changing the barcode design or assignment policy.",
            "",
            "## Findings",
            "",
            *(f"- {line}" for line in finding_lines),
            "",
            "## Workflow handoff",
            "",
            "Artifacts are written as stable TSV, JSON, FASTQ, HTML, and MultiQC custom-content inputs.",
            "",
            "## QC Checklist",
            "",
            "- Exact command provenance is recorded in `provenance.json`.",
            "- Offset evidence is recorded in `offset_scan.tsv`.",
            "- Barcode collision safety is recorded under `audit/` and summarized in `correction_safety.tsv`.",
            "- Ambiguous and unmatched reads are retained when requested instead of being silently assigned.",
            "- Benchmark notes stay tied to documented comparator settings.",
            "",
            "## Commands",
            "",
            *(f"- `{line}`" for line in command_lines if line),
            "",
            "## Artifacts",
            "",
            *(f"- `{line}`" for line in artifact_lines),
            "",
        ]
    )
    (out_dir / "report.md").write_text(md, encoding="utf-8")
    html_commands = "\n".join(f"<li><code>{html.escape(line)}</code></li>" for line in command_lines if line) or "<li>No commands recorded.</li>"
    html_artifacts = "\n".join(f"<li><code>{html.escape(line)}</code></li>" for line in artifact_lines) or "<li>No artifacts recorded.</li>"
    html_findings = "\n".join(f"<li>{html.escape(line)}</li>" for line in finding_lines) or "<li>No review findings.</li>"
    html_doc = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>DotMatch Barcode Troubleshooting Report</title>
  <style>
    body {{ margin: 2rem; color: #111816; font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }}
    .panel {{ border-top: 1px solid #dfe7e2; padding: 1rem 0; }}
    code {{ background: #eef6f2; padding: 0.1rem 0.25rem; }}
  </style>
</head>
<body>
  <h1>Barcode Troubleshooting Report</h1>
  <p>DotMatch reports how fixed-window barcode reads were assigned, rejected, or flagged for review.</p>
  <p><strong>Speed is reported only after the comparator settings are documented.</strong></p>
  <section class="panel"><h2>Decision Summary</h2><p>Inference status: <strong>{html.escape(str(inference_status))}</strong>; assignment rate: <strong>{html.escape(run_assignment_rate)}</strong>; top finding: <strong>{html.escape(findings[0]['finding'] if findings else diagnosis)}</strong>.</p></section>
  <section class="panel"><h2>Comparator check</h2><p>Use recorded comparator settings for the same barcode window, length, and correction policy before interpreting speed.</p></section>
  <section class="panel"><h2>Offset check</h2><p>Highest-scoring sampled window: start={html.escape(str(start))}, length={html.escape(str(length))}; exact assignment rate={html.escape(str(assignment_rate))}; status={html.escape(str(inference_status))}.</p></section>
  <section class="panel"><h2>Barcode safety audit</h2><p>Audit outputs identify duplicate and nearby barcode pairs before one-edit correction is enabled.</p></section>
  <section class="panel"><h2>Assignment summary</h2><p>Run assignment rate={html.escape(run_assignment_rate)}; ambiguous={html.escape(str(ambiguous))} ({html.escape(run_ambiguous_rate)}); unmatched={html.escape(str(unmatched))} ({html.escape(run_no_match_rate)}); invalid={html.escape(str(invalid))}. Top failure reason: {html.escape(diagnosis)}.</p><h3>What this means</h3><p>{html.escape(interpretation["meaning"])}</p><h3>Next action</h3><p>{html.escape(interpretation["next_action"])}</p><p>Do not rescue ambiguous reads into either sample without changing the barcode design or assignment policy.</p></section>
  <section class="panel"><h2>Findings</h2><ul>{html_findings}</ul></section>
  <section class="panel"><h2>Workflow handoff</h2><p>Stable TSV, JSON, FASTQ, HTML, and MultiQC custom-content artifacts are emitted for workflow systems.</p></section>
  <h2>QC Checklist</h2><ul><li>Exact command provenance is recorded in <code>provenance.json</code>.</li><li>Offset evidence is recorded in <code>offset_scan.tsv</code>.</li><li>Barcode collision safety is recorded under <code>audit/</code> and summarized in <code>correction_safety.tsv</code>.</li><li>Ambiguous and unmatched reads are retained when requested instead of being silently assigned.</li><li>Benchmark notes stay tied to documented comparator settings.</li></ul>
  <h2>Commands</h2><ul>{html_commands}</ul>
  <h2>Artifacts</h2><ul>{html_artifacts}</ul>
</body>
</html>
"""
    (out_dir / "report.html").write_text(html_doc, encoding="utf-8")


def _read_barcode_findings(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def _barcode_report_interpretation(diagnosis: str, inference_status: str) -> dict[str, str]:
    if inference_status == "review":
        return {
            "meaning": "The sampled fixed-window scan found a best window, but the evidence is weak enough that the assay specification should be reviewed before production use.",
            "next_action": "Check the barcode sheet, read side, offset, barcode length, and expected orientation; rerun inference after correcting the assay specification.",
        }
    if diagnosis == "ambiguous_collision":
        return {
            "meaning": "A material fraction of reads is compatible with more than one barcode target under the configured edit radius.",
            "next_action": "Audit the barcode library and lower k or redesign colliding barcodes before enabling correction.",
        }
    if diagnosis == "high_no_match":
        return {
            "meaning": "Many reads have no compatible barcode in the configured fixed window.",
            "next_action": "Inspect top unmatched sequences, scan offsets, and confirm the barcode sheet matches this run.",
        }
    if diagnosis == "invalid_window":
        return {
            "meaning": "The configured barcode window cannot be extracted from a notable fraction of reads.",
            "next_action": "Check read length, trimming, read side, and fixed-window coordinates.",
        }
    if diagnosis == "clean":
        return {
            "meaning": "The run has high unique assignment under the configured fixed-window semantics.",
            "next_action": "Review audit and unmatched artifacts, then use the manifest and report in the downstream workflow.",
        }
    return {
        "meaning": "The barcode run needs review before treating rescued or unmatched reads as final biological evidence.",
        "next_action": "Use the offset scan, collision audit, assignment table, and top unmatched reads to identify the assay or sample-sheet issue.",
    }


def _add_crispr_qc_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--counts", required=True, help="MAGeCK-style count matrix")
    parser.add_argument("--sample-qc", help="DotMatch sample_qc.tsv")
    parser.add_argument("--library", help="CRISPR guide library CSV/TSV")
    parser.add_argument("--k", type=int, default=1, help="edit radius used for library collision audit")
    parser.add_argument("--out", help="JSON QC report output; stdout is used when omitted")
    parser.add_argument("--summary-tsv", help="sample-level summary TSV output")
    parser.add_argument("--report", help="self-contained HTML report output")
    parser.add_argument("--fail-on-review", action="store_true", help="exit nonzero when any QC warning is emitted")


def _read_crispr_count_matrix(path: str | Path) -> dict[str, object]:
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or len(reader.fieldnames) < 3:
            raise ValueError("count matrix must have guide, gene, and at least one sample column")
        guide_col = reader.fieldnames[0]
        gene_col = reader.fieldnames[1]
        sample_cols = reader.fieldnames[2:]
        guides: list[dict[str, object]] = []
        sample_counts = {sample: [] for sample in sample_cols}
        for row in reader:
            guide_id = str(row.get(guide_col, "")).strip()
            gene = str(row.get(gene_col, "")).strip()
            if not guide_id:
                raise ValueError("count matrix contains an empty guide id")
            guide_counts: dict[str, int] = {}
            for sample in sample_cols:
                value = int(float(str(row.get(sample, "0") or "0")))
                if value < 0:
                    raise ValueError(f"negative count for {guide_id}/{sample}")
                guide_counts[sample] = value
                sample_counts[sample].append(value)
            guides.append({"id": guide_id, "gene": gene, "counts": guide_counts})
    if not guides:
        raise ValueError("count matrix contains no guides")
    return {
        "guide_col": guide_col,
        "gene_col": gene_col,
        "samples": sample_cols,
        "guides": guides,
        "sample_counts": sample_counts,
    }


def _read_sample_qc(path: str | Path) -> dict[str, dict[str, float | str]]:
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "sample_id" not in reader.fieldnames:
            raise ValueError("sample_qc.tsv must contain a sample_id column")
        rows: dict[str, dict[str, float | str]] = {}
        for row in reader:
            sample_id = str(row.get("sample_id", "")).strip()
            if not sample_id:
                raise ValueError("sample_qc.tsv contains an empty sample_id")
            parsed: dict[str, float | str] = {}
            for key, value in row.items():
                if key == "sample_id":
                    continue
                text = str(value or "")
                try:
                    parsed[key] = float(text)
                except ValueError:
                    parsed[key] = text
            rows[sample_id] = parsed
    return rows


def _read_crispr_library(path: str | Path) -> list[Target]:
    delimiter = "," if str(path).lower().endswith(".csv") else "\t"
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError("CRISPR library must contain a header")
        normalized = {name.lower().replace("_", "."): name for name in reader.fieldnames}
        id_col = _first_existing(normalized, ["id", "sgrna", "sgrna.id", "guide", "guide.id", "target.id"]) or reader.fieldnames[0]
        seq_col = _first_existing(
            normalized,
            ["grna.sequence", "sgrna.sequence", "guide.sequence", "target.seq", "sequence", "seq"],
        )
        if seq_col is None:
            if delimiter == "\t":
                return _read_targets(path)
            raise ValueError("CRISPR library must contain a guide sequence column")
        gene_col = _first_existing(normalized, ["gene", "gene.symbol", "target.gene"]) or ""
        targets: list[Target] = []
        for i, row in enumerate(reader):
            guide_id = str(row.get(id_col, "") or f"guide_{i}").strip()
            seq = str(row.get(seq_col, "") or "").strip().upper()
            gene = str(row.get(gene_col, "") or "").strip() if gene_col else ""
            if not guide_id or not seq:
                raise ValueError("CRISPR library contains empty guide id or sequence")
            targets.append(Target(guide_id, seq, gene))
    if not targets:
        raise ValueError("CRISPR library contains no guides")
    return targets


def _first_existing(normalized: dict[str, str], candidates: Sequence[str]) -> str | None:
    for candidate in candidates:
        if candidate in normalized:
            return normalized[candidate]
    return None


def _build_crispr_qc_report(
    counts: dict[str, object],
    sample_qc: dict[str, dict[str, float | str]],
    library: Sequence[Target],
    k: int,
    thresholds: dict[str, float],
) -> dict[str, object]:
    samples = list(counts["samples"])  # type: ignore[arg-type]
    sample_counts: dict[str, list[int]] = counts["sample_counts"]  # type: ignore[assignment]
    guide_count = len(counts["guides"])  # type: ignore[arg-type]
    warnings: list[dict[str, object]] = []
    sample_report: dict[str, dict[str, object]] = {}

    if not sample_qc:
        warnings.append(
            {
                "code": "sample_qc_not_provided",
                "severity": "review",
                "scope": "inputs",
                "message": "sample_qc.tsv was not provided; assignment, ambiguity, no-match, and invalid-rate QC were not evaluated.",
            }
        )
    if not library:
        warnings.append(
            {
                "code": "library_not_provided",
                "severity": "review",
                "scope": "inputs",
                "message": "Guide library was not provided; guide collision, duplicate, and sequence-content QC were not evaluated.",
            }
        )

    for sample in samples:
        values = sample_counts[sample]
        total = sum(values)
        zero_count = sum(1 for value in values if value == 0)
        nonzero = guide_count - zero_count
        top_n = max(1, math.ceil(guide_count * 0.01))
        top_fraction = sum(sorted(values, reverse=True)[:top_n]) / total if total else 0.0
        metrics: dict[str, object] = dict(sample_qc.get(sample, {}))
        metrics.update({
            "total_count": total,
            "nonzero_guides": nonzero,
            "zero_count_guides": zero_count,
            "coverage_fraction": nonzero / guide_count if guide_count else 0.0,
            "zero_count_fraction": zero_count / guide_count if guide_count else 0.0,
            "gini_index": _gini(values),
            "top_1pct_fraction": top_fraction,
        })
        if "invalid_rate" not in metrics and isinstance(metrics.get("invalid_reads"), (int, float)):
            total_reads = metrics.get("total_reads")
            metrics["invalid_rate"] = float(metrics["invalid_reads"]) / float(total_reads) if isinstance(total_reads, (int, float)) and total_reads else 0.0
        sample_warnings = _sample_crispr_warnings(sample, metrics, thresholds)
        metrics["qc_status"] = "review" if sample_warnings else "pass"
        sample_report[sample] = metrics
        warnings.extend(sample_warnings)

    library_report = _crispr_library_report(library, k)
    if library_report.get("one_edit_collision_pairs", 0):
        warnings.append(
            {
                "code": "guide_collision",
                "severity": "review",
                "scope": "library",
                "message": f"{library_report['one_edit_collision_pairs']} guide pairs are duplicate or within one edit; one-edit rescue can create ambiguity.",
            }
        )
    if k > 1:
        warnings.append(
            {
                "code": "library_collision_audit_radius",
                "severity": "review",
                "scope": "library",
                "message": "CRISPR QC reports duplicate and one-edit guide collisions; run dotmatch audit-targets for full k>1 target-collision auditing.",
            }
        )
    if library_report.get("non_acgt_sequences", 0):
        warnings.append(
            {
                "code": "non_acgt_guides",
                "severity": "review",
                "scope": "library",
                "message": f"{library_report['non_acgt_sequences']} guide sequences contain symbols outside A/C/G/T.",
            }
        )

    sample_correlation_report = _pairwise_sample_correlations(sample_counts)
    for pair in sample_correlation_report:
        pearson = pair.get("pearson_log2_count_plus_1")
        if isinstance(pearson, float) and pearson < thresholds["pairwise_sample_pearson_min"]:
            warnings.append(
                {
                    "code": "low_pairwise_sample_correlation",
                    "severity": "review",
                    "scope": f"{pair['sample_a']}:{pair['sample_b']}",
                    "message": f"Pairwise sample log2(count+1) Pearson correlation is {pearson:.3f}.",
                }
            )

    return {
        "schema_version": 1,
        "assay": "crispr_count_qc",
        "status": "review" if warnings else "pass",
        "guide_count": guide_count,
        "sample_count": len(samples),
        "thresholds": thresholds,
        "samples": sample_report,
        "library": library_report,
        "sample_correlations": sample_correlation_report,
        "replicates": sample_correlation_report,
        "warnings": warnings,
        "interpretation": "QC flags are conservative diagnostics for guide counting and representation; downstream screen statistics should be run with MAGeCK or another screen-analysis method.",
    }


def _sample_crispr_warnings(sample: str, metrics: dict[str, object], thresholds: dict[str, float]) -> list[dict[str, object]]:
    checks = [
        ("low_assignment_rate", "assignment_rate", "<", thresholds["assignment_rate_min"], "assignment rate is low"),
        ("high_ambiguous_rate", "ambiguous_rate", ">", thresholds["ambiguous_rate_max"], "ambiguous rate is high"),
        ("high_no_match_rate", "no_match_rate", ">", thresholds["no_match_rate_max"], "no-match rate is high"),
        ("high_invalid_rate", "invalid_rate", ">", thresholds["invalid_rate_max"], "invalid extraction rate is high"),
        ("low_guide_coverage", "coverage_fraction", "<", thresholds["coverage_fraction_min"], "guide coverage is low"),
        ("high_zero_count_fraction", "zero_count_fraction", ">", thresholds["zero_count_fraction_max"], "zero-count guide fraction is high"),
        ("high_gini_index", "gini_index", ">", thresholds["gini_index_max"], "guide representation is skewed"),
        ("high_top_1pct_fraction", "top_1pct_fraction", ">", thresholds["top_1pct_fraction_max"], "top guides dominate assigned guide counts"),
    ]
    warnings: list[dict[str, object]] = []
    for code, key, op, threshold, message in checks:
        value = metrics.get(key)
        if not isinstance(value, (int, float)):
            continue
        triggered = value < threshold if op == "<" else value > threshold
        if triggered:
            warnings.append(
                {
                    "code": code,
                    "severity": "review",
                    "scope": sample,
                    "metric": key,
                    "value": value,
                    "threshold": threshold,
                    "message": f"{sample}: {message} ({key}={value:.4g}, threshold {op} {threshold:.4g}).",
                }
            )
    return warnings


def _crispr_library_report(library: Sequence[Target], k: int) -> dict[str, object]:
    if not library:
        return {
            "provided": False,
            "guide_count": 0,
            "duplicate_ids": 0,
            "duplicate_sequences": 0,
            "non_acgt_sequences": 0,
            "collision_radius_audited": min(k, 1) if k >= 1 else 0,
            "one_edit_collision_pairs": 0,
            "safe_for_audited_radius": None,
            "safe_for_k": None,
        }
    duplicate_ids = _duplicate_count([target.target_id for target in library])
    duplicate_sequences = _duplicate_count([target.seq for target in library])
    non_acgt = sum(1 for target in library if any(base not in DNA for base in target.seq))
    pairs = list(_near_target_pairs(library, min(k, 1))) if k >= 1 else list(_near_target_pairs(library, 0))
    guides_per_gene: dict[str, int] = defaultdict(int)
    for target in library:
        if target.gene:
            guides_per_gene[target.gene] += 1
    per_gene_counts = sorted(guides_per_gene.values())
    return {
        "provided": True,
        "guide_count": len(library),
        "gene_count": len(guides_per_gene),
        "guides_per_gene_min": per_gene_counts[0] if per_gene_counts else 0,
        "guides_per_gene_median": _median(per_gene_counts) if per_gene_counts else 0,
        "guides_per_gene_max": per_gene_counts[-1] if per_gene_counts else 0,
        "duplicate_ids": duplicate_ids,
        "duplicate_sequences": duplicate_sequences,
        "non_acgt_sequences": non_acgt,
        "collision_radius_audited": min(k, 1) if k >= 1 else 0,
        "one_edit_collision_pairs": len(pairs),
        "safe_for_audited_radius": len(pairs) == 0,
        "safe_for_k": (len(pairs) == 0) if k <= 1 else None,
    }


def _pairwise_sample_correlations(sample_counts: dict[str, list[int]]) -> list[dict[str, object]]:
    samples = list(sample_counts)
    pairs: list[dict[str, object]] = []
    for i, sample_a in enumerate(samples):
        for sample_b in samples[i + 1 :]:
            a = [math.log2(value + 1) for value in sample_counts[sample_a]]
            b = [math.log2(value + 1) for value in sample_counts[sample_b]]
            pairs.append(
                {
                    "sample_a": sample_a,
                    "sample_b": sample_b,
                    "pearson_log2_count_plus_1": _pearson(a, b),
                    "spearman_count": _pearson(_ranks(sample_counts[sample_a]), _ranks(sample_counts[sample_b])),
                }
            )
    return pairs


def _gini(values: Sequence[int]) -> float:
    if not values:
        return 0.0
    sorted_values = sorted(max(0, value) for value in values)
    total = sum(sorted_values)
    if total == 0:
        return 0.0
    n = len(sorted_values)
    weighted = sum((2 * i - n - 1) * value for i, value in enumerate(sorted_values, start=1))
    return weighted / (n * total)


def _pearson(a: Sequence[float], b: Sequence[float]) -> float | None:
    if len(a) != len(b) or len(a) < 2:
        return None
    mean_a = sum(a) / len(a)
    mean_b = sum(b) / len(b)
    da = [value - mean_a for value in a]
    db = [value - mean_b for value in b]
    denom_a = math.sqrt(sum(value * value for value in da))
    denom_b = math.sqrt(sum(value * value for value in db))
    if denom_a == 0 or denom_b == 0:
        return None
    return sum(x * y for x, y in zip(da, db)) / (denom_a * denom_b)


def _ranks(values: Sequence[int]) -> list[float]:
    ordered = sorted((value, i) for i, value in enumerate(values))
    ranks = [0.0] * len(values)
    pos = 0
    while pos < len(ordered):
        end = pos + 1
        while end < len(ordered) and ordered[end][0] == ordered[pos][0]:
            end += 1
        rank = (pos + 1 + end) / 2
        for _value, index in ordered[pos:end]:
            ranks[index] = rank
        pos = end
    return ranks


def _duplicate_count(values: Sequence[str]) -> int:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return len(duplicates)


def _median(values: Sequence[int]) -> float:
    if not values:
        return 0.0
    mid = len(values) // 2
    if len(values) % 2:
        return float(values[mid])
    return (values[mid - 1] + values[mid]) / 2


def _write_crispr_qc_summary_tsv(report: dict[str, object], path: str | Path) -> None:
    samples: dict[str, dict[str, object]] = report["samples"]  # type: ignore[assignment]
    columns = [
        "sample_id",
        "qc_status",
        "total_count",
        "coverage_fraction",
        "zero_count_fraction",
        "gini_index",
        "top_1pct_fraction",
        "assignment_rate",
        "ambiguous_rate",
        "no_match_rate",
        "invalid_rate",
    ]
    with _open_text(path, "wt") as fh:
        fh.write("\t".join(columns) + "\n")
        for sample_id, metrics in samples.items():
            fh.write("\t".join(str(metrics.get(column, sample_id if column == "sample_id" else "")) for column in columns) + "\n")


def _write_crispr_qc_html(report: dict[str, object], path: str | Path) -> None:
    samples: dict[str, dict[str, object]] = report["samples"]  # type: ignore[assignment]
    warnings: list[dict[str, object]] = report["warnings"]  # type: ignore[assignment]
    library: dict[str, object] = report["library"]  # type: ignore[assignment]
    rows = "\n".join(
        "<tr>"
        f"<td>{html.escape(sample_id)}</td>"
        f"<td>{html.escape(str(metrics.get('qc_status', '')))}</td>"
        f"<td>{_fmt(metrics.get('assignment_rate'))}</td>"
        f"<td>{_fmt(metrics.get('coverage_fraction'))}</td>"
        f"<td>{_fmt(metrics.get('zero_count_fraction'))}</td>"
        f"<td>{_fmt(metrics.get('gini_index'))}</td>"
        f"<td>{_fmt(metrics.get('top_1pct_fraction'))}</td>"
        "</tr>"
        for sample_id, metrics in samples.items()
    )
    warning_items = "\n".join(
        f"<li><strong>{html.escape(str(item.get('code', 'warning')))}</strong>: {html.escape(str(item.get('message', '')))}</li>"
        for item in warnings
    ) or "<li>No QC warnings.</li>"
    library_items = "\n".join(
        f"<li><strong>{html.escape(str(key))}</strong>: {html.escape(str(value))}</li>"
        for key, value in library.items()
    )
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>DotMatch CRISPR QC Report</title>
  <style>
    body {{ margin: 2rem; color: #101513; font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #dfe7e2; padding: 0.55rem; text-align: left; }}
    th {{ background: #eef6f2; }}
    code {{ background: #eef6f2; padding: 0.1rem 0.25rem; }}
  </style>
</head>
<body>
  <h1>DotMatch CRISPR QC Report</h1>
  <p>Status: <strong>{html.escape(str(report["status"]))}</strong></p>
  <p>{html.escape(str(report["interpretation"]))}</p>
  <h2>Sample QC</h2>
  <table>
    <thead><tr><th>Sample</th><th>Status</th><th>Assignment</th><th>Coverage</th><th>Zero guides</th><th>Gini</th><th>Top 1%</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
  <h2>Guide Library Audit</h2>
  <ul>{library_items}</ul>
  <h2>Warnings</h2>
  <ul>{warning_items}</ul>
</body>
</html>
"""
    with _open_text(path, "wt") as fh:
        fh.write(document)


def _fmt(value: object) -> str:
    if isinstance(value, float):
        return html.escape(f"{value:.4g}")
    return html.escape(str(value if value is not None else ""))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dotmatch",
        description="Exact known-target short-DNA assignment and counting.",
    )
    parser.add_argument("--version", action="version", version=f"dotmatch {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    count = sub.add_parser("count", help="stream FASTQ/FASTQ.gz and emit target count tables")
    count.add_argument("--targets", required=True, help="TSV with target_id, target_seq, optional gene")
    count.add_argument("--reads", required=True, help="FASTQ or FASTQ.gz input")
    count.add_argument("--target-start", type=int, default=0)
    count.add_argument("--target-length", type=int, required=True)
    count.add_argument("--k", type=int, default=1)
    count.add_argument("--out", required=True, help="counts TSV output")
    count.add_argument("--assignments", help="optional per-read assignments TSV")
    count.add_argument("--summary", help="optional summary JSON output")
    count.add_argument("--ambiguous", choices=["discard", "report"], default="discard")
    count.add_argument("--batch-size", type=int, default=4096)
    count.set_defaults(func=command_count)

    audit = sub.add_parser("audit-targets", help="report target pairs that make k-edit correction ambiguous")
    audit.add_argument("--targets", required=True)
    audit.add_argument("--k", type=int, default=1)
    audit.add_argument("--out", help="optional nearby-pairs TSV")
    audit.set_defaults(func=command_audit_targets)

    validate = sub.add_parser("validate", help="compare indexed assignment against native exhaustive scan")
    validate.add_argument("--targets", required=True)
    validate.add_argument("--reads", required=True)
    validate.add_argument("--target-start", type=int, default=0)
    validate.add_argument("--target-length", type=int, required=True)
    validate.add_argument("--k", type=int, default=1)
    validate.add_argument("--sample", type=int, default=100000)
    validate.add_argument("--batch-size", type=int, default=4096)
    validate.add_argument("--show-mismatches", action="store_true")
    validate.set_defaults(func=command_validate)

    crispr_qc = sub.add_parser("crispr-qc", help="evaluate CRISPR guide-count QC and representation metrics")
    _add_crispr_qc_args(crispr_qc)
    crispr_qc.set_defaults(func=command_crispr_qc)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    raw_args = list(sys.argv[1:] if argv is None else argv)
    if raw_args and raw_args[0] == "assay":
        return command_assay(raw_args[1:])
    if raw_args and raw_args[0] == "crispr":
        return command_crispr_namespace(raw_args[1:])
    if raw_args and raw_args[0] == "barcode":
        return command_barcode_namespace(raw_args[1:])
    if raw_args and raw_args[0] == "panel":
        from .panel import command_panel_namespace

        return command_panel_namespace(raw_args[1:])
    if raw_args and raw_args[0] == "crispr-qc":
        parser = build_parser()
        args = parser.parse_args(raw_args)
        try:
            return int(args.func(args))
        except BrokenPipeError:
            return 1
        except Exception as exc:
            print(f"dotmatch: {exc}", file=sys.stderr)
            return 2
    if not os.environ.get("DOTMATCH_PYTHON_NO_DELEGATE"):
        try:
            return run_native_cli(raw_args)
        except FileNotFoundError:
            pass

    parser = build_parser()
    args = parser.parse_args(raw_args)
    if getattr(args, "k", 0) < 0:
        parser.error("--k must be non-negative")
    try:
        return int(args.func(args))
    except BrokenPipeError:
        return 1
    except Exception as exc:
        print(f"dotmatch: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
