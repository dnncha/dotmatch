from __future__ import annotations

import argparse
import gzip
import json
import os
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
from .native import run_native_cli


DNA = "ACGT"


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
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    raw_args = list(sys.argv[1:] if argv is None else argv)
    if raw_args and raw_args[0] == "assay":
        return command_assay(raw_args[1:])
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
