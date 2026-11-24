"""Declarative multi-segment read grammar for Perturb-seq, DEL, and dual-target FASTQs.

Each read (or R1/R2 pair) is decomposed into named fixed-window segments.
Target segments use the native DotMatch Matcher with explicit four-way outcomes
(unique, ambiguous, none, invalid). An optional UMI segment deduplicates
molecules within an explicit scope while surfacing UMI-to-target conflicts
rather than collapsing them silently.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal, Sequence

from ._types import AssignmentPolicy, Metric
from .core import (
    MATCH_AMBIGUOUS,
    MATCH_INVALID,
    MATCH_NONE,
    MATCH_UNIQUE,
    FastqRecord,
    Matcher,
    iter_fastq,
)
from .target_io import read_target_table

SegmentKind = Literal["target", "umi"]
ReadMate = Literal[1, 2]

SCHEMA_VERSION = "dotmatch.segments.v1"
_VALID_UMI_BASES = frozenset("ACGTacgt")


@dataclass(frozen=True)
class SegmentSpec:
    name: str
    kind: SegmentKind
    start: int
    length: int
    read: ReadMate = 1
    targets: str | Path | Sequence[tuple[str, str]] | None = None
    k: int = 0
    metric: Metric = "hamming"
    policy: AssignmentPolicy = "radius"
    scope: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.name or not self.name.strip():
            raise ValueError("segment name must be non-empty")
        if self.start < 0:
            raise ValueError(f"segment {self.name!r}: start must be non-negative")
        if self.length < 1:
            raise ValueError(f"segment {self.name!r}: length must be positive")
        if self.read not in (1, 2):
            raise ValueError(f"segment {self.name!r}: read must be 1 or 2")
        if self.kind == "target":
            if self.targets is None:
                raise ValueError(f"target segment {self.name!r} requires targets")
            if self.scope:
                raise ValueError(f"target segment {self.name!r} cannot declare umi scope")
            max_k = 3 if self.metric == "hamming" else (0 if self.metric == "exact" else 2)
            if self.k < 0 or self.k > max_k:
                raise ValueError(f"segment {self.name!r}: k must be in 0..{max_k} for metric {self.metric!r}")
        elif self.kind == "umi":
            if self.targets is not None:
                raise ValueError(f"umi segment {self.name!r} cannot specify targets")
        else:
            raise ValueError(f"unsupported segment kind: {self.kind!r}")


@dataclass(frozen=True)
class ReadGrammar:
    segments: tuple[SegmentSpec, ...]

    def __post_init__(self) -> None:
        if not self.segments:
            raise ValueError("ReadGrammar requires at least one segment")
        names = [seg.name for seg in self.segments]
        if len(set(names)) != len(names):
            raise ValueError("segment names in ReadGrammar must be unique")
        target_names = [seg.name for seg in self.segments if seg.kind == "target"]
        if not target_names:
            raise ValueError("ReadGrammar requires at least one target segment")
        umi_segments = [seg for seg in self.segments if seg.kind == "umi"]
        if len(umi_segments) > 1:
            raise ValueError("ReadGrammar supports at most one umi segment")
        if umi_segments:
            unknown_scope = [name for name in umi_segments[0].scope if name not in target_names]
            if unknown_scope:
                raise ValueError(f"umi scope references unknown target segment(s): {', '.join(unknown_scope)}")

    @property
    def target_segments(self) -> tuple[SegmentSpec, ...]:
        return tuple(seg for seg in self.segments if seg.kind == "target")

    @property
    def umi_segment(self) -> SegmentSpec | None:
        for seg in self.segments:
            if seg.kind == "umi":
                return seg
        return None

    @property
    def requires_r2(self) -> bool:
        return any(seg.read == 2 for seg in self.segments)


@dataclass(frozen=True)
class SegmentOutcomeSummary:
    name: str
    kind: SegmentKind
    unique: int
    ambiguous: int
    none: int
    invalid: int


@dataclass(frozen=True)
class DecodedSegmentsReport:
    schema_version: str
    total_reads: int
    joint_unique_reads: int
    joint_ambiguous_reads: int
    joint_unmatched_reads: int
    joint_invalid_reads: int
    umi_dedup_enabled: bool
    unique_molecules: int
    umi_conflict_molecules: int
    umi_conflict_reads: int
    segment_summaries: tuple[SegmentOutcomeSummary, ...]
    artifacts: dict[str, str]

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "total_reads": self.total_reads,
            "joint_unique_reads": self.joint_unique_reads,
            "joint_ambiguous_reads": self.joint_ambiguous_reads,
            "joint_unmatched_reads": self.joint_unmatched_reads,
            "joint_invalid_reads": self.joint_invalid_reads,
            "umi_dedup_enabled": self.umi_dedup_enabled,
            "unique_molecules": self.unique_molecules,
            "umi_conflict_molecules": self.umi_conflict_molecules,
            "umi_conflict_reads": self.umi_conflict_reads,
            "segment_summaries": [asdict(item) for item in self.segment_summaries],
            "artifacts": dict(self.artifacts),
        }


def _load_segment_targets(targets: str | Path | Sequence[tuple[str, str]]) -> tuple[list[str], list[str]]:
    if isinstance(targets, (str, Path)):
        rows = read_target_table(Path(targets))
        return [row.target_id for row in rows], [row.sequence for row in rows]
    ids = [str(tid) for tid, _seq in targets]
    seqs = [str(seq) for _tid, seq in targets]
    return ids, seqs


def _assign_batch(matcher: Matcher, seg: SegmentSpec, windows: list[str]) -> list[tuple[int, int]]:
    if not windows:
        return []
    if seg.metric == "exact":
        results = matcher.assign_exact(windows, policy=seg.policy)
    elif seg.metric == "hamming":
        results = matcher.assign_hamming(windows, k=seg.k, policy=seg.policy)
    else:
        results = matcher.assign(windows, k=seg.k, policy=seg.policy)
    return [(res.status, res.target_index) for res in results]


def _extract_window(seq: str, start: int, length: int) -> str | None:
    if start + length > len(seq):
        return None
    return seq[start : start + length]


def _iter_paired_records(
    reads_r1: Path, reads_r2: Path | None
) -> Any:
    if reads_r2 is None:
        for rec in iter_fastq(reads_r1):
            yield rec, None
        return
    it1 = iter(iter_fastq(reads_r1))
    it2 = iter(iter_fastq(reads_r2))
    while True:
        rec1: FastqRecord | None = next(it1, None)
        rec2: FastqRecord | None = next(it2, None)
        if rec1 is None and rec2 is None:
            break
        if rec1 is None or rec2 is None:
            raise ValueError("R1 and R2 FASTQ files have different record counts")
        yield rec1, rec2


def decode_segments(
    *,
    grammar: ReadGrammar,
    reads: str | Path,
    out_dir: str | Path,
    reads_r2: str | Path | None = None,
    batch_size: int = 4096,
) -> DecodedSegmentsReport:
    """Decode multi-segment reads into joint target counts and UMI molecule counts."""
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    if grammar.requires_r2 and reads_r2 is None:
        raise ValueError("ReadGrammar references read=2 segments; reads_r2 is required")
    reads_path = Path(reads)
    reads_r2_path = Path(reads_r2) if reads_r2 is not None else None
    if not reads_path.is_file():
        raise ValueError(f"FASTQ is not a regular file: {reads_path}")
    if reads_r2_path is not None and not reads_r2_path.is_file():
        raise ValueError(f"R2 FASTQ is not a regular file: {reads_r2_path}")

    output = Path(out_dir)
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError(f"out_dir must be a new or empty directory: {output}")
    output.mkdir(parents=True, exist_ok=True)

    target_segs = grammar.target_segments
    umi_seg = grammar.umi_segment
    target_names = [seg.name for seg in target_segs]
    scope_indices = tuple(target_names.index(name) for name in umi_seg.scope) if umi_seg is not None else ()

    matchers: list[Matcher] = []
    target_id_tables: list[list[str]] = []
    seg_counters: dict[str, Counter[str]] = {seg.name: Counter() for seg in grammar.segments}

    total_reads = 0
    joint_unique = 0
    joint_ambiguous = 0
    joint_unmatched = 0
    joint_invalid = 0

    tuple_read_counts: Counter[tuple[str, ...]] = Counter()
    # Maps (scope_values, umi) -> Counter of full target_tuples observed with that UMI
    umi_observations: dict[tuple[tuple[str, ...], str], Counter[tuple[str, ...]]] = defaultdict(Counter)

    try:
        for seg in target_segs:
            assert seg.targets is not None
            ids, seqs = _load_segment_targets(seg.targets)
            target_id_tables.append(ids)
            matchers.append(Matcher(seqs))

        batch: list[tuple[FastqRecord, FastqRecord | None]] = []

        def process_batch(records: list[tuple[FastqRecord, FastqRecord | None]]) -> None:
            nonlocal total_reads, joint_unique, joint_ambiguous, joint_unmatched, joint_invalid
            n = len(records)
            total_reads += n

            seg_statuses: list[list[int]] = []
            seg_target_ids: list[list[str]] = []

            for seg_idx, seg in enumerate(target_segs):
                valid_indices: list[int] = []
                windows: list[str] = []
                statuses: list[int] = [int(MATCH_INVALID)] * n
                resolved_ids = [""] * n
                for i, (rec1, rec2) in enumerate(records):
                    seq = rec1.seq if seg.read == 1 else (rec2.seq if rec2 is not None else "")
                    win = _extract_window(seq, seg.start, seg.length)
                    if win is None:
                        seg_counters[seg.name]["invalid"] += 1
                    else:
                        valid_indices.append(i)
                        windows.append(win)
                assigned = _assign_batch(matchers[seg_idx], seg, windows)
                ids_table = target_id_tables[seg_idx]
                for idx, (status, target_idx) in zip(valid_indices, assigned):
                    statuses[idx] = status
                    if status == MATCH_UNIQUE and 0 <= target_idx < len(ids_table):
                        resolved_ids[idx] = ids_table[target_idx]
                        seg_counters[seg.name]["unique"] += 1
                    elif status == MATCH_AMBIGUOUS:
                        seg_counters[seg.name]["ambiguous"] += 1
                    elif status == MATCH_NONE:
                        seg_counters[seg.name]["none"] += 1
                    else:
                        seg_counters[seg.name]["invalid"] += 1
                seg_statuses.append(statuses)
                seg_target_ids.append(resolved_ids)

            umi_values: list[str | None] = [None] * n
            if umi_seg is not None:
                for i, (rec1, rec2) in enumerate(records):
                    seq = rec1.seq if umi_seg.read == 1 else (rec2.seq if rec2 is not None else "")
                    win = _extract_window(seq, umi_seg.start, umi_seg.length)
                    if win is not None and all(ch in _VALID_UMI_BASES for ch in win):
                        umi_values[i] = win.upper()
                        seg_counters[umi_seg.name]["unique"] += 1
                    else:
                        seg_counters[umi_seg.name]["invalid"] += 1

            for i in range(n):
                statuses_i = [st[i] for st in seg_statuses]
                umi_bad = umi_seg is not None and umi_values[i] is None
                if umi_bad or any(st == MATCH_INVALID for st in statuses_i):
                    joint_invalid += 1
                elif any(st == MATCH_NONE for st in statuses_i):
                    joint_unmatched += 1
                elif any(st == MATCH_AMBIGUOUS for st in statuses_i):
                    joint_ambiguous += 1
                else:
                    joint_unique += 1
                    target_tuple = tuple(ids[i] for ids in seg_target_ids)
                    tuple_read_counts[target_tuple] += 1
                    if umi_seg is not None:
                        umi_val = umi_values[i]
                        assert umi_val is not None
                        scope_key = tuple(target_tuple[idx] for idx in scope_indices)
                        umi_observations[(scope_key, umi_val)][target_tuple] += 1

        for pair in _iter_paired_records(reads_path, reads_r2_path):
            batch.append(pair)
            if len(batch) >= batch_size:
                process_batch(batch)
                batch = []
        if batch:
            process_batch(batch)
    finally:
        for matcher in matchers:
            matcher.close()

    tuple_umi_counts: Counter[tuple[str, ...]] = Counter()
    unique_molecules = 0
    umi_conflict_molecules = 0
    umi_conflict_reads = 0
    conflict_rows: list[tuple[tuple[str, ...], str, dict[tuple[str, ...], int]]] = []

    if umi_seg is not None:
        for (scope_key, umi_val), target_counts in sorted(umi_observations.items()):
            if len(target_counts) == 1:
                only_tuple = next(iter(target_counts))
                tuple_umi_counts[only_tuple] += 1
                unique_molecules += 1
            else:
                umi_conflict_molecules += 1
                umi_conflict_reads += sum(target_counts.values())
                conflict_rows.append((scope_key, umi_val, dict(target_counts)))

    counts_path = output / "segment_counts.tsv"
    header = [*target_names, "read_count"] + (["umi_count"] if umi_seg is not None else [])
    with counts_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for target_tuple in sorted(tuple_read_counts):
            row: list[Any] = [*target_tuple, tuple_read_counts[target_tuple]]
            if umi_seg is not None:
                row.append(tuple_umi_counts[target_tuple])
            writer.writerow(row)

    artifacts: dict[str, str] = {"segment_counts": str(counts_path)}
    if umi_seg is not None:
        conflicts_path = output / "umi_conflicts.tsv"
        with conflicts_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([*umi_seg.scope, umi_seg.name, "conflicting_targets", "read_count"])
            for scope_key, umi_val, conflict_counts in conflict_rows:
                targets_str = ",".join(
                    f"{':'.join(t)}={c}" for t, c in sorted(conflict_counts.items())
                )
                writer.writerow([*scope_key, umi_val, targets_str, sum(conflict_counts.values())])
        artifacts["umi_conflicts"] = str(conflicts_path)

    summaries = tuple(
        SegmentOutcomeSummary(
            name=seg.name,
            kind=seg.kind,
            unique=seg_counters[seg.name]["unique"],
            ambiguous=seg_counters[seg.name]["ambiguous"],
            none=seg_counters[seg.name]["none"],
            invalid=seg_counters[seg.name]["invalid"],
        )
        for seg in grammar.segments
    )

    summary_path = output / "segments_summary.json"
    artifacts["segments_summary"] = str(summary_path)
    report = DecodedSegmentsReport(
        schema_version=SCHEMA_VERSION,
        total_reads=total_reads,
        joint_unique_reads=joint_unique,
        joint_ambiguous_reads=joint_ambiguous,
        joint_unmatched_reads=joint_unmatched,
        joint_invalid_reads=joint_invalid,
        umi_dedup_enabled=umi_seg is not None,
        unique_molecules=unique_molecules,
        umi_conflict_molecules=umi_conflict_molecules,
        umi_conflict_reads=umi_conflict_reads,
        segment_summaries=summaries,
        artifacts=artifacts,
    )
    summary_path.write_text(json.dumps(report.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report
