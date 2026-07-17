"""Streaming, bounded-memory assignment quality monitoring."""

from __future__ import annotations

import json
import math
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping, TextIO


STATUSES = {"unique", "ambiguous", "none", "invalid"}


@dataclass(frozen=True)
class WatchThresholds:
    min_assignment_rate: float = 0.80
    max_ambiguous_rate: float = 0.05
    max_unmatched_rate: float = 0.15
    max_invalid_rate: float = 0.02
    min_reads: int = 1000


@dataclass(frozen=True)
class RunSnapshot:
    total_reads: int
    unique: int
    ambiguous: int
    unmatched: int
    invalid: int
    assignment_rate: float
    ambiguous_rate: float
    unmatched_rate: float
    invalid_rate: float
    assignment_rate_interval95: tuple[float, float]
    distinct_targets: int
    decision: str
    findings: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


class SequentialMonitor:
    def __init__(self, thresholds: WatchThresholds | None = None) -> None:
        self.thresholds = thresholds or WatchThresholds()
        if self.thresholds.min_reads <= 0:
            raise ValueError("min_reads must be positive")
        self.counts: Counter[str] = Counter()
        self.targets: Counter[str] = Counter()

    def update(self, status: str, target: str | None = None) -> None:
        if status not in STATUSES:
            raise ValueError(f"unsupported assignment status: {status}")
        self.counts[status] += 1
        if status == "unique" and target:
            self.targets[target] += 1

    def snapshot(self) -> RunSnapshot:
        total = sum(self.counts.values())
        unique = self.counts["unique"]
        ambiguous = self.counts["ambiguous"]
        unmatched = self.counts["none"]
        invalid = self.counts["invalid"]
        assignment_rate = unique / total if total else 0.0
        ambiguous_rate = ambiguous / total if total else 0.0
        unmatched_rate = unmatched / total if total else 0.0
        invalid_rate = invalid / total if total else 0.0
        findings: list[str] = []
        if total >= self.thresholds.min_reads:
            if assignment_rate < self.thresholds.min_assignment_rate:
                findings.append("assignment_rate_below_min")
            if ambiguous_rate > self.thresholds.max_ambiguous_rate:
                findings.append("ambiguous_rate_above_max")
            if unmatched_rate > self.thresholds.max_unmatched_rate:
                findings.append("unmatched_rate_above_max")
            if invalid_rate > self.thresholds.max_invalid_rate:
                findings.append("invalid_rate_above_max")
            decision = "review" if findings else "on_track"
        else:
            decision = "insufficient_data"
        return RunSnapshot(
            total_reads=total,
            unique=unique,
            ambiguous=ambiguous,
            unmatched=unmatched,
            invalid=invalid,
            assignment_rate=assignment_rate,
            ambiguous_rate=ambiguous_rate,
            unmatched_rate=unmatched_rate,
            invalid_rate=invalid_rate,
            assignment_rate_interval95=_wilson_interval(unique, total),
            distinct_targets=len(self.targets),
            decision=decision,
            findings=tuple(findings),
        )


def monitor_events(
    events: Iterable[Mapping[str, object]],
    *,
    every: int,
    thresholds: WatchThresholds | None = None,
) -> Iterator[RunSnapshot]:
    if every <= 0:
        raise ValueError("every must be positive")
    monitor = SequentialMonitor(thresholds)
    emitted_at = 0
    for event in events:
        status = event.get("status")
        if not isinstance(status, str):
            raise ValueError("assignment event requires a string status")
        target_value = event.get("target")
        target = target_value if isinstance(target_value, str) and target_value else None
        monitor.update(status, target)
        total = sum(monitor.counts.values())
        if total % every == 0:
            emitted_at = total
            yield monitor.snapshot()
    total = sum(monitor.counts.values())
    if total != emitted_at:
        yield monitor.snapshot()


def iter_jsonl(handle: TextIO) -> Iterator[Mapping[str, object]]:
    for line_number, line in enumerate(handle, start=1):
        if not line.strip():
            continue
        try:
            event = json.loads(line)
        except json.JSONDecodeError as exc:
            raise ValueError(f"invalid JSON on line {line_number}: {exc}") from exc
        if not isinstance(event, dict):
            raise ValueError(f"line {line_number} must contain a JSON object")
        yield event


def watch_jsonl(
    source: str | Path,
    output: str | Path,
    *,
    every: int = 100000,
    thresholds: WatchThresholds | None = None,
) -> RunSnapshot | None:
    input_handle = sys.stdin if str(source) == "-" else Path(source).open("rt", encoding="utf-8")
    output_handle = sys.stdout if str(output) == "-" else Path(output).open("wt", encoding="utf-8")
    latest: RunSnapshot | None = None
    try:
        for snapshot in monitor_events(iter_jsonl(input_handle), every=every, thresholds=thresholds):
            latest = snapshot
            output_handle.write(json.dumps(snapshot.to_dict(), sort_keys=True) + "\n")
            output_handle.flush()
    finally:
        if input_handle is not sys.stdin:
            input_handle.close()
        if output_handle is not sys.stdout:
            output_handle.close()
    return latest


def _wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total == 0:
        return (0.0, 1.0)
    proportion = successes / total
    denominator = 1.0 + z * z / total
    centre = proportion + z * z / (2.0 * total)
    margin = z * math.sqrt(
        proportion * (1.0 - proportion) / total + z * z / (4.0 * total * total)
    )
    return (
        max(0.0, (centre - margin) / denominator),
        min(1.0, (centre + margin) / denominator),
    )
