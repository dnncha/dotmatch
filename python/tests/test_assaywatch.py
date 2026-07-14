from __future__ import annotations

import json
from pathlib import Path

import pytest

from assaycode import cli
from dotmatch.assaywatch import (
    SequentialMonitor,
    WatchThresholds,
    monitor_events,
    watch_jsonl,
)


def test_monitor_waits_for_minimum_evidence() -> None:
    monitor = SequentialMonitor(WatchThresholds(min_reads=3))
    monitor.update("unique", "g1")
    monitor.update("none")

    snapshot = monitor.snapshot()
    assert snapshot.decision == "insufficient_data"
    assert snapshot.total_reads == 2
    assert snapshot.assignment_rate_interval95[0] < snapshot.assignment_rate
    assert snapshot.assignment_rate_interval95[1] > snapshot.assignment_rate


def test_monitor_flags_actionable_rate_failures() -> None:
    monitor = SequentialMonitor(
        WatchThresholds(
            min_reads=4,
            min_assignment_rate=0.8,
            max_ambiguous_rate=0.1,
            max_unmatched_rate=0.1,
            max_invalid_rate=0.1,
        )
    )
    for status in ["unique", "ambiguous", "none", "invalid"]:
        monitor.update(status, "g1" if status == "unique" else None)

    snapshot = monitor.snapshot()
    assert snapshot.decision == "review"
    assert set(snapshot.findings) == {
        "assignment_rate_below_min",
        "ambiguous_rate_above_max",
        "unmatched_rate_above_max",
        "invalid_rate_above_max",
    }


def test_monitor_reports_on_track_and_target_coverage() -> None:
    monitor = SequentialMonitor(WatchThresholds(min_reads=4))
    for target in ["g1", "g1", "g2", "g2"]:
        monitor.update("unique", target)

    snapshot = monitor.snapshot()
    assert snapshot.decision == "on_track"
    assert snapshot.distinct_targets == 2
    assert snapshot.assignment_rate == 1.0


def test_monitor_events_emits_intervals_and_final_partial_batch() -> None:
    snapshots = list(
        monitor_events(
            [
                {"status": "unique", "target": "g1"},
                {"status": "unique", "target": "g1"},
                {"status": "none"},
            ],
            every=2,
            thresholds=WatchThresholds(min_reads=2),
        )
    )

    assert [snapshot.total_reads for snapshot in snapshots] == [2, 3]


def test_watch_jsonl_writes_machine_readable_snapshots(tmp_path: Path) -> None:
    events = tmp_path / "events.jsonl"
    events.write_text(
        "\n".join(
            [
                json.dumps({"status": "unique", "target": "g1"}),
                json.dumps({"status": "none"}),
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output = tmp_path / "snapshots.jsonl"

    latest = watch_jsonl(
        events,
        output,
        every=1,
        thresholds=WatchThresholds(min_reads=2),
    )

    assert latest is not None
    rows = [json.loads(line) for line in output.read_text(encoding="utf-8").splitlines()]
    assert [row["total_reads"] for row in rows] == [1, 2]
    assert rows[-1]["decision"] == "review"


def test_watch_rejects_unknown_status(tmp_path: Path) -> None:
    events = tmp_path / "events.jsonl"
    events.write_text('{"status":"fabricated"}\n', encoding="utf-8")

    with pytest.raises(ValueError, match="unsupported assignment status"):
        watch_jsonl(events, tmp_path / "out.jsonl", every=1)


def test_assaycode_watch_exit_gate(tmp_path: Path) -> None:
    events = tmp_path / "events.jsonl"
    events.write_text(
        '{"status":"none"}\n{"status":"none"}\n',
        encoding="utf-8",
    )
    output = tmp_path / "snapshots.jsonl"

    rc = cli.main(
        [
            "watch",
            str(events),
            "--out",
            str(output),
            "--every",
            "1",
            "--min-reads",
            "2",
            "--fail-on-review",
        ]
    )

    assert rc == 1
    assert json.loads(output.read_text(encoding="utf-8").splitlines()[-1])["decision"] == "review"
