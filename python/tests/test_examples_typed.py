"""Example tests for the typed DotMatch API.

Every example is executable, uses the same public entry points a user would,
and pins the result with an inline snapshot. Run with ``UPDATE_EXPECT=1`` to
re-record a deliberate behaviour change.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import dotmatch as dm
from dotmatch import assignment_summary, load_targets
from dotmatch.panel import (
    AssignmentConfig,
    PanelConstraints,
    PanelRecord,
    check_panel,
    hamming_distance,
    levenshtein_distance,
    reverse_complement,
)
from expect import expect_snapshot


def _json(value: object) -> str:
    return json.dumps(value, indent=2, sort_keys=True, default=str)


def test_example_scan_assignment_statuses() -> None:
    reads = ["ACGT", "ACGC", "TTTT", "ACGG"]
    targets = ["ACGT", "AGGT", "ACGA"]

    results = dm.assign(reads, targets, k=1)
    observed = [
        {
            "read": read,
            "status": dm.status_name(result.status),
            "status_value": int(result.status),
            "target_index": result.target_index,
            "best_distance": result.best_distance,
            "match_count": result.match_count,
        }
        for read, result in zip(reads, results)
    ]

    expect_snapshot(_json(observed), name="scan_assignment_statuses")


def test_example_radius_and_best_policy_are_explicit() -> None:
    reads = ["ACGT"]
    targets = ["ACGT", "AGGT", "ACGA"]

    radius = dm.assign(reads, targets, k=1, policy="radius")[0]
    best = dm.assign(reads, targets, k=1, policy="best")[0]

    expect_snapshot(
        _json(
            {
                "radius": {"status": dm.status_name(radius.status), "match_count": radius.match_count},
                "best": {"status": dm.status_name(best.status), "match_count": best.match_count},
            }
        ),
        name="radius_vs_best_policy",
    )


def test_example_hamming_and_exact_paths() -> None:
    targets = ["ACGT", "TTTT"]
    reads = ["ACGT", "ACGA", "ACGTTA"]

    hamming = [int(result.status) for result in dm.assign_hamming(reads, targets, k=1)]
    exact = [int(result.status) for result in dm.assign_exact(reads, targets)]

    expect_snapshot(
        _json({"hamming_status": hamming, "exact_status": exact}),
        name="hamming_and_exact_statuses",
    )


def test_example_stats_are_candidate_counts() -> None:
    matcher = dm.Matcher(["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT"])
    results, stats = matcher.assign_with_stats(["AAAAAAAT", "CCCCCCCA"], k=1)
    matcher.close()

    expect_snapshot(
        _json(
            {
                "statuses": [dm.status_name(result.status) for result in results],
                "candidates_considered": stats.candidates_considered,
                "candidates_verified": stats.candidates_verified,
            }
        ),
        name="assignment_stats",
    )


def test_example_posterior_uses_quality() -> None:
    result = dm.assign_posterior("ACGT", ["ACGT", "AGGT"], "IIII")
    low_quality = dm.assign_posterior("ACGT", ["ACGT", "AGGT"], "!!!!")

    expect_snapshot(
        _json(
            {
                "posterior": round(result.posterior, 6),
                "status": dm.status_name(result.status),
                "low_quality_status": dm.status_name(low_quality.status),
            }
        ),
        name="posterior_quality",
    )


def test_example_load_targets_is_typed(tmp_path: Path) -> None:
    panel = tmp_path / "targets.tsv"
    panel.write_text("guide_id\tguide_seq\nG1\tACGT\nG2\tAGGT\n", encoding="utf-8")

    targets = load_targets(panel)

    assert targets[0].target_id == "G1"
    assert targets[0].sequence == "ACGT"
    expect_snapshot(_json([target._asdict() for target in targets]), name="typed_targets")


def test_example_stream_fastq_end_to_end(tmp_path: Path) -> None:
    fastq = tmp_path / "reads.fastq"
    fastq.write_text(
        "@r1\nTTTTTTTT\n+\nIIIIIIII\n"
        "@r2\nACGTACGT\n+\nIIIIIIII\n"
        "@r3\nTTTT\n+\nIIII\n",
        encoding="ascii",
    )
    out = tmp_path / "assignments.tsv"

    streamed = list(
        dm.stream_assign(
            fastq,
            ["ACGTACGT", "TTTTTTTT", "ACGTACGA"],
            metric="hamming",
            k=1,
        )
    )
    summary = dm.write_assignments_tsv(streamed, out)

    expect_snapshot(_json(summary), name="stream_summary")
    expect_snapshot(out.read_text(encoding="utf-8"), name="stream_assignments")
    assert assignment_summary(streamed) == summary


def test_example_panel_certificate(tmp_path: Path) -> None:
    records = [PanelRecord("BC001", "ACGTACGT"), PanelRecord("BC002", "ACGTACGA")]
    certificate = check_panel(
        records,
        AssignmentConfig(metric="hamming", k=1, ambiguity_policy="radius"),
        PanelConstraints(length=8, min_hamming_distance=5, min_levenshtein_distance=4),
        tmp_path,
        source_panel=Path("barcodes.tsv"),
    )

    selected = {
        "status": certificate["status"],
        "n_barcodes": certificate["n_barcodes"],
        "lengths": certificate["lengths"],
        "assignment_metric": certificate["assignment_metric"],
        "safe_for_k1_hamming": certificate["safe_for_k1_hamming"],
        "minimum_hamming_distance": certificate["minimum_hamming_distance"],
        "ambiguous_error_spheres": certificate["ambiguous_error_spheres"],
        "certified_dotmatch_command": certificate["certified_dotmatch_command"],
    }
    expect_snapshot(_json(selected), name="panel_certificate")
    assert (tmp_path / "panel_summary.json").exists()


def test_example_distance_primitives() -> None:
    assert hamming_distance("ACGT", "AGGT") == 1
    assert levenshtein_distance("ACGT", "ACGGT") == 1
    assert reverse_complement("AACG") == "CGTT"


def test_example_typed_status_remains_int_compatible() -> None:
    result = dm.assign(["ACGT"], ["ACGT"], k=0)[0]

    assert isinstance(result.status, dm.MatchStatus)
    assert result.status == 1
    assert int(result.status) == 1
    assert json.loads(json.dumps(result.status)) == 1
    assert dm.status_name(result.status) == "unique"


def test_example_script_runs_from_a_clean_process(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[2]
    script = root / "examples" / "python" / "typed_quickstart.py"
    extension = "dylib" if sys.platform == "darwin" else "so"
    env = {
        **os.environ,
        "PYTHONPATH": str(root / "python"),
        "DOTMATCH_LIB": str(root / f"libdotmatch.{extension}"),
    }

    completed = subprocess.run(
        [sys.executable, str(script)],
        cwd=tmp_path,
        env=env,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert completed.returncode == 0, completed.stderr
    payload = json.loads(completed.stdout)
    assert payload["statuses"] == ["ambiguous", "ambiguous", "none"]
    assert payload["typed"]["summary"]["total_reads"] == 3
