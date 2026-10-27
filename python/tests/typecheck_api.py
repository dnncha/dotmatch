"""Static-only API contract tests.

This file is analyzed by mypy and pyright; it is not collected by pytest.
Every ``assert_type`` call is a no-op at runtime and a hard check in the type
checker. The deliberately invalid calls use ``type: ignore[arg-type]`` so that
mypy's ``warn_unused_ignores`` proves the checker rejects them.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING, Iterator

if sys.version_info >= (3, 11):
    from typing import assert_type
else:
    from typing_extensions import assert_type

import dotmatch as dm
from dotmatch import (
    AssignmentSummary,
    FastqRecord,
    Matcher,
    MatchResult,
    MatchStatus,
    PosteriorAssignment,
    StreamAssignment,
    TargetRow,
)

assert_type(dm.MATCH_UNIQUE, MatchStatus)

results = dm.assign(["ACGT"], ["ACGT"], k=1, policy="best")
assert_type(results, list[MatchResult])
assert_type(results[0].status, MatchStatus)
assert_type(dm.status_name(results[0].status), str)

hamming = dm.assign_hamming(["ACGT"], ["ACGT"], k=0)
assert_type(hamming, list[MatchResult])

exact = dm.assign_exact(["ACGT"], ["ACGT"], ascii_fold=True)
assert_type(exact, list[MatchResult])

with Matcher(["ACGT"]) as matcher:
    with_stats = matcher.assign_with_stats(["ACGT"], k=1, policy="radius")
    assert_type(with_stats, tuple[list[MatchResult], dm.AssignmentStats])
    assert_type(with_stats[1].candidates_verified, int)

posterior = dm.assign_posterior("ACGT", ["ACGT", "AGGT"], "IIII")
assert_type(posterior, PosteriorAssignment)
assert_type(posterior.status, MatchStatus)

distance: int = dm.distance("ACGT", "AGGT")
assert_type(distance, int)
assert_type(dm.distance_leq("ACGT", "AGGT", 1), bool)

path = Path("targets.tsv")
assert_type(dm.load_targets(path), list[TargetRow])
assert_type(dm.iter_fastq(path), Iterator[FastqRecord])
assert_type(dm.targets_from_dataframe(object()), list[TargetRow])

summary = dm.assignment_summary([])
assert_type(summary, AssignmentSummary)
assert_type(summary["total_reads"], int)
assert_type(summary["assignment_rate"], float)

with_targets = dm.stream_assign(path, ["ACGT"], metric="hamming", policy="best")
assert_type(with_targets, Iterator[StreamAssignment])

# --- Negative type contracts -------------------------------------------------
# Literal policies and metrics are closed sets: a typo must fail the type check.
if TYPE_CHECKING:
    dm.assign(["ACGT"], ["ACGT"], policy="radios")  # type: ignore[arg-type]
    dm.assign(["ACGT"], ["ACGT"], k="1")  # type: ignore[arg-type]
    dm.assign_hamming(["ACGT"], ["ACGT"], policy="nearest")  # type: ignore[arg-type]
    dm.stream_assign(path, ["ACGT"], metric="edit")  # type: ignore[arg-type]
    dm.assignment_summary([])["assignment_rate"] = "high"  # type: ignore[typeddict-item]
