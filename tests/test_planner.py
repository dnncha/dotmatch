import itertools
import random

from editwitness.models import Assay, Interval
from editwitness.planner import plan_panel


def assay(name, cost=1):
    return Assay(id=name, left_primer=Interval(start=0, end=2), right_primer=Interval(start=8, end=10), cost_units=cost)


def test_minimum_cost_not_maximum_single_assay_coverage():
    candidates = (assay("all", 10), assay("left", 2), assay("right", 2))
    result = plan_panel(candidates, {"all": {"x", "y"}, "left": {"x"}, "right": {"y"}}, {"x", "y", "z"})
    assert result.selected_assays == ("left", "right")
    assert result.cost_units == 4
    assert result.unresolved_hypotheses == ("z",)
    assert result.optimality == "proven_within_declared_candidates"


def test_exact_planner_against_independent_combinations_oracle():
    rng = random.Random(444)
    for _ in range(150):
        candidates = tuple(assay(f"a{i}", rng.randint(1, 15)) for i in range(rng.randint(1, 7)))
        coverage = {a.id: {str(j) for j in range(7) if rng.random() < .4} for a in candidates}
        goal = set().union(*coverage.values())
        expected = None
        for size in range(len(candidates) + 1):
            for subset in itertools.combinations(candidates, size):
                covered = set().union(*(coverage[a.id] for a in subset)) if subset else set()
                if covered == goal:
                    key = (sum(a.cost_units for a in subset), len(subset), tuple(sorted(a.id for a in subset)))
                    if expected is None or key < expected:
                        expected = key
        actual = plan_panel(candidates, coverage, {str(j) for j in range(7)})
        assert (actual.cost_units, len(actual.selected_assays), actual.selected_assays) == expected


def test_greedy_fallback_does_not_claim_optimality():
    candidates = tuple(assay(f"a{i:02d}") for i in range(19))
    coverage = {a.id: {f"h{i}"} for i, a in enumerate(candidates)}
    result = plan_panel(candidates, coverage, set().union(*coverage.values()))
    assert result.algorithm == "greedy_weighted_cover"
    assert result.optimality == "not_proven"
    assert len(result.selected_assays) == 19


def test_duplicate_coverage_lexicographic_tie_and_no_candidates():
    result = plan_panel((assay("z"), assay("a")), {"z": {"h"}, "a": {"h"}}, {"h"})
    assert result.selected_assays == ("a",)
    empty = plan_panel((), {}, {"h"})
    assert empty.selected_assays == ()
    assert empty.unresolved_hypotheses == ("h",)
