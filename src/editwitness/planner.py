"""Minimum-cost or explicitly nonoptimal weighted set cover of counterexamples."""
from __future__ import annotations

from fractions import Fraction
from typing import Literal

from .models import Assay, PanelPlan


def plan_panel(candidates: tuple[Assay, ...], coverage: dict[str, set[str]],
               alternatives: set[str]) -> PanelPlan:
    useful = sorted((a for a in candidates if coverage.get(a.id)), key=lambda a: a.id)
    resolvable: set[str] = set().union(*(coverage[a.id] for a in useful)) if useful else set()
    unresolved = alternatives - resolvable
    if not useful:
        return PanelPlan(
            algorithm="not_needed", optimality="not_applicable", selected_assays=(),
            cost_units=0, resolved_hypotheses=(), unresolved_hypotheses=tuple(sorted(unresolved)),
            note="No candidate distinguishes a remaining alternative." if alternatives else
                 "No equivalent alternative was declared. This is not evidence of assay completeness.",
        )
    names = sorted(resolvable)
    bits = {name: 1 << i for i, name in enumerate(names)}
    masks = [sum(bits[name] for name in coverage[a.id]) for a in useful]
    goal = (1 << len(names)) - 1
    algorithm: Literal["exhaustive_minimum_cost", "greedy_weighted_cover"]
    optimality: Literal["proven_within_declared_candidates", "not_proven"]
    if len(useful) <= 18:
        # Enumerate all subsets of useful assays. This proves minimum total declared cost.
        # Tie-break by number of assays, then lexicographic IDs; never float comparisons.
        n = 1 << len(useful)
        covers = [0] * n
        costs = [0] * n
        best: tuple[int, int, tuple[str, ...]] | None = None
        for subset in range(1, n):
            low = subset & -subset
            index = low.bit_length() - 1
            prior = subset ^ low
            covers[subset] = covers[prior] | masks[index]
            costs[subset] = costs[prior] + useful[index].cost_units
            if covers[subset] == goal:
                ids = tuple(a.id for j, a in enumerate(useful) if subset & (1 << j))
                key = (costs[subset], len(ids), ids)
                if best is None or key < best:
                    best = key
        assert best is not None
        total, _, chosen = best
        algorithm = "exhaustive_minimum_cost"
        optimality = "proven_within_declared_candidates"
    else:
        remaining = set(resolvable)
        selected: list[Assay] = []
        while remaining:
            choices = [a for a in useful if coverage[a.id] & remaining]
            chosen_assay = min(choices, key=lambda a: (
                -Fraction(len(coverage[a.id] & remaining), a.cost_units), a.cost_units, a.id
            ))
            selected.append(chosen_assay)
            remaining -= coverage[chosen_assay.id]
        chosen = tuple(sorted(a.id for a in selected))
        total = sum(a.cost_units for a in selected)
        algorithm = "greedy_weighted_cover"
        optimality = "not_proven"
    return PanelPlan(
        algorithm=algorithm, optimality=optimality, selected_assays=chosen,
        cost_units=total, resolved_hypotheses=tuple(names),
        unresolved_hypotheses=tuple(sorted(unresolved)),
        note="Separates the expected hypothesis from all alternatives separable by these candidates "
             "under this response model. Does not certify real assay performance or resolve unmodeled events.",
    )
