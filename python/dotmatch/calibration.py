"""Calibrated quality-aware and joint known-target decoding.

This module is experimental. It is deliberately separate from the deterministic
DotMatch default until public calibration and throughput gates pass.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Mapping, Sequence


DNA = "ACGT"


@dataclass(frozen=True)
class ErrorModel:
    cycle_totals: tuple[int, ...]
    cycle_errors: tuple[int, ...]
    substitution_counts: Mapping[str, Mapping[str, int]]
    prior_strength: float = 100.0

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ProbabilisticCall:
    status: str
    target: str | None
    posterior: float
    second_posterior: float
    likelihood_ratio: float
    candidates: Mapping[str, float]


@dataclass(frozen=True)
class JointCall:
    status: str
    combination: Mapping[str, str] | None
    posterior: float
    second_posterior: float
    candidates: tuple[tuple[Mapping[str, str], float], ...]


@dataclass(frozen=True)
class CalibrationMetrics:
    count: int
    accuracy: float
    mean_confidence: float
    expected_calibration_error: float
    brier_score: float


def fit_error_model(
    observations: Sequence[tuple[str, str, str]],
    *,
    prior_strength: float = 100.0,
) -> ErrorModel:
    """Fit cycle and substitution errors from (observed, expected, qualities).

    Training pairs must already be independently trusted, for example exact or
    uniquely assigned spike-ins. The function does not bootstrap labels from its
    own probabilistic calls.
    """
    if prior_strength <= 0:
        raise ValueError("prior_strength must be positive")
    max_length = max((len(expected) for _observed, expected, _quality in observations), default=0)
    totals = [0] * max_length
    errors = [0] * max_length
    substitutions = {base: {other: 0 for other in DNA if other != base} for base in DNA}
    for observed, expected, quality in observations:
        if len(observed) != len(expected) or len(quality) != len(observed):
            raise ValueError("training observed, expected, and quality strings must have equal length")
        for cycle, (called, truth) in enumerate(zip(observed.upper(), expected.upper())):
            if called not in DNA or truth not in DNA:
                continue
            totals[cycle] += 1
            if called != truth:
                errors[cycle] += 1
                substitutions[truth][called] += 1
    return ErrorModel(
        cycle_totals=tuple(totals),
        cycle_errors=tuple(errors),
        substitution_counts=substitutions,
        prior_strength=float(prior_strength),
    )


def decode(
    observed: str,
    quality: str,
    targets: Sequence[str],
    model: ErrorModel,
    *,
    priors: Mapping[str, float] | None = None,
    posterior_min: float = 0.99,
    likelihood_ratio_min: float = 10.0,
) -> ProbabilisticCall:
    if len(observed) != len(quality):
        raise ValueError("observed and quality strings must have equal length")
    if not 0.0 < posterior_min <= 1.0:
        raise ValueError("posterior_min must be in (0, 1]")
    if likelihood_ratio_min < 1.0:
        raise ValueError("likelihood_ratio_min must be at least 1")
    unique_targets = list(dict.fromkeys(target.upper() for target in targets))
    if not unique_targets:
        raise ValueError("at least one target is required")

    log_scores: dict[str, float] = {}
    default_prior = 1.0 / len(unique_targets)
    for target in unique_targets:
        if len(target) != len(observed):
            continue
        prior = default_prior if priors is None else float(priors.get(target, 0.0))
        if prior <= 0:
            continue
        score = math.log(prior)
        for cycle, (called, truth, qchar) in enumerate(zip(observed.upper(), target, quality)):
            score += math.log(_emission_probability(called, truth, ord(qchar) - 33, cycle, model))
        log_scores[target] = score
    if not log_scores:
        return ProbabilisticCall("none", None, 0.0, 0.0, 0.0, {})

    probabilities = _normalize_logs(log_scores)
    ranked = sorted(probabilities.items(), key=lambda item: (-item[1], item[0]))
    best_target, best = ranked[0]
    second = ranked[1][1] if len(ranked) > 1 else 0.0
    ratio = math.inf if second == 0.0 else best / second
    status = "unique" if best >= posterior_min and ratio >= likelihood_ratio_min else "ambiguous"
    return ProbabilisticCall(
        status=status,
        target=best_target if status == "unique" else None,
        posterior=best,
        second_posterior=second,
        likelihood_ratio=ratio,
        candidates=dict(ranked),
    )


def decode_joint(
    components: Mapping[str, Mapping[str, float]],
    allowed_combinations: Sequence[Mapping[str, str]],
    *,
    priors: Sequence[float] | None = None,
    posterior_min: float = 0.99,
    likelihood_ratio_min: float = 10.0,
) -> JointCall:
    """Decode an allowed tuple from independently calibrated component masses."""
    if not components:
        raise ValueError("at least one component posterior table is required")
    if not allowed_combinations:
        raise ValueError("at least one allowed combination is required")
    if priors is not None and len(priors) != len(allowed_combinations):
        raise ValueError("joint priors must match allowed combinations")

    scores: list[tuple[Mapping[str, str], float]] = []
    default_prior = 1.0 / len(allowed_combinations)
    for index, combination in enumerate(allowed_combinations):
        if set(combination) != set(components):
            raise ValueError("each combination must define exactly the component names")
        prior = default_prior if priors is None else float(priors[index])
        if prior <= 0:
            continue
        log_score = math.log(prior)
        possible = True
        for component, target in combination.items():
            mass = float(components[component].get(target, 0.0))
            if mass <= 0:
                possible = False
                break
            log_score += math.log(mass)
        if possible:
            scores.append((dict(combination), log_score))
    if not scores:
        return JointCall("none", None, 0.0, 0.0, ())

    maximum = max(score for _combination, score in scores)
    weights = [(combination, math.exp(score - maximum)) for combination, score in scores]
    total = sum(weight for _combination, weight in weights)
    ranked = sorted(
        ((combination, weight / total) for combination, weight in weights),
        key=lambda item: (-item[1], tuple(sorted(item[0].items()))),
    )
    best_combination, best = ranked[0]
    second = ranked[1][1] if len(ranked) > 1 else 0.0
    ratio = math.inf if second == 0.0 else best / second
    status = "unique" if best >= posterior_min and ratio >= likelihood_ratio_min else "ambiguous"
    return JointCall(
        status=status,
        combination=best_combination if status == "unique" else None,
        posterior=best,
        second_posterior=second,
        candidates=tuple(ranked),
    )


def calibration_metrics(
    calls: Sequence[tuple[float, bool]],
    *,
    bins: int = 10,
) -> CalibrationMetrics:
    if bins <= 0:
        raise ValueError("bins must be positive")
    if not calls:
        return CalibrationMetrics(0, 0.0, 0.0, 0.0, 0.0)
    for confidence, _correct in calls:
        if not 0.0 <= confidence <= 1.0:
            raise ValueError("confidence values must be between 0 and 1")
    count = len(calls)
    accuracy = sum(correct for _confidence, correct in calls) / count
    mean_confidence = sum(confidence for confidence, _correct in calls) / count
    brier = sum((confidence - float(correct)) ** 2 for confidence, correct in calls) / count
    ece = 0.0
    for bin_index in range(bins):
        lower = bin_index / bins
        upper = (bin_index + 1) / bins
        bucket = [
            (confidence, correct)
            for confidence, correct in calls
            if lower <= confidence < upper or (bin_index == bins - 1 and confidence == 1.0)
        ]
        if not bucket:
            continue
        bucket_confidence = sum(confidence for confidence, _correct in bucket) / len(bucket)
        bucket_accuracy = sum(correct for _confidence, correct in bucket) / len(bucket)
        ece += len(bucket) / count * abs(bucket_confidence - bucket_accuracy)
    return CalibrationMetrics(count, accuracy, mean_confidence, ece, brier)


def threshold_for_fdr(
    calls: Sequence[tuple[float, bool]],
    *,
    max_fdr: float,
) -> float | None:
    """Return the lowest confidence cutoff whose accepted prefix meets max FDR.

    The selected prefix maximizes accepted calls. Ties are deterministic.
    """
    if not 0.0 <= max_fdr < 1.0:
        raise ValueError("max_fdr must be in [0, 1)")
    ranked = sorted(calls, key=lambda item: (-item[0], not item[1]))
    errors = 0
    accepted_cutoff: float | None = None
    for index, (confidence, correct) in enumerate(ranked, start=1):
        if not 0.0 <= confidence <= 1.0:
            raise ValueError("confidence values must be between 0 and 1")
        errors += int(not correct)
        if errors / index <= max_fdr:
            accepted_cutoff = confidence
    return accepted_cutoff


def smoothed_abundance_priors(
    counts: Mapping[str, int],
    *,
    alpha: float = 1.0,
) -> dict[str, float]:
    if alpha <= 0:
        raise ValueError("alpha must be positive")
    if not counts:
        raise ValueError("counts must not be empty")
    if any(isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in counts.values()):
        raise ValueError("counts must be non-negative integers")
    total = sum(counts.values()) + alpha * len(counts)
    return {target: (count + alpha) / total for target, count in counts.items()}


def _emission_probability(
    called: str,
    truth: str,
    quality: int,
    cycle: int,
    model: ErrorModel,
) -> float:
    if called not in DNA or truth not in DNA:
        return 0.25
    quality = min(60, max(0, quality))
    phred_error = 10.0 ** (-quality / 10.0)
    if cycle < len(model.cycle_totals):
        total = model.cycle_totals[cycle]
        errors = model.cycle_errors[cycle]
    else:
        total = errors = 0
    empirical_error = (errors + 0.5) / (total + 1.0)
    weight = total / (total + model.prior_strength)
    error_rate = min(0.75, max(1e-9, weight * empirical_error + (1.0 - weight) * phred_error))
    if called == truth:
        return max(1e-12, 1.0 - error_rate)

    substitution_row = model.substitution_counts.get(truth, {})
    substitution_total = sum(substitution_row.values())
    substitution_probability = (substitution_row.get(called, 0) + 1.0) / (
        substitution_total + 3.0
    )
    return max(1e-12, error_rate * substitution_probability)


def _normalize_logs(scores: Mapping[str, float]) -> dict[str, float]:
    maximum = max(scores.values())
    weights = {target: math.exp(score - maximum) for target, score in scores.items()}
    total = sum(weights.values())
    return {target: weight / total for target, weight in weights.items()}
