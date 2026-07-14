from __future__ import annotations

import math

import pytest

from dotmatch.calibration import (
    calibration_metrics,
    decode,
    decode_joint,
    fit_error_model,
    smoothed_abundance_priors,
    threshold_for_fdr,
)


def test_fit_error_model_records_cycle_and_substitution_errors() -> None:
    model = fit_error_model(
        [
            ("ACGT", "ACGT", "IIII"),
            ("AGGT", "ACGT", "IIII"),
        ],
        prior_strength=10.0,
    )

    assert model.cycle_totals == (2, 2, 2, 2)
    assert model.cycle_errors == (0, 1, 0, 0)
    assert model.substitution_counts["C"]["G"] == 1


def test_decode_returns_selective_unique_call() -> None:
    model = fit_error_model([("ACGT", "ACGT", "IIII")] * 20)

    call = decode(
        "ACGT",
        "IIII",
        ["ACGT", "AGGT", "TTTT"],
        model,
        posterior_min=0.99,
        likelihood_ratio_min=10.0,
    )

    assert call.status == "unique"
    assert call.target == "ACGT"
    assert call.posterior > 0.999
    assert call.likelihood_ratio > 10.0


def test_decode_abstains_when_targets_are_not_separated() -> None:
    model = fit_error_model([])
    call = decode(
        "ANGT",
        "!!!!",
        ["ACGT", "AGGT"],
        model,
        posterior_min=0.9,
        likelihood_ratio_min=10.0,
    )

    assert call.status == "ambiguous"
    assert call.target is None
    assert call.posterior == pytest.approx(0.5)


def test_decode_uses_smoothed_abundance_priors_without_zeroing_targets() -> None:
    priors = smoothed_abundance_priors({"ACGT": 99, "AGGT": 0}, alpha=1.0)

    assert sum(priors.values()) == pytest.approx(1.0)
    assert priors["AGGT"] > 0.0


def test_joint_decode_uses_allowed_combinations_to_resolve_evidence() -> None:
    call = decode_joint(
        {
            "sample": {"s1": 0.9, "s2": 0.1},
            "guide": {"g1": 0.45, "g2": 0.55},
        },
        [
            {"sample": "s1", "guide": "g2"},
            {"sample": "s2", "guide": "g1"},
        ],
        posterior_min=0.8,
        likelihood_ratio_min=5.0,
    )

    assert call.status == "unique"
    assert call.combination == {"sample": "s1", "guide": "g2"}
    assert call.posterior > 0.8


def test_joint_decode_abstains_on_symmetric_evidence() -> None:
    call = decode_joint(
        {
            "left": {"a": 0.5, "b": 0.5},
            "right": {"x": 0.5, "y": 0.5},
        },
        [
            {"left": "a", "right": "x"},
            {"left": "b", "right": "y"},
        ],
        posterior_min=0.8,
    )

    assert call.status == "ambiguous"
    assert call.combination is None
    assert call.posterior == pytest.approx(0.5)


def test_calibration_metrics_report_ece_and_brier() -> None:
    metrics = calibration_metrics(
        [(0.9, True), (0.8, True), (0.2, False), (0.1, False)],
        bins=2,
    )

    assert metrics.count == 4
    assert metrics.accuracy == 0.5
    assert metrics.mean_confidence == pytest.approx(0.5)
    assert metrics.expected_calibration_error == pytest.approx(0.15)
    assert metrics.brier_score == pytest.approx(0.025)


def test_threshold_for_fdr_maximizes_valid_prefix() -> None:
    cutoff = threshold_for_fdr(
        [(0.99, True), (0.95, True), (0.90, False), (0.80, True)],
        max_fdr=0.25,
    )

    assert cutoff == 0.8


def test_threshold_for_fdr_returns_none_when_no_call_is_safe() -> None:
    assert threshold_for_fdr([(0.9, False)], max_fdr=0.0) is None


@pytest.mark.parametrize(
    "observed, expected, quality",
    [
        ("AC", "A", "II"),
        ("AC", "AC", "I"),
    ],
)
def test_fit_error_model_rejects_misaligned_training_rows(
    observed: str,
    expected: str,
    quality: str,
) -> None:
    with pytest.raises(ValueError, match="equal length"):
        fit_error_model([(observed, expected, quality)])


def test_joint_decode_rejects_incomplete_constraints() -> None:
    with pytest.raises(ValueError, match="exactly the component names"):
        decode_joint(
            {"sample": {"s1": 1.0}, "guide": {"g1": 1.0}},
            [{"sample": "s1"}],
        )
