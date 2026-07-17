"""Deterministic design-time simulation for known-target panels.

This module is experimental. It estimates assignment behavior under an explicit
substitution-error model; it is not a substitute for held-out empirical
validation or a promise about a sequencing platform.
"""

from __future__ import annotations

import random
from dataclasses import asdict, dataclass
from typing import Mapping, Sequence

DNA = "ACGT"


@dataclass(frozen=True)
class SimulationResult:
    seed: int
    total_reads: int
    correct_unique: int
    misassigned_unique: int
    ambiguous: int
    none: int
    usable_yield: float
    ambiguity_rate: float
    no_call_rate: float
    false_discovery_rate: float
    confusion: Mapping[str, Mapping[str, int]]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def simulate_panel(
    targets: Mapping[str, str],
    *,
    reads_per_target: int = 1000,
    k: int = 1,
    error_rate: float | Sequence[float] = 0.01,
    seed: int = 1,
) -> SimulationResult:
    """Simulate substitutions and assign each read by bounded Hamming distance.

    Every target receives the same number of reads. A call is unique only when
    exactly one target is within the configured distance. Duplicate target
    sequences are rejected because their truth labels are not identifiable.
    """
    if not targets:
        raise ValueError("targets must not be empty")
    if isinstance(reads_per_target, bool) or not isinstance(reads_per_target, int) or reads_per_target <= 0:
        raise ValueError("reads_per_target must be a positive integer")
    if isinstance(k, bool) or not isinstance(k, int) or k < 0:
        raise ValueError("k must be a non-negative integer")

    normalized = {str(name): str(sequence).upper() for name, sequence in targets.items()}
    if len(normalized) != len(targets) or any(not name for name in normalized):
        raise ValueError("target names must be unique and non-empty")
    lengths = {len(sequence) for sequence in normalized.values()}
    if len(lengths) != 1 or not lengths or 0 in lengths:
        raise ValueError("all target sequences must have the same positive length")
    if any(set(sequence) - set(DNA) for sequence in normalized.values()):
        raise ValueError("target sequences must contain only A, C, G, and T")
    if len(set(normalized.values())) != len(normalized):
        raise ValueError("duplicate target sequences are not identifiable")

    length = next(iter(lengths))
    rates = _cycle_rates(error_rate, length)
    rng = random.Random(seed)
    counts = {"correct": 0, "misassigned": 0, "ambiguous": 0, "none": 0}
    confusion: dict[str, dict[str, int]] = {name: {} for name in normalized}

    for truth_name, truth_sequence in normalized.items():
        for _ in range(reads_per_target):
            observed = _mutate(truth_sequence, rates, rng)
            candidates = [
                name
                for name, target_sequence in normalized.items()
                if _hamming(observed, target_sequence) <= k
            ]
            if not candidates:
                counts["none"] += 1
                continue
            if len(candidates) > 1:
                counts["ambiguous"] += 1
                continue
            called = candidates[0]
            confusion[truth_name][called] = confusion[truth_name].get(called, 0) + 1
            if called == truth_name:
                counts["correct"] += 1
            else:
                counts["misassigned"] += 1

    total = reads_per_target * len(normalized)
    unique = counts["correct"] + counts["misassigned"]
    return SimulationResult(
        seed=seed,
        total_reads=total,
        correct_unique=counts["correct"],
        misassigned_unique=counts["misassigned"],
        ambiguous=counts["ambiguous"],
        none=counts["none"],
        usable_yield=unique / total,
        ambiguity_rate=counts["ambiguous"] / total,
        no_call_rate=counts["none"] / total,
        false_discovery_rate=counts["misassigned"] / unique if unique else 0.0,
        confusion=confusion,
    )


def _cycle_rates(error_rate: float | Sequence[float], length: int) -> tuple[float, ...]:
    if isinstance(error_rate, (int, float)) and not isinstance(error_rate, bool):
        rates = (float(error_rate),) * length
    else:
        rates = tuple(float(value) for value in error_rate)
        if len(rates) != length:
            raise ValueError("per-cycle error rates must match target length")
    if any(rate < 0.0 or rate > 1.0 for rate in rates):
        raise ValueError("error rates must be between 0 and 1")
    return rates


def _mutate(sequence: str, rates: Sequence[float], rng: random.Random) -> str:
    called: list[str] = []
    for base, rate in zip(sequence, rates):
        if rng.random() < rate:
            alternatives = [candidate for candidate in DNA if candidate != base]
            called.append(alternatives[rng.randrange(3)])
        else:
            called.append(base)
    return "".join(called)


def _hamming(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right))
