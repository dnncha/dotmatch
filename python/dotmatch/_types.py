"""Strictly typed vocabulary shared by DotMatch's public API.

The native library speaks in signed integers. This module gives those integers a
closed ``IntEnum``, gives every stringly-typed policy a ``Literal``, and gives
every tabular result a ``TypedDict``. Parsing helpers live here too so that
untrusted or loosely typed input is narrowed once, at the boundary, instead of
being cast at each use site.
"""

from __future__ import annotations

from enum import IntEnum
from typing import Final, Literal, Mapping, NamedTuple, Optional, TypedDict, Union


class MatchStatus(IntEnum):
    """Stable assignment status returned by the native kernels."""

    INVALID = -1
    NONE = 0
    UNIQUE = 1
    AMBIGUOUS = 2

    def __str__(self) -> str:
        """Render like the integer it replaces, on every supported Python version."""
        return str(int(self))


MATCH_INVALID: Final[MatchStatus] = MatchStatus.INVALID
MATCH_NONE: Final[MatchStatus] = MatchStatus.NONE
MATCH_UNIQUE: Final[MatchStatus] = MatchStatus.UNIQUE
MATCH_AMBIGUOUS: Final[MatchStatus] = MatchStatus.AMBIGUOUS

StatusName = Literal["invalid", "none", "unique", "ambiguous"]
STATUS_NAMES: Final[Mapping[MatchStatus, StatusName]] = {
    MATCH_INVALID: "invalid",
    MATCH_NONE: "none",
    MATCH_UNIQUE: "unique",
    MATCH_AMBIGUOUS: "ambiguous",
}

Metric = Literal["levenshtein", "hamming", "exact"]
AssignmentPolicy = Literal["radius", "best"]
DistanceMetric = Literal["hamming", "levenshtein", "sequence-levenshtein"]
PanelMetric = DistanceMetric
PanelStatus = Literal["pass", "warn", "fail"]
AssayMode = Literal["count", "demux", "pair-count"]
AssayType = Literal[
    "crispr",
    "feature_barcode",
    "inline_barcode",
    "amplicon_panel",
    "oligo_adapter",
    "generic",
]
AssayStatus = Literal["ready", "draft"]
ReadInput = Union[str, bytes]


class TargetRow(NamedTuple):
    """One known target: the identifier written to output plus its DNA sequence."""

    target_id: str
    sequence: str


class AssignmentCounts(TypedDict):
    """Integer counters accumulated while streaming assignments."""

    total_reads: int
    assigned_unique: int
    assigned_exact: int
    assigned_corrected: int
    ambiguous: int
    unmatched: int
    invalid: int


class AssignmentSummary(AssignmentCounts):
    """Counters plus the four rates derived from ``total_reads``."""

    assignment_rate: float
    ambiguous_rate: float
    no_match_rate: float
    invalid_rate: float


class PosteriorSummary(TypedDict):
    """Serializable view of one posterior assignment call."""

    target_index: int
    posterior: float
    second_posterior: float
    status: int
    posteriors: list[float]


class PanelCertificate(TypedDict):
    """The safety certificate written by ``dotmatch panel check``."""

    schema_version: int
    status: PanelStatus
    panel_grade: str
    n_barcodes: int
    lengths: list[int]
    assignment_metric: PanelMetric
    configured_assignment_k: int
    exact_error_sphere_radius: int
    maximum_exact_certificate_radius: int
    ambiguity_policy: AssignmentPolicy
    safe_for_k0: bool
    safe_for_k1_hamming: bool
    safe_for_k1_levenshtein: bool
    safe_for_k2_hamming: bool
    safe_for_k2_levenshtein: bool
    minimum_target_safe_hamming_radius: int
    maximum_target_safe_hamming_radius: int
    minimum_hamming_distance: Optional[int]
    minimum_levenshtein_distance: Optional[int]
    minimum_sequence_levenshtein_distance: Optional[int]
    guaranteed_substitution_correction: int
    guaranteed_substitution_detection: int
    ambiguous_error_spheres: int
    unsafe_k1_variants: int
    unsafe_configured_radius_variants: int
    silent_assignment_risk: int
    duplicate_pairs: int
    collision_pairs: int
    prefix_collisions: int
    suffix_collisions: int
    reverse_complement_warnings: int
    cycle_balance_warnings: int
    cycle_balance_score: float
    forbidden_motif_hits: int
    maximum_homopolymer: int
    certified_dotmatch_command: str
    scope: str


class PanelTraceRow(TypedDict):
    """One design-search step."""

    restart: Union[int, str]
    step: int
    candidate_id: str
    sequence: str
    score: str
    selected_count: int
    reason: str


def require_str(value: object, *, field: str) -> str:
    """Narrow a boundary value to ``str`` or fail with its field name."""
    if isinstance(value, str):
        return value
    raise TypeError(f"{field} must be a string, got {type(value).__name__}")


def require_int(value: object, *, field: str) -> int:
    """Narrow a boundary value to an integer without silent truncation.

    Examples
    --------
    >>> require_int("42", field="k")
    42
    >>> require_int(True, field="k")
    Traceback (most recent call last):
        ...
    TypeError: k must be an integer, got bool
    """
    if isinstance(value, bool):
        raise TypeError(f"{field} must be an integer, got bool")
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, str):
        try:
            return int(value)
        except ValueError as exc:
            raise TypeError(f"{field} must be an integer, got {value!r}") from exc
    raise TypeError(f"{field} must be an integer, got {type(value).__name__}")


def require_float(value: object, *, field: str) -> float:
    """Narrow a boundary value to ``float`` without accepting booleans.

    Examples
    --------
    >>> require_float("0.5", field="rate")
    0.5
    """
    if isinstance(value, bool):
        raise TypeError(f"{field} must be a number, got bool")
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError as exc:
            raise TypeError(f"{field} must be a number, got {value!r}") from exc
    raise TypeError(f"{field} must be a number, got {type(value).__name__}")

