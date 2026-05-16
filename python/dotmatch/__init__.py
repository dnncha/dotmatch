from importlib.metadata import PackageNotFoundError, version as _metadata_version

from .core import (
    MATCH_AMBIGUOUS,
    MATCH_INVALID,
    MATCH_NONE,
    MATCH_UNIQUE,
    AssignmentStats,
    Matcher,
    MatchResult,
    alphabet_policy,
    assign,
    distance,
    distance_leq,
)

try:
    __version__ = _metadata_version("dotmatch")
except PackageNotFoundError:
    __version__ = "0.1.1"

__all__ = [
    "__version__",
    "MATCH_AMBIGUOUS",
    "MATCH_INVALID",
    "MATCH_NONE",
    "MATCH_UNIQUE",
    "AssignmentStats",
    "Matcher",
    "MatchResult",
    "alphabet_policy",
    "assign",
    "distance",
    "distance_leq",
]
