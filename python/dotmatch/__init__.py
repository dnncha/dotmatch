import re
from importlib.metadata import PackageNotFoundError, version as _metadata_version
from pathlib import Path
from typing import Optional

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

def _source_tree_version() -> Optional[str]:
    pyproject = Path(__file__).resolve().parents[2] / "pyproject.toml"
    if not pyproject.exists():
        return None
    match = re.search(r'^version\s*=\s*"([^"]+)"', pyproject.read_text(encoding="utf-8"), flags=re.MULTILINE)
    return match.group(1) if match else None


try:
    __version__ = _source_tree_version() or _metadata_version("dotmatch")
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
