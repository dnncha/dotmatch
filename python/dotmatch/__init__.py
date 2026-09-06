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
    FastqRecord,
    Matcher,
    MatchResult,
    PosteriorAssignment,
    StreamAssignment,
    alphabet_policy,
    assign,
    assign_dataframe,
    assign_exact,
    assign_hamming,
    assignment_summary,
    assign_posterior,
    assignments_to_anndata,
    counts_tsv_to_anndata,
    distance,
    distance_leq,
    iter_fastq,
    load_targets,
    results_to_dataframe,
    status_name,
    stream_assign,
    targets_from_dataframe,
    write_assignments_tsv,
)
from .feature_matrix import FeatureMatrixResult, build_feature_matrix

# Advanced / optional integrations (import submodules to avoid heavy dep cost)
# from . import anndata as anndata  # if you have the extra
# from . import multiqc as multiqc

from . import tl as tl  # scverse-style tools (dotmatch.tl.assign_features, dotmatch.tl.feature_counts, ...)


def _source_tree_version() -> Optional[str]:
    here = Path(__file__).resolve()
    # A --target/vendor install can sit beneath somebody else's pyproject.
    # Only our actual source layout is allowed to override installed metadata.
    if here.parent.name != "dotmatch" or here.parents[1].name != "python":
        return None
    pyproject = here.parents[2] / "pyproject.toml"
    if not pyproject.exists():
        return None
    project = re.search(r"(?ms)^\[project\]\s*$(.*?)(?=^\[|\Z)", pyproject.read_text(encoding="utf-8"))
    if not project or not re.search(r'^name\s*=\s*"dotmatch"\s*$', project.group(1), re.MULTILINE):
        return None
    match = re.search(r'^version\s*=\s*"([^"]+)"', project.group(1), flags=re.MULTILINE)
    return match.group(1) if match else None


try:
    __version__ = _source_tree_version() or _metadata_version("dotmatch")
except PackageNotFoundError:
    __version__ = "0.5.0"

__all__ = [
    "__version__",
    "MATCH_AMBIGUOUS",
    "MATCH_INVALID",
    "MATCH_NONE",
    "MATCH_UNIQUE",
    "AssignmentStats",
    "FastqRecord",
    "Matcher",
    "MatchResult",
    "PosteriorAssignment",
    "StreamAssignment",
    "alphabet_policy",
    "assign",
    "assign_dataframe",
    "assign_exact",
    "assign_hamming",
    "assignment_summary",
    "assign_posterior",
    "assignments_to_anndata",
    "counts_tsv_to_anndata",
    "distance",
    "distance_leq",
    "iter_fastq",
    "load_targets",
    "results_to_dataframe",
    "status_name",
    "stream_assign",
    "targets_from_dataframe",
    "write_assignments_tsv",
    "FeatureMatrixResult",
    "build_feature_matrix",
    "tl",
]
