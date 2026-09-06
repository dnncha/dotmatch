"""EditWitness: bounded, inspectable CRISPR assay observability analysis.

No network activity, telemetry, alignment, or clinical interpretation.
"""

from ._version import (
    EXACT_MODEL_VERSION,
    MODEL_VERSION,
    SCHEMA_VERSION,
    SUPPORTED_MODELS,
    __version__,
)
from .design import expand_deletions
from .engine import analyze
from .io import load_manifest
from .models import Manifest

__all__ = [
    "EXACT_MODEL_VERSION",
    "MODEL_VERSION",
    "SCHEMA_VERSION",
    "SUPPORTED_MODELS",
    "Manifest",
    "__version__",
    "analyze",
    "expand_deletions",
    "load_manifest",
]
