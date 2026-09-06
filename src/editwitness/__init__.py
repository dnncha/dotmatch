"""EditWitness: bounded, inspectable CRISPR assay observability analysis.

No network activity, telemetry, alignment, or clinical interpretation.
"""

from ._version import EXACT_MODEL_VERSION, MODEL_VERSION, SCHEMA_VERSION, __version__
from .compare import compare_models
from .generate import expand_deletions
from .engine import analyze
from .io import load_manifest
from .models import Manifest

__all__ = ["EXACT_MODEL_VERSION", "MODEL_VERSION", "SCHEMA_VERSION", "compare_models", "expand_deletions", "Manifest", "__version__", "analyze", "load_manifest"]
