"""EditWitness: bounded, inspectable CRISPR assay observability analysis.

No network activity, telemetry, alignment, or clinical interpretation.
"""

from ._version import MODEL_VERSION, SCHEMA_VERSION, __version__
from .engine import analyze
from .io import load_manifest
from .models import Manifest

__all__ = ["MODEL_VERSION", "SCHEMA_VERSION", "Manifest", "__version__", "analyze", "load_manifest"]
