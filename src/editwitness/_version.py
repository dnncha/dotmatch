"""Package, schema and scientific-model identifiers."""

from typing import Final

__version__ = "0.2.0a2"
SCHEMA_VERSION = "1.1"
LEGACY_MODEL_VERSION: Final = "original-sites-presence-v1"
EXACT_MODEL_VERSION: Final = "exact-local-sequence-presence-v2"
# Keep the original constant for downstream users of the original-site primitive.
MODEL_VERSION = LEGACY_MODEL_VERSION
SUPPORTED_MODELS = (LEGACY_MODEL_VERSION, EXACT_MODEL_VERSION)
