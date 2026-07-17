"""AssayCode platform identity backed by the DotMatch assignment engine.

AssayCode is the assay-level product surface. The published dotmatch package,
native library, scientific citation, and compatibility contracts remain the
engine of record.
"""

from __future__ import annotations

import dotmatch as engine

__version__ = engine.__version__
PLATFORM_NAME = "AssayCode"
SPEC_NAME = "AssayScript"
ENGINE_NAME = "DotMatch"

__all__ = [
    "ENGINE_NAME",
    "PLATFORM_NAME",
    "SPEC_NAME",
    "__version__",
    "engine",
]
