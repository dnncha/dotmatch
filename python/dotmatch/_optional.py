"""Keep optional dataframe stacks off the native matching and CLI import path."""
from __future__ import annotations

from importlib import import_module
from importlib.util import find_spec
from types import ModuleType
from typing import Any


class _DeferredModule:
    """Resolve an optional module only when its API is used.

    Python's import lock handles concurrent initialization. Broken optional
    environments remain visible when requested instead of being hidden.
    """

    def __init__(self, name: str) -> None:
        self._name = name
        self._module: ModuleType | None = None

    def __getattr__(self, name: str) -> Any:
        if self._module is None:
            self._module = import_module(self._name)
        return getattr(self._module, name)


def optional_module(name: str) -> tuple[Any, bool]:
    try:
        available = find_spec(name) is not None
    except (ImportError, ValueError):
        available = False
    return (_DeferredModule(name), True) if available else (None, False)
