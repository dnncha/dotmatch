"""Run the executable examples embedded in public modules."""

from __future__ import annotations

import doctest
import types

import dotmatch._types
import dotmatch.core
import dotmatch.panel
import pytest

MODULES = [
    dotmatch._types,
    dotmatch.core,
    dotmatch.panel,
]


@pytest.mark.parametrize("module", MODULES, ids=[module.__name__ for module in MODULES])
def test_module_doctests(module: types.ModuleType) -> None:
    failures, attempted = doctest.testmod(module, verbose=False, raise_on_error=False)
    assert attempted > 0, f"{module!r} has no doctests"
    assert failures == 0, f"{failures} doctest example(s) failed in {module!r}"
