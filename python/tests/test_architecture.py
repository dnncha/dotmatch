"""Mutation-style checks: architectural regressions must make the gate fail."""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "architecture", ROOT / "scripts/check_architecture.py"
)
assert SPEC is not None and SPEC.loader is not None
architecture = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(architecture)


def write_module(root, name, source):
    path = root / "python" / "dotmatch" / f"{name}.py"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(source, encoding="utf-8")


def test_repository_boundaries():
    assert architecture.check(ROOT) == []


@pytest.mark.parametrize(
    "source",
    [
        "from .core import Matcher\n",
        "from dotmatch import core\n",
        "def parse():\n    from . import core\n",
    ],
)
def test_leaf_cannot_import_matching_even_lazily(tmp_path, source):
    write_module(tmp_path, "target_io", source)
    assert any("low-level module" in error for error in architecture.check(tmp_path))


@pytest.mark.parametrize(
    "source,reason",
    [
        ("import quickdna\n", "compatibility"),
        ("from scripts import report\n", "repository tooling"),
        ("from .cli import main\n", "CLI orchestration"),
        ("try:\n    import pandas\nexcept ImportError:\n    pass\n", "imported lazily"),
        ("import ctypes as ct\nct.CDLL('untrusted.so')\n", "_bindings.py"),
        ("from ctypes import CDLL as load\nload('untrusted.so')\n", "_bindings.py"),
    ],
)
def test_forbidden_dependencies(tmp_path, source, reason):
    write_module(tmp_path, "core", source)
    assert any(reason in error for error in architecture.check(tmp_path))


def test_eager_cycle_reports_path(tmp_path):
    write_module(tmp_path, "alpha", "from .beta import value\n")
    write_module(tmp_path, "beta", "from .alpha import value\n")
    assert any(
        "dotmatch.alpha -> dotmatch.beta -> dotmatch.alpha" in error
        for error in architecture.check(tmp_path)
    )


def test_lazy_integration_and_allowed_leaf_dependency(tmp_path):
    write_module(tmp_path, "dataframes", "def convert():\n    import pandas\n")
    write_module(tmp_path, "count_io", "from ._validation import identifier\n")
    write_module(tmp_path, "_validation", "import operator\n")
    assert architecture.check(tmp_path) == []


def test_missing_package_fails(tmp_path):
    assert architecture.check(tmp_path) == ["python/dotmatch: no runtime modules found"]
