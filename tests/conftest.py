from pathlib import Path

import pytest

from editwitness import load_manifest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def demo():
    return load_manifest(ROOT / "examples/demo.json")


@pytest.fixture
def paired():
    return load_manifest(ROOT / "examples/paired_end.json")
