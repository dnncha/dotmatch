from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_agent_discovery.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_agent_discovery", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_agent_discovery_manifest_schema_copies_and_surfaces_are_valid():
    checker = _load_checker()

    assert checker.validate_schema(ROOT) == []
    assert checker.validate_manifest(ROOT) == []
    assert checker.validate_copies(ROOT) == []
    assert checker.validate_surfaces(ROOT) == []


def test_agent_discovery_measurement_reaches_full_score():
    checker = _load_checker()

    report = checker.local_measure(ROOT)

    assert report["score"] == report["maximum"] == 12
    assert all(item["passed"] for item in report["checks"])
