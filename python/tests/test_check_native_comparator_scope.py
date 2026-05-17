import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_native_comparator_scope.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_native_comparator_scope", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_native_comparator_scope_accepts_present_document(tmp_path):
    checker = _load_checker()
    path = tmp_path / "docs" / "native-comparator-scope.md"
    path.parent.mkdir(parents=True)
    path.write_text("# Native Comparator Scope\n", encoding="utf-8")

    assert checker.check(tmp_path) == []


def test_native_comparator_scope_reports_missing_document(tmp_path):
    checker = _load_checker()

    failures = checker.check(tmp_path)

    assert any("missing native comparator scope document" in failure for failure in failures)


def test_native_comparator_scope_reports_empty_document(tmp_path):
    checker = _load_checker()
    path = tmp_path / "docs" / "native-comparator-scope.md"
    path.parent.mkdir(parents=True)
    path.write_text("", encoding="utf-8")

    failures = checker.check(tmp_path)

    assert any("empty native comparator scope document" in failure for failure in failures)
