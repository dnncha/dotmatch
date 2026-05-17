import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "check_barcode_public_surface.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_barcode_public_surface", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_barcode_public_surface_is_easy_and_trustworthy():
    checker = _load_checker()

    result = checker.audit(ROOT)

    assert result.ok, result.failures
    assert any("one-command barcode autopsy" in item for item in result.passed)


def test_barcode_public_surface_rejects_hype(tmp_path):
    checker = _load_checker()
    (tmp_path / "docs").mkdir()
    (tmp_path / "docs" / "barcode-science-readiness.md").write_text(
        "DotMatch will dominate and replace Cutadapt.\n",
        encoding="utf-8",
    )
    (tmp_path / "README.md").write_text(
        "DotMatch will dominate and replace Cutadapt.\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert not result.ok
    assert any("overbroad or hype wording" in failure for failure in result.failures)
