from __future__ import annotations

import importlib.util
import tarfile
import sys
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_python_wheel.py"


VALID_METADATA = """\
Metadata-Version: 2.4
Name: dotmatch
Version: 0.1.0
Summary: Deterministic known-target short-DNA assignment for CRISPR guide counting, barcode demultiplexing, and FASTQ workflows
License-Expression: Apache-2.0
Keywords: bioinformatics,computational biology,CRISPR,FASTQ,known-target assignment,barcode demultiplexing,edit distance
Classifier: Intended Audience :: Science/Research
Classifier: Topic :: Scientific/Engineering :: Bio-Informatics
Project-URL: Homepage, https://github.com/dnncha/dotmatch
Project-URL: Repository, https://github.com/dnncha/dotmatch
Project-URL: Issues, https://github.com/dnncha/dotmatch/issues
Project-URL: Documentation, https://github.com/dnncha/dotmatch#readme

Long description.
"""


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_python_wheel", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_sdist(path: Path, metadata: str = VALID_METADATA) -> Path:
    sdist = path / "dotmatch-0.1.0.tar.gz"
    metadata_path = path / "dotmatch-0.1.0" / "PKG-INFO"
    metadata_path.parent.mkdir(parents=True)
    metadata_path.write_text(metadata, encoding="utf-8")
    with tarfile.open(sdist, "w:gz") as archive:
        archive.add(metadata_path, arcname="dotmatch-0.1.0/PKG-INFO")
    return sdist


def _write_sdist_with_egg_info_copy(path: Path, metadata: str = VALID_METADATA) -> Path:
    sdist = path / "dotmatch-0.1.0.tar.gz"
    root_metadata = path / "dotmatch-0.1.0" / "PKG-INFO"
    egg_metadata = path / "dotmatch-0.1.0" / "python" / "dotmatch.egg-info" / "PKG-INFO"
    root_metadata.parent.mkdir(parents=True)
    egg_metadata.parent.mkdir(parents=True)
    root_metadata.write_text(metadata, encoding="utf-8")
    egg_metadata.write_text(metadata.replace("Name: dotmatch", "Name: stale-copy"), encoding="utf-8")
    with tarfile.open(sdist, "w:gz") as archive:
        archive.add(root_metadata, arcname="dotmatch-0.1.0/PKG-INFO")
        archive.add(egg_metadata, arcname="dotmatch-0.1.0/python/dotmatch.egg-info/PKG-INFO")
    return sdist


def _write_wheel(path: Path, metadata: str = VALID_METADATA) -> Path:
    wheel = path / "dotmatch-0.1.0-py3-none-macosx_11_0_arm64.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr("dotmatch-0.1.0.dist-info/METADATA", metadata)
    return wheel


def test_distribution_metadata_accepts_valid_sdist_and_wheel_metadata(tmp_path):
    checker = _load_checker()
    sdist = _write_sdist(tmp_path)
    wheel = _write_wheel(tmp_path)

    checker.check_distribution_metadata(sdist, "0.1.0")
    checker.check_distribution_metadata(wheel, "0.1.0")


def test_distribution_metadata_uses_top_level_sdist_pkg_info(tmp_path):
    checker = _load_checker()
    sdist = _write_sdist_with_egg_info_copy(tmp_path)

    checker.check_distribution_metadata(sdist, "0.1.0")


def test_distribution_metadata_rejects_version_mismatch(tmp_path):
    checker = _load_checker()
    sdist = _write_sdist(tmp_path, VALID_METADATA.replace("Version: 0.1.0", "Version: 0.2.0"))

    try:
        checker.check_distribution_metadata(sdist, "0.1.0")
    except SystemExit as exc:
        assert "Version must be 0.1.0" in str(exc)
    else:
        raise AssertionError("expected metadata version mismatch to fail")


def test_distribution_metadata_requires_discovery_fields(tmp_path):
    checker = _load_checker()
    weak = (
        VALID_METADATA
        .replace("known-target short-DNA assignment", "fast matching")
        .replace("known-target assignment,", "")
        .replace("Classifier: Topic :: Scientific/Engineering :: Bio-Informatics\n", "")
        .replace("Project-URL: Repository, https://github.com/dnncha/dotmatch\n", "")
    )
    wheel = _write_wheel(tmp_path, weak)

    try:
        checker.check_distribution_metadata(wheel, "0.1.0")
    except SystemExit as exc:
        message = str(exc)
        assert "Summary must mention known-target short-DNA assignment" in message
        assert "Keywords must include known-target assignment" in message
        assert "Classifier must include Topic :: Scientific/Engineering :: Bio-Informatics" in message
        assert "Project-URL must include Repository" in message
    else:
        raise AssertionError("expected weak discovery metadata to fail")


def test_python_package_verifier_calls_metadata_checks() -> None:
    verifier = CHECKER.read_text(encoding="utf-8")

    assert "check_distribution_metadata(sdist, expected_version)" in verifier
    assert "check_distribution_metadata(wheel, expected_version)" in verifier


def test_build_and_verify_sdist_builds_sdist_without_wheel(tmp_path, monkeypatch):
    checker = _load_checker()
    calls: list[list[str]] = []

    def fake_run(cmd, *, cwd=None, env=None):
        calls.append(cmd)
        if cmd == [sys.executable, "-m", "build", "--sdist", "--outdir", str(tmp_path / "dist")]:
            (tmp_path / "dist").mkdir(exist_ok=True)
            _write_sdist(tmp_path / "dist")

    checked: list[str] = []
    monkeypatch.setattr(checker, "run", fake_run)
    monkeypatch.setattr(checker, "check_sdist_members", lambda artifact: checked.append(f"members:{artifact.name}"))
    monkeypatch.setattr(
        checker,
        "check_distribution_metadata",
        lambda artifact, version: checked.append(f"metadata:{artifact.name}:{version}"),
    )
    monkeypatch.setattr(
        checker,
        "verify_clean_install",
        lambda artifact, install_root, version: checked.append(f"install:{artifact.name}:{version}"),
    )

    checker.build_and_verify_sdist(tmp_path / "dist", tmp_path / "install", "0.1.0")

    assert calls == [[sys.executable, "-m", "build", "--sdist", "--outdir", str(tmp_path / "dist")]]
    assert checked == [
        "members:dotmatch-0.1.0.tar.gz",
        "metadata:dotmatch-0.1.0.tar.gz:0.1.0",
        "install:dotmatch-0.1.0.tar.gz:0.1.0",
    ]
