from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from scripts.check_assaycode_bioconda_recipe import audit
from scripts.prepare_bioconda_handoff import PLACEHOLDER, render


def test_assaycode_metapackage_contract() -> None:
    assert audit() == []


def test_handoff_renders_both_recipes_with_real_checksum(tmp_path: Path) -> None:
    archive = tmp_path / "v0.2.2.tar.gz"
    archive.write_bytes(b"immutable release fixture")

    dotmatch_dir, assaycode_dir, digest = render(archive, tmp_path / "handoff")

    expected = hashlib.sha256(archive.read_bytes()).hexdigest()
    dotmatch_meta = (dotmatch_dir / "meta.yaml").read_text(encoding="utf-8")
    assaycode_meta = (assaycode_dir / "meta.yaml").read_text(encoding="utf-8")
    assert digest == expected
    assert expected in dotmatch_meta
    assert PLACEHOLDER not in dotmatch_meta
    assert (dotmatch_dir / "build.sh").is_file()
    assert "name: assaycode" in assaycode_meta
    assert "- dotmatch =={{ version }}" in assaycode_meta


def test_handoff_rejects_archive_for_another_version(tmp_path: Path) -> None:
    archive = tmp_path / "v9.9.9.tar.gz"
    archive.write_bytes(b"wrong release")
    with pytest.raises(ValueError, match="must identify 0.2.2"):
        render(archive, tmp_path / "handoff")
