"""Offline checks for publication preflights. No GitHub calls or remote writes."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import publish_github as publisher


def assets(directory):
    paths = [directory / "editwitness-0.2.0a1-py3-none-any.whl",
             directory / "editwitness-0.2.0a1.tar.gz"]
    for p in paths:
        p.write_bytes(b"synthetic test asset")
    (directory / "SHA256SUMS").write_text("".join(
        f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n" for p in paths
    ), encoding="ascii")
    return paths


def test_release_artifact_hashes_and_inventory_must_match(tmp_path):
    paths = assets(tmp_path)
    assert len(publisher.verify_artifacts(tmp_path)) == 3
    paths[0].write_bytes(b"tampered")
    with pytest.raises(ValueError, match="checksum mismatch"):
        publisher.verify_artifacts(tmp_path)


def test_release_refuses_extra_assets_and_unsafe_paths(tmp_path):
    assets(tmp_path)
    extra = tmp_path / "unreviewed.txt"
    extra.write_text("no", encoding="utf-8")
    with pytest.raises(ValueError, match="unlisted"):
        publisher.verify_artifacts(tmp_path)
    extra.unlink()
    (tmp_path / "SHA256SUMS").write_text("a" * 64 + "  ../unsafe.whl\n", encoding="ascii")
    with pytest.raises(ValueError, match="invalid"):
        publisher.verify_artifacts(tmp_path)


def test_publisher_accepts_only_explicit_alpha_versions(tmp_path):
    path = tmp_path / "pyproject.toml"
    path.write_text('[project]\nversion="0.2.0a1"\n', encoding="utf-8")
    assert publisher.checked_version(tmp_path) == "0.2.0a1"
    path.write_text('[project]\nversion="1.0.0"\n', encoding="utf-8")
    with pytest.raises(ValueError, match="alpha"):
        publisher.checked_version(tmp_path)


def test_ci_gate_watches_only_matching_push_commit(monkeypatch):
    calls = []
    def run(*args, **kwargs):
        calls.append(args)
        if args[1:3] == ("run", "list"):
            data = [{"databaseId": 1, "headSha": "other", "event": "push"},
                    {"databaseId": 2, "headSha": "target", "event": "pull_request"},
                    {"databaseId": 3, "headSha": "target", "event": "push"}]
            return subprocess.CompletedProcess(args, 0, json.dumps(data), "")
        return subprocess.CompletedProcess(args, 0, "", "")
    monkeypatch.setattr(publisher, "command", run)
    assert publisher.wait_for_ci("example/editwitness", "target") == 3
    assert calls[-1][:4] == ("gh", "run", "watch", "3")
    assert "--exit-status" in calls[-1]


def test_failed_ci_aborts_release_gate(monkeypatch):
    def run(*args, **kwargs):
        if args[1:3] == ("run", "list"):
            return subprocess.CompletedProcess(args, 0, json.dumps([
                {"databaseId": 3, "headSha": "target", "event": "push"}]), "")
        raise subprocess.CalledProcessError(1, args, stderr="tests failed")
    monkeypatch.setattr(publisher, "command", run)
    with pytest.raises(subprocess.CalledProcessError):
        publisher.wait_for_ci("example/editwitness", "target")


def test_all_actions_are_pinned_and_permissions_are_read_only():
    import re
    root = SCRIPTS.parent / ".github/workflows"
    for path in root.glob("*.yml"):
        text = path.read_text(encoding="utf-8")
        for action in re.findall(r"uses:\s*(\S+)", text):
            assert re.fullmatch(r"[\w.-]+/[\w.-]+@[a-f0-9]{40}", action)
        assert "contents: read" in text
        assert "contents: write" not in text


def test_release_artifact_version_matches_reviewed_source(tmp_path):
    assets(tmp_path)
    assert len(publisher.verify_artifacts(tmp_path, expected_version="0.2.0a1")) == 3
    with pytest.raises(ValueError, match="reviewed package version"):
        publisher.verify_artifacts(tmp_path, expected_version="0.2.0a2")


def test_release_refuses_multiple_wheels_even_with_valid_checksums(tmp_path):
    assets(tmp_path)
    extra = tmp_path / "other.whl"
    extra.write_bytes(b"different")
    checksum = tmp_path / "SHA256SUMS"
    with checksum.open("a", encoding="ascii") as handle:
        handle.write(f"{hashlib.sha256(extra.read_bytes()).hexdigest()}  {extra.name}\n")
    with pytest.raises(ValueError, match="exactly one"):
        publisher.verify_artifacts(tmp_path)
