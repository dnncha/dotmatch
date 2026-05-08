import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_distribution_record.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_distribution_record", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _channel(channel_id: str, **overrides) -> dict:
    expected_urls = {
        "pypi": "https://pypi.org/project/dotmatch/0.1.0/",
        "bioconda": "https://anaconda.org/bioconda/dotmatch",
        "ghcr": "https://github.com/dnncha/dotmatch/pkgs/container/dotmatch",
        "biocontainers": "https://quay.io/repository/biocontainers/dotmatch",
        "zenodo": "https://zenodo.org/records/123456",
    }
    item = {
        "id": channel_id,
        "status": "prepared",
        "expected_url": expected_urls.get(channel_id, f"https://github.com/dnncha/dotmatch/{channel_id}"),
        "verification_command": "make distribution-channels",
        "blocker": f"{channel_id} is not public yet.",
        "next_action": f"Publish {channel_id} and rerun make distribution-channels.",
    }
    if channel_id == "pypi":
        item["blocker"] = "The source distribution and repaired manylinux/musllinux wheels are not public yet."
        item["next_action"] = "Publish the source distribution and repaired Linux wheels, then rerun make distribution-channels."
    item.update(overrides)
    return item


def _manifest(status: str = "not_released", channels=None) -> dict:
    return {
        "schema_version": 1,
        "status": status,
        "release_version": "0.1.0",
        "post_release_gate": "make distribution-channels",
        "channels": channels
        if channels is not None
        else [
            _channel("pypi"),
            _channel("bioconda"),
            _channel("ghcr"),
            _channel("biocontainers"),
            _channel("zenodo"),
        ],
        "blockers": ["Public release has not been executed yet."],
        "next_action": "Tag and publish the release, then replace expected URLs with verified public links.",
    }


def _write_repo(root: Path, manifest=None) -> None:
    files = {
        "pyproject.toml": '[project]\nname = "dotmatch"\nversion = "0.1.0"\n',
        "docs/distribution-release.json": json.dumps(manifest or _manifest(), indent=2) + "\n",
        "Makefile": "distribution-channels:\n\ttrue\n",
    }
    for path, text in files.items():
        full = root / path
        full.parent.mkdir(parents=True, exist_ok=True)
        full.write_text(text, encoding="utf-8")


def test_distribution_record_accepts_not_released_manifest(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)

    result = checker.audit(tmp_path)

    assert result.failures == []
    assert any("distribution release record valid" in item for item in result.passed)


def test_distribution_record_accepts_released_manifest_with_verified_links(tmp_path):
    checker = _load_checker()
    manifest = _manifest(
        status="released",
        channels=[
            _channel(
                channel_id,
                status="verified",
                public_url=f"https://github.com/dnncha/dotmatch/releases/tag/v0.1.0-{channel_id}",
                evidence_url=f"https://github.com/dnncha/dotmatch/actions/runs/123456#{channel_id}",
                verified_date="2026-05-07",
                blocker="",
                next_action="",
            )
            for channel_id in ["pypi", "bioconda", "ghcr", "biocontainers", "zenodo"]
        ],
    )
    manifest["blockers"] = []
    manifest["next_action"] = ""
    _write_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert result.failures == []


def test_distribution_record_rejects_missing_required_channel(tmp_path):
    checker = _load_checker()
    manifest = _manifest(channels=[_channel("pypi")])
    _write_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("missing required distribution channel: bioconda" in failure for failure in result.failures)


def test_distribution_record_rejects_released_manifest_without_verified_links(tmp_path):
    checker = _load_checker()
    manifest = _manifest(status="released")
    manifest["blockers"] = []
    manifest["next_action"] = ""
    _write_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("released channel pypi must be verified" in failure for failure in result.failures)
    assert any("pypi must declare public_url" in failure for failure in result.failures)


def test_distribution_record_rejects_released_channel_with_stale_next_action(tmp_path):
    checker = _load_checker()
    manifest = _manifest(
        status="released",
        channels=[
            _channel(
                channel_id,
                status="verified",
                public_url=f"https://github.com/dnncha/dotmatch/releases/tag/v0.1.0-{channel_id}",
                evidence_url=f"https://github.com/dnncha/dotmatch/actions/runs/123456#{channel_id}",
                verified_date="2026-05-07",
                blocker="",
                next_action="Publish this channel later.",
            )
            for channel_id in ["pypi", "bioconda", "ghcr", "biocontainers", "zenodo"]
        ],
    )
    manifest["blockers"] = []
    manifest["next_action"] = ""
    _write_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("released channel pypi must not keep next_action text" in failure for failure in result.failures)


def test_distribution_record_rejects_not_released_manifest_without_blockers(tmp_path):
    checker = _load_checker()
    manifest = _manifest()
    manifest["blockers"] = []
    manifest["next_action"] = ""
    manifest["channels"][0]["blocker"] = ""
    _write_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("not_released distribution record must declare blockers" in failure for failure in result.failures)
    assert any("not_released distribution record must declare next_action" in failure for failure in result.failures)
    assert any("pypi must declare blocker while not released" in failure for failure in result.failures)


def test_distribution_record_requires_pypi_repaired_wheel_publication_scope(tmp_path):
    checker = _load_checker()
    manifest = _manifest()
    manifest["channels"][0]["blocker"] = "The source distribution is not public yet."
    manifest["channels"][0]["next_action"] = "Publish the source distribution."
    _write_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("source distribution and repaired wheel publication" in failure for failure in result.failures)


def test_distribution_record_rejects_version_mismatch(tmp_path):
    checker = _load_checker()
    manifest = _manifest()
    manifest["release_version"] = "0.2.0"
    _write_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("release_version 0.2.0 does not match pyproject version 0.1.0" in failure for failure in result.failures)


def test_distribution_record_rejects_duplicate_channel_ids(tmp_path):
    checker = _load_checker()
    manifest = _manifest(channels=[_channel("pypi"), _channel("pypi")])
    _write_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("duplicate distribution channel id: pypi" in failure for failure in result.failures)


def test_distribution_record_rejects_placeholder_public_and_evidence_urls(tmp_path):
    checker = _load_checker()
    manifest = _manifest(
        status="released",
        channels=[
            _channel(
                channel_id,
                status="verified",
                public_url=f"https://example.org/{channel_id}/dotmatch/0.1.0",
                evidence_url=f"https://example.org/{channel_id}/evidence",
                verified_date="2026-05-07",
                blocker="",
                next_action="",
            )
            for channel_id in ["pypi", "bioconda", "ghcr", "biocontainers", "zenodo"]
        ],
    )
    manifest["blockers"] = []
    manifest["next_action"] = ""
    _write_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("pypi public_url must not use placeholder domains" in failure for failure in result.failures)
    assert any("pypi evidence_url must not use placeholder domains" in failure for failure in result.failures)


def test_distribution_record_rejects_invalid_verified_date(tmp_path):
    checker = _load_checker()
    manifest = _manifest(
        status="released",
        channels=[
            _channel(
                channel_id,
                status="verified",
                public_url=f"https://github.com/dnncha/dotmatch/releases/tag/v0.1.0-{channel_id}",
                evidence_url=f"https://github.com/dnncha/dotmatch/actions/runs/123456#{channel_id}",
                verified_date="May 7 2026",
                blocker="",
                next_action="",
            )
            for channel_id in ["pypi", "bioconda", "ghcr", "biocontainers", "zenodo"]
        ],
    )
    manifest["blockers"] = []
    manifest["next_action"] = ""
    _write_repo(tmp_path, manifest)

    result = checker.audit(tmp_path)

    assert any("pypi must declare verified_date as YYYY-MM-DD" in failure for failure in result.failures)
