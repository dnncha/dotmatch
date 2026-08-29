from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_oci_manifest.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_oci_manifest", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _manifest(*platforms: str) -> dict[str, object]:
    descriptors: list[dict[str, object]] = []
    for platform in platforms:
        operating_system, architecture = platform.split("/", 1)
        descriptors.append(
            {
                "mediaType": "application/vnd.oci.image.manifest.v1+json",
                "digest": f"sha256:{architecture}",
                "platform": {"os": operating_system, "architecture": architecture},
            }
        )
    return {"schemaVersion": 2, "manifests": descriptors}


def test_manifest_checker_accepts_required_linux_platforms() -> None:
    checker = _load_checker()
    manifest = _manifest("linux/amd64", "linux/arm64", "unknown/unknown")

    platforms = checker.check_manifest(manifest, ["linux/amd64", "linux/arm64"])

    assert platforms == {"linux/amd64", "linux/arm64", "unknown/unknown"}


def test_manifest_checker_rejects_missing_required_platform() -> None:
    checker = _load_checker()

    try:
        checker.check_manifest(_manifest("linux/amd64"), ["linux/amd64", "linux/arm64"])
    except ValueError as exc:
        assert "linux/arm64" in str(exc)
    else:
        raise AssertionError("expected missing linux/arm64 platform to fail")


def test_manifest_checker_requires_an_image_index_or_manifest_list() -> None:
    checker = _load_checker()

    try:
        checker.check_manifest({"schemaVersion": 2, "config": {}}, ["linux/amd64"])
    except ValueError as exc:
        assert "image index or manifest list" in str(exc)
    else:
        raise AssertionError("expected single-image manifest to fail")
