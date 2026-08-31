#!/usr/bin/env python3
"""Validate the structured public distribution release record."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from check_common import AuditResult, check_https_url, check_simple_make_target, read_json


MANIFEST = Path("docs") / "distribution-release.json"
REQUIRED_CHANNELS = ["pypi", "bioconda", "ghcr", "biocontainers", "zenodo"]
VALID_OVERALL_STATUSES = {"not_released", "partially_verified", "released"}
VALID_CHANNEL_STATUSES = {"prepared", "blocked", "manifest_verified", "verified"}
VERIFIED_CHANNEL_STATUSES = {"manifest_verified", "verified"}
VALID_CANDIDATE_STATUSES = {"prepared_not_published"}
SUPPORTED_PYPI_LINUX_WHEEL_ARCHITECTURES = {"x86_64", "aarch64"}
SUPPORTED_GHCR_PLATFORMS = {"linux/amd64", "linux/arm64"}


def _project_version(root: Path) -> str:
    text = (root / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    return match.group(1) if match else ""


def _check_https_url(channel_id: str, field: str, value: str, result: AuditResult) -> bool:
    return check_https_url(channel_id, field, value, result)


def _check_declared_values(
    channel_id: str,
    item: dict[str, object],
    field: str,
    supported_values: set[str],
    required_value: str,
    result: AuditResult,
) -> None:
    value = item.get(field)
    if not isinstance(value, list) or not value or not all(isinstance(entry, str) and entry for entry in value):
        result.failures.append(f"{channel_id} must declare {field} as a non-empty list of strings")
        return
    entries = [str(entry) for entry in value]
    if len(entries) != len(set(entries)):
        result.failures.append(f"{channel_id} {field} must not contain duplicates")
    invalid = sorted(set(entries) - supported_values)
    if invalid:
        result.failures.append(f"{channel_id} {field} contains unsupported value(s): {', '.join(invalid)}")
    if required_value not in entries:
        result.failures.append(f"{channel_id} {field} must include {required_value}")


def _check_channel(item: object, overall_status: str, result: AuditResult) -> str:
    if not isinstance(item, dict):
        result.failures.append("distribution channels must be objects")
        return ""
    channel_id = str(item.get("id") or "").strip()
    if not channel_id:
        result.failures.append("distribution channel missing id")
        return ""
    status = str(item.get("status") or "").strip()
    if status not in VALID_CHANNEL_STATUSES:
        result.failures.append(f"{channel_id} has invalid distribution channel status: {status}")
    if status in VERIFIED_CHANNEL_STATUSES:
        if channel_id == "pypi":
            _check_declared_values(
                channel_id,
                item,
                "linux_wheel_architectures",
                SUPPORTED_PYPI_LINUX_WHEEL_ARCHITECTURES,
                "x86_64",
                result,
            )
        elif channel_id == "ghcr":
            _check_declared_values(
                channel_id,
                item,
                "platforms",
                SUPPORTED_GHCR_PLATFORMS,
                "linux/amd64",
                result,
            )

    verification = str(item.get("verification_command") or "").strip()
    if verification != "make distribution-channels":
        result.failures.append(f"{channel_id} must use make distribution-channels as verification_command")

    expected_url = str(item.get("expected_url") or "").strip()
    _check_https_url(channel_id, "expected_url", expected_url, result)

    if overall_status == "released":
        if status != "verified":
            result.failures.append(f"released channel {channel_id} must be verified")
        public_url = str(item.get("public_url") or "").strip()
        evidence_url = str(item.get("evidence_url") or "").strip()
        verified_date = str(item.get("verified_date") or "").strip()
        _check_https_url(channel_id, "public_url", public_url, result)
        _check_https_url(channel_id, "evidence_url", evidence_url, result)
        if not re.match(r"^\d{4}-\d{2}-\d{2}$", verified_date):
            result.failures.append(f"{channel_id} must declare verified_date as YYYY-MM-DD")
        if item.get("blocker"):
            result.failures.append(f"released channel {channel_id} must not keep blocker text")
        if item.get("next_action"):
            result.failures.append(f"released channel {channel_id} must not keep next_action text")
    elif overall_status == "partially_verified":
        if status == "verified":
            public_url = str(item.get("public_url") or item.get("expected_url") or "").strip()
            evidence_url = str(item.get("evidence_url") or item.get("expected_url") or "").strip()
            _check_https_url(channel_id, "public_url", public_url, result)
            _check_https_url(channel_id, "evidence_url", evidence_url, result)
            if item.get("blocker"):
                result.failures.append(f"verified channel {channel_id} must not keep blocker text")
        else:
            if not str(item.get("blocker") or "").strip():
                result.failures.append(f"{channel_id} must declare blocker while partially verified")
            if not str(item.get("next_action") or "").strip():
                result.failures.append(f"{channel_id} must declare next_action while partially verified")
    else:
        if status == "verified":
            result.failures.append(f"not_released channel {channel_id} must not be marked verified")
        if not str(item.get("blocker") or "").strip():
            result.failures.append(f"{channel_id} must declare blocker while not released")
        if not str(item.get("next_action") or "").strip():
            result.failures.append(f"{channel_id} must declare next_action while not released")
        if channel_id == "pypi":
            pypi_text = f"{item.get('blocker') or ''} {item.get('next_action') or ''}"
            if "source distribution" not in pypi_text or "repaired" not in pypi_text or "wheel" not in pypi_text:
                result.failures.append("pypi not_released record must mention source distribution and repaired wheel publication")
    return channel_id


def audit(root: Path) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    try:
        manifest = read_json(root / MANIFEST)
    except Exception as exc:
        result.failures.append(f"{MANIFEST.as_posix()} could not be read: {exc}")
        return result

    if manifest.get("schema_version") != 1:
        result.failures.append("distribution release record must declare schema_version 1")
    status = str(manifest.get("status") or "").strip()
    if status not in VALID_OVERALL_STATUSES:
        result.failures.append("distribution release status must be not_released or released")

    release_version = str(manifest.get("release_version") or "").strip()
    project_version = _project_version(root)
    if not release_version:
        result.failures.append("distribution release record must declare release_version")
    elif project_version and release_version != project_version:
        candidate_version = str(manifest.get("candidate_version") or "").strip()
        candidate_status = str(manifest.get("candidate_status") or "").strip()
        publication_authorized = manifest.get("publication_authorized")
        if candidate_version != project_version:
            result.failures.append(
                f"release_version {release_version} does not match pyproject version {project_version}; "
                "candidate_version must identify the local package candidate"
            )
        if candidate_status not in VALID_CANDIDATE_STATUSES:
            result.failures.append("an unreleased package candidate must declare candidate_status prepared_not_published")
        if publication_authorized is not False:
            result.failures.append("an unreleased package candidate must declare publication_authorized false")

    gate = str(manifest.get("post_release_gate") or "").strip()
    if gate != "make distribution-channels":
        result.failures.append("distribution release record must use make distribution-channels as post_release_gate")
    elif root.joinpath("Makefile").exists():
        check_simple_make_target(root, gate, "distribution release post_release_gate", result)

    channels = manifest.get("channels")
    if not isinstance(channels, list):
        result.failures.append("distribution release record must contain channels list")
        return result
    seen: set[str] = set()
    for item in channels:
        channel_id = _check_channel(item, status, result)
        if channel_id and channel_id in seen:
            result.failures.append(f"duplicate distribution channel id: {channel_id}")
        seen.add(channel_id)
    seen.discard("")
    for channel_id in REQUIRED_CHANNELS:
        if channel_id not in seen:
            result.failures.append(f"missing required distribution channel: {channel_id}")

    if status not in VALID_OVERALL_STATUSES:
        return result

    if status == "released":
        if manifest.get("blockers"):
            result.failures.append("released distribution record must not declare blockers")
        if manifest.get("next_action"):
            result.failures.append("released distribution record must not declare next_action")
    elif status == "partially_verified":
        if not manifest.get("blockers"):
            result.failures.append("partially_verified distribution record must declare blockers")
        if not manifest.get("next_action"):
            result.failures.append("partially_verified distribution record must declare next_action")
    elif status == "not_released":
        if not manifest.get("blockers"):
            result.failures.append("not_released distribution record must declare blockers")
        if not manifest.get("next_action"):
            result.failures.append("not_released distribution record must declare next_action")

    if result.ok:
        result.passed.append("distribution release record valid")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="repository root")
    args = parser.parse_args()

    result = audit(Path(args.root))
    for item in result.passed:
        print(f"PASS: {item}")
    for item in result.failures:
        print(f"FAIL: {item}")
    if result.ok:
        print("DISTRIBUTION RECORD: PASS")
        return 0
    print("DISTRIBUTION RECORD: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
