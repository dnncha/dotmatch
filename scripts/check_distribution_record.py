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
VALID_OVERALL_STATUSES = {"not_released", "released"}
VALID_CHANNEL_STATUSES = {"prepared", "blocked", "verified"}


def _project_version(root: Path) -> str:
    text = (root / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    return match.group(1) if match else ""


def _check_https_url(channel_id: str, field: str, value: str, result: AuditResult) -> bool:
    return check_https_url(channel_id, field, value, result)


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
        result.failures.append(f"release_version {release_version} does not match pyproject version {project_version}")

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

    if status == "released":
        if manifest.get("blockers"):
            result.failures.append("released distribution record must not declare blockers")
        if manifest.get("next_action"):
            result.failures.append("released distribution record must not declare next_action")
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
