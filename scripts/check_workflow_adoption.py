#!/usr/bin/env python3
"""Verify recorded external workflow adoption links for DotMatch.

This is a post-adoption checker. It is expected to fail until at least one real
external workflow integration has landed and been recorded in
docs/workflow-adoption.json.
"""

from __future__ import annotations

import argparse
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from check_common import AuditResult, is_https_url, read_json, uses_placeholder_host


MANIFEST = Path("docs") / "workflow-adoption.json"
VALID_STATUSES = {"not_ready", "ready"}
VALID_INTEGRATION_STATUSES = {"accepted", "released", "published"}
VALID_INTEGRATION_TYPES = {
    "nf_core_module",
    "nextflow_pipeline",
    "snakemake_workflow",
    "galaxy_toolshed",
    "multiqc_plugin",
    "independent_workflow",
    "bioconda_recipe",
}

def url_ok(url: str) -> bool:
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "DotMatch workflow adoption verifier"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return 200 <= int(response.status) < 400
    except urllib.error.HTTPError as exc:
        if exc.code not in {405, 501}:
            return False
    except Exception:
        return False
    request = urllib.request.Request(url, headers={"User-Agent": "DotMatch workflow adoption verifier"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return 200 <= int(response.status) < 400
    except Exception:
        return False


def _load_manifest(root: Path) -> dict:
    return read_json(root / MANIFEST)


def _check_required_text(item: dict, field: str, result: AuditResult) -> str:
    value = str(item.get(field) or "").strip()
    if not value:
        result.failures.append(f"{item.get('id', '<missing id>')} must declare {field}")
    return value


def _check_external_url(item: dict, field: str, result: AuditResult) -> None:
    value = _check_required_text(item, field, result)
    if not value:
        return
    if not is_https_url(value):
        result.failures.append(f"{item.get('id', '<missing id>')} {field} must be an https URL")
        return
    if uses_placeholder_host(value):
        result.failures.append(f"{item.get('id', '<missing id>')} {field} must not use placeholder domains")
        return
    if not url_ok(value):
        result.failures.append(f"{item.get('id', '<missing id>')} {field} is not reachable: {value}")


def check_integration(item: object, result: AuditResult) -> None:
    if not isinstance(item, dict):
        result.failures.append("workflow adoption integrations must be objects")
        return
    integration_id = _check_required_text(item, "id", result)
    integration_type = _check_required_text(item, "type", result)
    _check_required_text(item, "name", result)
    integration_status = _check_required_text(item, "status", result)
    _check_required_text(item, "validation_notes", result)
    recorded_date = _check_required_text(item, "recorded_date", result)
    if integration_type and integration_type not in VALID_INTEGRATION_TYPES:
        result.failures.append(f"{integration_id or '<missing id>'} has unsupported integration type: {integration_type}")
    if integration_status and integration_status not in VALID_INTEGRATION_STATUSES:
        result.failures.append(f"{integration_id or '<missing id>'} has invalid integration status: {integration_status}")
    if recorded_date and not re.match(r"^\d{4}-\d{2}-\d{2}$", recorded_date):
        result.failures.append(f"{integration_id or '<missing id>'} recorded_date must be YYYY-MM-DD")
    _check_external_url(item, "adoption_url", result)
    _check_external_url(item, "evidence_url", result)


def audit(root: Path) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    try:
        manifest = _load_manifest(root)
    except Exception as exc:
        result.failures.append(f"{MANIFEST.as_posix()} could not be read: {exc}")
        return result

    if manifest.get("schema_version") != 1:
        result.failures.append("workflow adoption manifest must declare schema_version 1")
    status = manifest.get("status")
    if status not in VALID_STATUSES:
        result.failures.append("workflow adoption status must be not_ready or ready")

    integrations = manifest.get("integrations")
    if not isinstance(integrations, list):
        result.failures.append("workflow adoption manifest must contain an integrations list")
        return result
    seen_ids: set[str] = set()
    for item in integrations:
        if isinstance(item, dict):
            integration_id = str(item.get("id") or "").strip()
            if integration_id and integration_id in seen_ids:
                result.failures.append(f"duplicate workflow adoption integration id: {integration_id}")
            seen_ids.add(integration_id)
        check_integration(item, result)

    if status == "ready" and not integrations:
        result.failures.append("ready status requires at least one external integration")
    if status == "ready":
        if manifest.get("blockers"):
            result.failures.append("ready workflow adoption manifest must not declare blockers")
        if manifest.get("next_action"):
            result.failures.append("ready workflow adoption manifest must not declare next_action")
    if status == "not_ready":
        if not manifest.get("blockers"):
            result.failures.append("not_ready workflow adoption manifest must declare blockers")
        if not manifest.get("next_action"):
            result.failures.append("not_ready workflow adoption manifest must declare next_action")
        result.failures.append("external workflow adoption is not recorded")

    if status == "ready" and not result.failures:
        result.passed.append("external workflow adoption recorded")
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
        print("WORKFLOW ADOPTION: PASS")
        return 0
    print("WORKFLOW ADOPTION: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
