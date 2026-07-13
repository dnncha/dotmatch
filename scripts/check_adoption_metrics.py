#!/usr/bin/env python3
"""Validate the adoption measurement contract without collecting user data."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from check_common import AuditResult, is_https_url, read_json


MANIFEST = Path("docs/adoption-metrics.json")
REQUIRED_METRIC_IDS = {
    "anaconda_downloads_6m",
    "pypi_downloads_30d",
    "completed_external_evaluations_30d",
    "repeat_workflows_90d",
    "accepted_external_integrations",
    "docs_to_install_intent",
}
VALID_ROLES = {"distribution_health", "north_star", "retention", "ecosystem", "funnel"}


def audit(root: Path) -> AuditResult:
    result = AuditResult()
    path = root / MANIFEST
    try:
        data = read_json(path)
    except Exception as exc:
        result.failures.append(f"{MANIFEST.as_posix()} could not be read: {exc}")
        return result

    if data.get("schema_version") != 1:
        result.failures.append("adoption metrics must declare schema_version 1")
    if data.get("status") != "instrumentation_ready":
        result.failures.append("adoption metrics must remain instrumentation_ready until a live dashboard exists")

    north_star = data.get("north_star")
    if not isinstance(north_star, dict):
        result.failures.append("adoption metrics must declare a north_star object")
    else:
        for field in ["id", "definition", "target_next_90_days", "why"]:
            if not str(north_star.get(field) or "").strip() and north_star.get(field) != 0:
                result.failures.append(f"north_star missing {field}")
        if north_star.get("id") != "completed_external_evaluations_30d":
            result.failures.append("north_star must be completed_external_evaluations_30d")

    metrics = data.get("metrics")
    if not isinstance(metrics, list) or len(metrics) != len(REQUIRED_METRIC_IDS):
        result.failures.append(f"adoption metrics must contain exactly {len(REQUIRED_METRIC_IDS)} metrics")
        metrics = []
    seen: set[str] = set()
    for metric in metrics:
        if not isinstance(metric, dict):
            result.failures.append("adoption metrics entries must be objects")
            continue
        metric_id = str(metric.get("id") or "").strip()
        if metric_id in seen:
            result.failures.append(f"duplicate adoption metric id: {metric_id}")
        seen.add(metric_id)
        for field in ["id", "name", "definition", "source", "cadence", "target", "limitation"]:
            if not str(metric.get(field) or "").strip():
                result.failures.append(f"{metric_id or '<missing id>'} missing {field}")
        role = str(metric.get("role") or "").strip()
        if role not in VALID_ROLES:
            result.failures.append(f"{metric_id or '<missing id>'} has unsupported role: {role}")
        source = str(metric.get("source") or "").strip()
        if source.startswith("http") and not is_https_url(source):
            result.failures.append(f"{metric_id or '<missing id>'} source must use https")
    missing = REQUIRED_METRIC_IDS - seen
    extra = seen - REQUIRED_METRIC_IDS
    if missing:
        result.failures.append(f"adoption metrics missing ids: {sorted(missing)}")
    if extra:
        result.failures.append(f"adoption metrics has unexpected ids: {sorted(extra)}")

    rules = data.get("decision_rules")
    if not isinstance(rules, list) or len(rules) < 3:
        result.failures.append("adoption metrics must declare at least three decision_rules")
    joined = " ".join(str(rule) for rule in rules or []).lower()
    for phrase in ["download count", "unique-user", "public use"]:
        if phrase not in joined:
            result.failures.append(f"adoption metrics decision_rules must mention {phrase}")
    if not result.failures:
        result.passed.append("adoption measurement contract is explicit and non-identifying")
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
        print("ADOPTION METRICS: PASS")
        return 0
    print("ADOPTION METRICS: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
