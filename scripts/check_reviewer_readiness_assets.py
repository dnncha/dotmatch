#!/usr/bin/env python3
"""Verify DotMatch reviewer materials are concrete, scoped, and public-facing."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from check_common import AuditResult, is_https_url, read_json, uses_placeholder_host


ROOT = Path(".")
READINESS = Path("docs/reviewer-readiness.json")
WORKFLOW_PLAN = Path("docs/workflow-integration-plan.json")
INTEGRATIONS = Path("docs/integration-targets.json")
BIOTOOLS = Path("docs/registries/biotools.yml")

EXPECTED_READINESS_IDS = {
    "bioinformatics_evaluation_packet",
    "external_review_packet",
    "integration_target_tracker",
    "biotools_registry_draft",
    "evaluation_protocol",
    "public_use_record_template",
    "workflow_integration_issue_template",
    "evaluation_feedback_issue_template",
    "pull_request_language_checklist",
    "reviewer_readiness_gate",
}

EXPECTED_WORKFLOW_PLAN_IDS = {
    "reviewer_decision_tree",
    "persona_one_pagers",
    "integration_tracker",
    "reviewer_packet",
    "conference_abstracts",
    "technical_communication_pack",
    "maintainer_issue_templates",
    "evaluation_scorecard",
    "integration_tracking_metrics",
    "release_publication_checklist",
}

EXPECTED_INTEGRATION_IDS = {
    "nf_core_modules",
    "multiqc_module",
    "galaxy_iuc",
    "snakemake_wrapper",
    "bio_tools_record",
}

ALLOWED_INTEGRATION_STATUSES = {
    "external_pr_open",
    "local_payload_ready",
    "local_parser_ready",
    "local_wrapper_ready",
    "local_workflow_ready",
    "draft_metadata_ready",
}

PUBLIC_LANGUAGE_FILES = [
    Path("README.md"),
    Path("app/page.tsx"),
    Path("app/layout.tsx"),
    Path("docs/bioinformatics-evaluation.md"),
    Path("docs/external-review-packet.md"),
    Path("docs/workflow-integration-kit.md"),
    Path("docs/workflow-integration-roadmap.md"),
    Path("docs/pilot-program.md"),
    Path("docs/adopters/README.md"),
    Path("docs/adopters/record-template.md"),
]

REQUIRED_FILES = [
    READINESS,
    WORKFLOW_PLAN,
    INTEGRATIONS,
    BIOTOOLS,
    Path("docs/bioinformatics-evaluation.md"),
    Path("docs/external-review-packet.md"),
    Path("docs/workflow-integration-kit.md"),
    Path("docs/workflow-integration-roadmap.md"),
    Path("docs/pilot-program.md"),
    Path("docs/adopters/record-template.md"),
    Path(".github/ISSUE_TEMPLATE/workflow_integration.yml"),
    Path(".github/ISSUE_TEMPLATE/pilot_feedback.yml"),
    Path("scripts/check_reviewer_readiness_assets.py"),
]

HYPE_PHRASES = [
    "massive industry impact",
    "ai slop",
    "game-changing",
    "revolutionary",
    "world-class",
    "best-in-class",
    "enterprise-grade",
    "production-ready",
    "just works",
    "magic",
]

INTERNAL_PROCESS_PHRASES = [
    "adoption evidence",
    "adoption trust",
    "evidence-bounded",
    "industry exposure",
    "next wins",
    "pilot conversations",
    "private feedback",
    "private pilot",
    "quote-approved",
    "turning private",
    "without turning",
    "community traction",
    "evaluator campaign",
]


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _require_file(path: Path, result: AuditResult) -> None:
    if not path.exists():
        result.failures.append(f"missing reviewer readiness asset: {path.as_posix()}")


def _check_item_plan(
    result: AuditResult,
    path: Path,
    expected_ids: set[str],
    label: str,
) -> None:
    data = read_json(path)
    if data.get("schema_version") != 1:
        result.failures.append(f"{label} must declare schema_version 1")
    if data.get("status") != "ready_to_execute":
        result.failures.append(f"{label} must stay ready_to_execute while public materials are complete")
    items = data.get("items")
    if not isinstance(items, list) or len(items) != len(expected_ids):
        result.failures.append(f"{label} must contain exactly {len(expected_ids)} items")
        return
    seen: set[str] = set()
    for item in items:
        if not isinstance(item, dict):
            result.failures.append(f"{label} items must be objects")
            continue
        item_id = str(item.get("id") or "").strip()
        asset = Path(str(item.get("asset") or "").split("#", 1)[0].strip())
        if item_id in seen:
            result.failures.append(f"duplicate {label} item id: {item_id}")
        seen.add(item_id)
        for field in ["id", "title", "asset", "done_when"]:
            if not str(item.get(field) or "").strip():
                result.failures.append(f"{label} item {item_id or '<missing>'} missing {field}")
        if asset and not asset.exists():
            result.failures.append(f"{label} item asset missing: {asset.as_posix()}")
    missing = expected_ids - seen
    extra = seen - expected_ids
    if missing:
        result.failures.append(f"{label} missing ids: {sorted(missing)}")
    if extra:
        result.failures.append(f"{label} has unexpected ids: {sorted(extra)}")


def _check_integrations(result: AuditResult) -> None:
    data = read_json(INTEGRATIONS)
    if data.get("schema_version") != 1:
        result.failures.append("integration target tracker must declare schema_version 1")
    if data.get("status") != "planned":
        result.failures.append("integration target tracker must remain planned until external records exist")
    if data.get("adoption_record") != "docs/workflow-adoption.json":
        result.failures.append("integration tracker must point at docs/workflow-adoption.json")
    targets = data.get("targets")
    if not isinstance(targets, list) or len(targets) != 5:
        result.failures.append("integration target tracker must contain five targets")
        return
    seen: set[str] = set()
    for target in targets:
        if not isinstance(target, dict):
            result.failures.append("integration targets must be objects")
            continue
        target_id = str(target.get("id") or "").strip()
        seen.add(target_id)
        status = str(target.get("status") or "").strip()
        if status not in ALLOWED_INTEGRATION_STATUSES:
            result.failures.append(f"{target_id} has unsupported planned status: {status}")
        url = str(target.get("external_target") or "").strip()
        if not is_https_url(url) or uses_placeholder_host(url):
            result.failures.append(f"{target_id} must declare a real https external_target")
        if target.get("public_record_required") is not True:
            result.failures.append(f"{target_id} must require a public record before integration claims")
        if not str(target.get("next_action") or "").strip():
            result.failures.append(f"{target_id} must declare next_action")
        for asset in target.get("source_assets") or []:
            asset_path = Path(str(asset))
            if not asset_path.exists():
                result.failures.append(f"{target_id} source asset missing: {asset_path.as_posix()}")
    missing = EXPECTED_INTEGRATION_IDS - seen
    extra = seen - EXPECTED_INTEGRATION_IDS
    if missing:
        result.failures.append(f"integration tracker missing ids: {sorted(missing)}")
    if extra:
        result.failures.append(f"integration tracker has unexpected ids: {sorted(extra)}")


def _check_biotools(result: AuditResult) -> None:
    text = _read(BIOTOOLS)
    for phrase in [
        "name: DotMatch",
        "biotoolsID: dotmatch",
        "documentation: https://dotmatch.readthedocs.io/",
        "license: Apache-2.0",
        "known-target short-DNA assignment",
        "Draft registry metadata",
    ]:
        if phrase not in text:
            result.failures.append(f"bio.tools draft missing: {phrase}")
    if "accepted bio.tools record" in text and "Do not describe" not in text:
        result.failures.append("bio.tools draft must not imply accepted registry status")


def _check_templates(result: AuditResult) -> None:
    workflow_issue = _read(Path(".github/ISSUE_TEMPLATE/workflow_integration.yml"))
    evaluation_issue = _read(Path(".github/ISSUE_TEMPLATE/pilot_feedback.yml"))
    pr_template = _read(Path(".github/PULL_REQUEST_TEMPLATE.md"))
    for phrase in ["Workflow Manager", "Expected Outputs", "Validated Scope"]:
        if phrase not in workflow_issue:
            result.failures.append(f"workflow integration template missing: {phrase}")
    for phrase in ["Assay Context", "Outputs Reviewed", "Public Use Permission"]:
        if phrase not in evaluation_issue:
            result.failures.append(f"evaluation feedback template missing: {phrase}")
    for phrase in [
        "Public Language and Reviewer Readiness",
        "This PR does not imply accepted external workflow integration",
        "Broad replacement wording or launch copy was removed or avoided",
    ]:
        if phrase not in pr_template:
            result.failures.append(f"PR template missing reviewer readiness checkbox: {phrase}")


def _check_docs_wiring(result: AuditResult) -> None:
    readme = _read(Path("README.md"))
    index = _read(Path("docs/index.md"))
    page = _read(Path("app/page.tsx"))
    public_text = readme + index + page
    for path in [
        "getting-started",
        "command-reference",
        "schemas",
        "methods-and-citation",
    ]:
        if path not in public_text:
            result.failures.append(f"public docs do not provide a user route to {path}")


def _check_public_language(result: AuditResult) -> None:
    forbidden = HYPE_PHRASES + INTERNAL_PROCESS_PHRASES
    for path in PUBLIC_LANGUAGE_FILES:
        text = _read(path).lower()
        for phrase in forbidden:
            if phrase in text:
                result.failures.append(f"{path.as_posix()} contains public-facing internal or hype phrase: {phrase}")


def audit(root: Path) -> AuditResult:
    result = AuditResult()
    for path in REQUIRED_FILES:
        _require_file(path, result)
    if result.failures:
        return result
    _check_item_plan(result, READINESS, EXPECTED_READINESS_IDS, "reviewer readiness record")
    _check_item_plan(result, WORKFLOW_PLAN, EXPECTED_WORKFLOW_PLAN_IDS, "workflow integration plan")
    _check_integrations(result)
    _check_biotools(result)
    _check_templates(result)
    _check_docs_wiring(result)
    _check_public_language(result)
    if result.ok:
        result.passed.append("reviewer materials are concrete, scoped, and public-facing")
    return result


def main() -> int:
    result = audit(ROOT)
    for item in result.passed:
        print(f"PASS: {item}")
    for item in result.failures:
        print(f"FAIL: {item}")
    if result.ok:
        print("REVIEWER READINESS: PASS")
        return 0
    print("REVIEWER READINESS: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
