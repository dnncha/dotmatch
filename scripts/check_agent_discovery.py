#!/usr/bin/env python3
"""Validate and measure DotMatch's agent discovery and task-routing surfaces."""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REQUIRED_INTENTS = {
    "crispr-guide-counting",
    "inline-barcode-demultiplexing",
    "feature-barcode-assignment",
    "perturb-seq-guide-capture",
    "barcode-panel-design",
    "known-target-fastq-matching",
}
REQUIRED_OUTCOMES = {"unique", "ambiguous", "none", "invalid"}
REQUIRED_SEARCH_TERMS = {
    "barcode-demultiplexing",
    "crispr",
    "fastq",
    "feature-barcodes",
    "guide-capture",
    "guide-counting",
    "perturb-seq",
}
PUBLISHED_COPIES = {
    "llms.txt": ["docs/llms.txt", "public/llms.txt"],
    "llms-full.txt": ["docs/llms-full.txt", "public/llms-full.txt"],
    "agent-capabilities.json": [
        "docs/agent-capabilities.json",
        "public/agent-capabilities.json",
        "python/dotmatch/data/agent-capabilities.json",
    ],
    "agent-capabilities.schema.json": [
        "docs/agent-capabilities.schema.json",
        "public/agent-capabilities.schema.json",
        "python/dotmatch/data/agent-capabilities.schema.json",
    ],
    "agent-capabilities.v1.json": [
        "docs/agent-capabilities.v1.json",
        "public/agent-capabilities.v1.json",
    ],
    "agent-capabilities.v1.schema.json": [
        "docs/agent-capabilities.v1.schema.json",
        "public/agent-capabilities.v1.schema.json",
    ],
    "agent-tools.json": [
        "docs/agent-tools.json",
        "public/agent-tools.json",
        "python/dotmatch/data/agent-tools.json",
    ],
    "agent-tools.schema.json": [
        "docs/agent-tools.schema.json",
        "public/agent-tools.schema.json",
        "python/dotmatch/data/agent-tools.schema.json",
    ],
    "agent-reference-crispr.json": [
        "docs/agent-reference-crispr.json",
        "public/agent-reference-crispr.json",
    ],
}


def _read(root: Path, relative: str) -> str:
    return (root / relative).read_text(encoding="utf-8")


def _json(root: Path, relative: str) -> dict[str, Any]:
    value = json.loads(_read(root, relative))
    if not isinstance(value, dict):
        raise ValueError(f"{relative} must contain a JSON object")
    return value


def _project_version(root: Path) -> str:
    match = re.search(r'^version\s*=\s*"([^"]+)"', _read(root, "pyproject.toml"), flags=re.M)
    return match.group(1) if match else ""


def validate_manifest(root: Path) -> list[str]:
    failures: list[str] = []
    try:
        manifest = _json(root, "agent-capabilities.json")
    except Exception as exc:
        return [f"agent-capabilities.json is invalid: {exc}"]

    required_top = {
        "$schema",
        "schema_version",
        "generated_for_version",
        "project",
        "scope",
        "install",
        "outcomes",
        "intents",
        "error_recovery",
        "interfaces",
    }
    missing_top = sorted(required_top - set(manifest))
    failures.extend(f"agent-capabilities.json missing top-level field: {field}" for field in missing_top)

    if manifest.get("schema_version") != "1.1":
        failures.append("agent capability schema_version must be 1.1")
    version = _project_version(root)
    if manifest.get("generated_for_version") != version:
        failures.append("agent capability version must match pyproject.toml")
    if manifest.get("$schema") != "https://dnncha.github.io/dotmatch/agent-capabilities.schema.json":
        failures.append("agent capability $schema must use the public schema URL")

    project = manifest.get("project")
    if not isinstance(project, dict) or project.get("name") != "DotMatch" or project.get("package") != "dotmatch":
        failures.append("agent capability project identity must be DotMatch/dotmatch")
    if not isinstance(project, dict) or project.get("license") != "Apache-2.0":
        failures.append("agent capability project license must be Apache-2.0")

    scope = manifest.get("scope")
    if not isinstance(scope, dict) or "fixed-window known-target short-DNA assignment" not in str(scope.get("summary", "")):
        failures.append("agent capability scope must state fixed-window known-target short-DNA assignment")
    if not isinstance(scope, dict) or len(scope.get("does_not_do", [])) < 5:
        failures.append("agent capability scope must state at least five explicit non-goals")

    outcomes = manifest.get("outcomes")
    if not isinstance(outcomes, dict) or set(outcomes) != REQUIRED_OUTCOMES:
        failures.append("agent capability outcomes must be exactly unique, ambiguous, none, and invalid")

    intents = manifest.get("intents")
    if not isinstance(intents, list):
        failures.append("agent capability intents must be a list")
        intents = []
    ids = [intent.get("id") for intent in intents if isinstance(intent, dict)]
    if len(ids) != len(set(ids)):
        failures.append("agent capability intent ids must be unique")
    failures.extend(f"agent capability is missing required intent: {intent_id}" for intent_id in sorted(REQUIRED_INTENTS - set(ids)))

    required_intent_fields = {
        "id",
        "task",
        "queries",
        "entrypoint",
        "command",
        "inputs",
        "outputs",
        "limitations",
        "documentation",
        "evidence",
    }
    for index, intent in enumerate(intents):
        if not isinstance(intent, dict):
            failures.append(f"intent {index} must be an object")
            continue
        intent_id = str(intent.get("id") or index)
        missing = sorted(required_intent_fields - set(intent))
        failures.extend(f"intent {intent_id} missing field: {field}" for field in missing)
        if len(intent.get("queries", [])) < 2:
            failures.append(f"intent {intent_id} must provide at least two search queries")
        entrypoint = str(intent.get("entrypoint", ""))
        command = str(intent.get("command", ""))
        if not entrypoint.startswith("dotmatch ") or not command.startswith(entrypoint):
            failures.append(f"intent {intent_id} command must start with its DotMatch entrypoint")
        if not intent.get("inputs") or not intent.get("outputs") or not intent.get("limitations"):
            failures.append(f"intent {intent_id} must state inputs, outputs, and limitations")
        documentation = str(intent.get("documentation", ""))
        if not documentation.startswith("https://dotmatch.readthedocs.io/"):
            failures.append(f"intent {intent_id} documentation must use the rendered DotMatch docs")
        for evidence in intent.get("evidence", []):
            evidence_path = root / str(evidence)
            if not evidence_path.is_file():
                failures.append(f"intent {intent_id} evidence path does not exist: {evidence}")

    recovery = manifest.get("error_recovery")
    if not isinstance(recovery, list) or len(recovery) < 3:
        failures.append("agent capability error_recovery must provide at least three routes")
    else:
        for item in recovery:
            if not isinstance(item, dict) or not str(item.get("next_command", "")).startswith("dotmatch "):
                failures.append("every error recovery route must provide a DotMatch next_command")

    interfaces = manifest.get("interfaces")
    if not isinstance(interfaces, dict) or interfaces.get("machine_help") != "dotmatch capabilities --json":
        failures.append("agent capability machine_help must be dotmatch capabilities --json")
    for key in ("agent_tools", "agent_tools_schema", "codex_skill"):
        if not isinstance(interfaces, dict) or not str(interfaces.get(key, "")).startswith("https://"):
            failures.append(f"agent capability interfaces must link {key}")
    return failures


def validate_agent_tools(root: Path) -> list[str]:
    failures: list[str] = []
    try:
        contract = _json(root, "agent-tools.json")
        schema = _json(root, "agent-tools.schema.json")
    except Exception as exc:
        return [f"agent tool contract is invalid: {exc}"]
    if contract.get("tool_contract_version") != "1.0" or contract.get("schema_version") != "1.0":
        failures.append("agent tool contract versions must be 1.0")
    if contract.get("generated_for_version") != _project_version(root):
        failures.append("agent tool contract version must match pyproject.toml")
    names = [item.get("name") for item in contract.get("tools", []) if isinstance(item, dict)]
    required = {"discover", "prepare_assay", "inspect_assay", "run_assay", "review_assay", "handoff_assay"}
    if set(names) != required or len(names) != len(required):
        failures.append("agent tool contract must define exactly the six supported tools")
    if schema.get("$id") != "https://dnncha.github.io/dotmatch/agent-tools.schema.json":
        failures.append("agent tool schema must use the public $id")
    if not bool((contract.get("safety") or {}).get("ambiguity_is_never_forced_into_counts")):
        failures.append("agent tool contract must preserve ambiguous assignments")
    extension = root / "extensions/codex/dotmatch-agent"
    packaged = root / "python/dotmatch/data/codex-skill"
    if not extension.is_dir() or not packaged.is_dir():
        failures.append("Codex skill must exist in extension and installed-package locations")
    else:
        ext_files = {p.relative_to(extension).as_posix(): p.read_bytes() for p in extension.rglob("*") if p.is_file()}
        pkg_files = {p.relative_to(packaged).as_posix(): p.read_bytes() for p in packaged.rglob("*") if p.is_file()}
        if ext_files != pkg_files:
            failures.append("packaged Codex skill is stale")
    return failures


def validate_schema(root: Path) -> list[str]:
    failures: list[str] = []
    try:
        schema = _json(root, "agent-capabilities.schema.json")
    except Exception as exc:
        return [f"agent-capabilities.schema.json is invalid: {exc}"]
    if schema.get("$schema") != "https://json-schema.org/draft/2020-12/schema":
        failures.append("agent capability schema must declare JSON Schema draft 2020-12")
    if schema.get("$id") != "https://dnncha.github.io/dotmatch/agent-capabilities.schema.json":
        failures.append("agent capability schema $id must use the public schema URL")
    if schema.get("type") != "object" or not isinstance(schema.get("$defs"), dict):
        failures.append("agent capability schema must define an object and reusable definitions")
    return failures


def validate_legacy_schema(root: Path) -> list[str]:
    failures: list[str] = []
    try:
        manifest = _json(root, "agent-capabilities.v1.json")
        schema = _json(root, "agent-capabilities.v1.schema.json")
    except Exception as exc:
        return [f"legacy agent capability snapshot is invalid: {exc}"]
    if manifest.get("schema_version") != "1.0":
        failures.append("legacy agent capability snapshot must remain schema 1.0")
    if manifest.get("$schema") != "https://dnncha.github.io/dotmatch/agent-capabilities.v1.schema.json":
        failures.append("legacy agent capability snapshot must use its immutable schema URL")
    if schema.get("$id") != "https://dnncha.github.io/dotmatch/agent-capabilities.v1.schema.json":
        failures.append("legacy agent capability schema must use its immutable $id")
    if ((schema.get("properties") or {}).get("schema_version") or {}).get("const") != "1.0":
        failures.append("legacy agent capability schema must remain 1.0")
    return failures


def validate_copies(root: Path) -> list[str]:
    failures: list[str] = []
    for source_name, copy_names in PUBLISHED_COPIES.items():
        source = root / source_name
        if not source.is_file():
            failures.append(f"missing canonical agent discovery file: {source_name}")
            continue
        expected = source.read_bytes()
        for copy_name in copy_names:
            copy = root / copy_name
            if not copy.is_file():
                failures.append(f"missing published agent discovery copy: {copy_name}")
            elif copy.read_bytes() != expected:
                failures.append(f"published agent discovery copy is stale: {copy_name}")
    return failures


def validate_surfaces(root: Path) -> list[str]:
    failures: list[str] = []
    required_files = [
        "docs/agent-guide.md",
        "scripts/sync_agent_discovery.py",
        "scripts/check_agent_discovery.py",
    ]
    failures.extend(f"missing agent discovery surface: {path}" for path in required_files if not (root / path).is_file())
    if failures:
        return failures

    pyproject = _read(root, "pyproject.toml")
    if '"Agent guide" = "https://dotmatch.readthedocs.io/en/latest/agent-guide.html"' not in pyproject:
        failures.append("pyproject.toml must publish the Agent guide project URL")
    release = re.search(r'^version\s*=\s*"([^"]+)"', pyproject, re.M).group(1)
    readme = _read(root, "README.md")
    readme_surface = " ".join(readme.split())
    for phrase in ["Choose by task", "Agent guide", "dotmatch capabilities --json"]:
        if phrase not in readme:
            failures.append(f"README.md missing agent discovery phrase: {phrase}")
    for phrase in [
        f"Release {release}",
        "six `dotmatch agent` tools",
        "includes the six",
    ]:
        if phrase not in readme_surface:
            failures.append(f"README.md must state the public/agent-tools version boundary: {phrase}")
    docs_index = _read(root, "docs/index.md")
    docs_index_surface = " ".join(docs_index.split())
    if "agent-guide" not in docs_index:
        failures.append("docs/index.md must route to the agent guide")
    for phrase in [
        f"Release {release}",
        "six `dotmatch agent` tools",
        "includes the six",
    ]:
        if phrase not in docs_index_surface:
            failures.append(f"docs/index.md must state the public/agent-tools version boundary: {phrase}")
    agent_guide = _read(root, "docs/agent-guide.md")
    agent_guide_surface = " ".join(agent_guide.split())
    for phrase in [
        f"This six-tool interface is included in release {release}.",
        "older package",
        f"until that channel reaches {release}",
    ]:
        if phrase not in agent_guide_surface:
            failures.append(f"docs/agent-guide.md must state the public/agent-tools version boundary: {phrase}")
    docs_conf = _read(root, "docs/conf.py")
    for name in ["llms.txt", "llms-full.txt", "agent-capabilities.json", "agent-capabilities.schema.json", "agent-tools.json", "agent-tools.schema.json"]:
        if name not in docs_conf:
            failures.append(f"docs/conf.py html_extra_path must publish {name}")
    help_source = _read(root, "python/dotmatch/cli.py")
    for phrase in ["Choose by task:", "dotmatch capabilities --json", "Feature-barcode assignment", "Perturb-seq guide capture"]:
        if phrase not in help_source:
            failures.append(f"CLI help source missing agent route: {phrase}")
    layout = _read(root, "app/layout.tsx")
    if 'rel="describedby"' not in layout or "llms.txt" not in layout:
        failures.append("homepage metadata must point agents to llms.txt with rel=describedby")
    page = _read(root, "app/page.tsx")
    if "agent-capabilities.json" not in page or "featureList" not in page:
        failures.append("homepage structured data must expose the capability manifest and feature list")
    manifest = _read(root, "MANIFEST.in")
    for name in ["agent-capabilities.json", "agent-capabilities.schema.json", "agent-tools.json", "agent-tools.schema.json"]:
        if f"include {name}" not in manifest:
            failures.append(f"MANIFEST.in must include {name}")
    wheel_checker = _read(root, "scripts/check_python_wheel.py")
    for phrase in ["agent-capabilities.json", '"capabilities"', '"--json"', "agent_smoke_summary.json"]:
        if phrase not in wheel_checker:
            failures.append(f"clean-install package gate missing proof hook: {phrase}")
    return failures


def local_measure(root: Path) -> dict[str, Any]:
    checks: list[tuple[str, bool, str]] = []

    def add(check_id: str, condition: bool, evidence: str) -> None:
        checks.append((check_id, bool(condition), evidence))

    manifest_path = root / "agent-capabilities.json"
    manifest_valid = manifest_path.is_file() and not validate_manifest(root)
    copies_valid = not validate_copies(root) if manifest_path.is_file() else False
    readme = _read(root, "README.md") if (root / "README.md").is_file() else ""
    cli = _read(root, "python/dotmatch/cli.py") if (root / "python/dotmatch/cli.py").is_file() else ""
    docs_conf = _read(root, "docs/conf.py") if (root / "docs/conf.py").is_file() else ""
    layout = _read(root, "app/layout.tsx") if (root / "app/layout.tsx").is_file() else ""
    page = _read(root, "app/page.tsx") if (root / "app/page.tsx").is_file() else ""
    pyproject = _read(root, "pyproject.toml") if (root / "pyproject.toml").is_file() else ""
    wheel_checker = _read(root, "scripts/check_python_wheel.py") if (root / "scripts/check_python_wheel.py").is_file() else ""

    add("readme_task_routing", "Choose by task" in readme and "Agent guide" in readme, "README task table and agent guide link")
    add("cli_task_routing", "Choose by task:" in cli and "Perturb-seq guide capture" in cli, "installed help task routes")
    add("machine_cli", "dotmatch capabilities --json" in cli, "installed JSON capability command")
    add("capability_manifest", manifest_path.is_file(), "canonical agent-capabilities.json")
    add("capability_manifest_valid", manifest_valid, "required intents, inputs, outputs, limitations, and evidence")
    add("capability_schema", (root / "agent-capabilities.schema.json").is_file() and not validate_schema(root), "JSON Schema draft 2020-12")
    add("llms_index", (root / "llms.txt").is_file(), "concise llms.txt")
    add("llms_full", (root / "llms-full.txt").is_file(), "self-contained llms-full.txt")
    add("published_copies", copies_valid, "GitHub Pages, Read the Docs, and installed-package copies aligned")
    add("structured_web_discovery", 'rel="describedby"' in layout and "featureList" in page, "describedby link and SoftwareApplication feature list")
    add("pypi_agent_route", "Agent guide" in pyproject, "PyPI project URL metadata")
    add("clean_install_workflow", "agent_smoke_summary.json" in wheel_checker, "fresh-venv wheel and sdist FASTQ count gate")

    passed = sum(1 for _check_id, ok, _evidence in checks if ok)
    return {
        "mode": "local",
        "root": str(root),
        "score": passed,
        "maximum": len(checks),
        "checks": [
            {"id": check_id, "passed": ok, "evidence": evidence}
            for check_id, ok, evidence in checks
        ],
    }


def _fetch(url: str) -> tuple[int, str, str]:
    request = urllib.request.Request(url, headers={"User-Agent": "dotmatch-agent-discovery-audit/1"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return int(response.status), str(response.headers.get("Content-Type", "")), response.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        return int(exc.code), str(exc.headers.get("Content-Type", "")), exc.read().decode("utf-8", errors="replace")
    except Exception as exc:
        return 0, "", str(exc)


def live_measure() -> dict[str, Any]:
    checks: list[dict[str, Any]] = []

    def add(check_id: str, passed: bool, evidence: str) -> None:
        checks.append({"id": check_id, "passed": bool(passed), "evidence": evidence})

    gh_status, _gh_type, gh_text = _fetch("https://api.github.com/repos/dnncha/dotmatch")
    try:
        gh = json.loads(gh_text)
    except Exception:
        gh = {}
    topics = set(gh.get("topics", [])) if isinstance(gh, dict) else set()
    add(
        "github_search_metadata",
        gh_status == 200 and gh.get("visibility") == "public" and REQUIRED_SEARCH_TERMS <= topics,
        f"HTTP {gh_status}; required topics {sorted(REQUIRED_SEARCH_TERMS)}",
    )

    pypi_status, _pypi_type, pypi_text = _fetch("https://pypi.org/pypi/dotmatch/json")
    try:
        pypi = json.loads(pypi_text).get("info", {})
    except Exception:
        pypi = {}
    pypi_blob = " ".join([str(pypi.get("summary", "")), str(pypi.get("keywords", "")), json.dumps(pypi.get("project_urls", {}))]).lower()
    add(
        "pypi_scope_metadata",
        pypi_status == 200 and all(term in pypi_blob for term in ["known-target", "crispr", "fastq", "barcode", "documentation"]),
        f"HTTP {pypi_status}; version {pypi.get('version', 'unknown')}",
    )

    rtd_status, _rtd_type, rtd_text = _fetch("https://dotmatch.readthedocs.io/en/latest/")
    add(
        "readthedocs_onboarding",
        rtd_status == 200 and all(term in rtd_text for term in ["Getting started", "Choose a workflow", "Scope and limitations"]),
        f"HTTP {rtd_status}; rendered onboarding terms",
    )

    pages_status, _pages_type, pages_text = _fetch("https://dnncha.github.io/dotmatch/")
    add(
        "pages_structured_scope",
        pages_status == 200 and "application/ld+json" in pages_text and all(term in pages_text for term in ["CRISPR guides", "inline barcodes", "feature tags"]),
        f"HTTP {pages_status}; JSON-LD and task vocabulary",
    )

    pages_llms_status, pages_llms_type, _pages_llms_text = _fetch("https://dnncha.github.io/dotmatch/llms.txt")
    add("pages_llms", pages_llms_status == 200 and "text/plain" in pages_llms_type, f"HTTP {pages_llms_status}; {pages_llms_type or 'no content type'}")

    rtd_llms_status, rtd_llms_type, _rtd_llms_text = _fetch("https://dotmatch.readthedocs.io/en/latest/llms.txt")
    add("readthedocs_llms", rtd_llms_status == 200 and "text/plain" in rtd_llms_type, f"HTTP {rtd_llms_status}; {rtd_llms_type or 'no content type'}")

    cap_status, cap_type, cap_text = _fetch("https://dnncha.github.io/dotmatch/agent-capabilities.json")
    try:
        cap = json.loads(cap_text)
    except Exception:
        cap = {}
    cap_ids = {item.get("id") for item in cap.get("intents", []) if isinstance(item, dict)} if isinstance(cap, dict) else set()
    add(
        "public_capability_manifest",
        cap_status == 200 and "json" in cap_type and REQUIRED_INTENTS <= cap_ids,
        f"HTTP {cap_status}; {cap_type or 'no content type'}; required intent ids",
    )

    passed = sum(1 for check in checks if check["passed"])
    return {"mode": "live", "score": passed, "maximum": len(checks), "checks": checks}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT, help="repository root to inspect")
    parser.add_argument("--measure", action="store_true", help="print a deterministic local score without enforcing the final gate")
    parser.add_argument("--live", action="store_true", help="measure the currently deployed public surfaces")
    parser.add_argument("--json", action="store_true", help="emit JSON")
    args = parser.parse_args()

    if args.live:
        report = live_measure()
        print(json.dumps(report, indent=2) if args.json else f"Live agent discoverability: {report['score']}/{report['maximum']}")
        return 0

    root = args.root.resolve()
    if args.measure:
        report = local_measure(root)
        print(json.dumps(report, indent=2) if args.json else f"Local agent discoverability: {report['score']}/{report['maximum']}")
        return 0

    failures = validate_schema(root) + validate_legacy_schema(root) + validate_manifest(root) + validate_agent_tools(root) + validate_copies(root) + validate_surfaces(root)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}", file=sys.stderr)
        return 1
    report = local_measure(root)
    print(f"Agent discovery gate passed ({report['score']}/{report['maximum']} checks)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
