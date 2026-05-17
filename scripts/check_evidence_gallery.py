#!/usr/bin/env python3
"""Validate the manifest-driven DotMatch evidence gallery."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import generate_evidence_gallery  # noqa: E402


MANIFEST = Path("docs") / "evidence-gallery" / "manifest.json"
VALID_STATUSES = {"supported", "gated", "smoke"}
REQUIRED_SCENARIOS = {
    "public_crispr_yusa",
    "barcode_autopsy_review",
    "barcode_wrong_offset_fixture",
    "barcode_unsafe_correction",
    "feature_barcode_10x",
    "perturb_seq_10x_guide_capture",
    "amplicon_artic_primer_start",
    "oligo_adapter_truseq_prefix",
    "bcl_tiny_classic",
}
REQUIRED_ROLES = {"known_good", "low_confidence", "wrong_offset", "unsafe_correction", "gated_parser_milestone"}


class AuditResult:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def _load_manifest(root: Path, result: AuditResult) -> dict[str, Any]:
    try:
        return json.loads((root / MANIFEST).read_text(encoding="utf-8"))
    except Exception as exc:
        result.failures.append(f"{MANIFEST.as_posix()} could not be read: {exc}")
        return {}


def _make_targets(root: Path) -> set[str]:
    path = root / "Makefile"
    if not path.is_file():
        return set()
    targets: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"^([A-Za-z0-9_.-]+):(?:\s|$)", line)
        if match:
            targets.add(match.group(1))
    return targets


def _as_list(value: object) -> list[Any]:
    return value if isinstance(value, list) else []


def _check_repo_file(root: Path, field: str, value: object, result: AuditResult) -> None:
    text = str(value or "")
    path = Path(text)
    if not text or path.is_absolute() or ".." in path.parts:
        result.failures.append(f"{field} must be a repository-relative file path: {text}")
        return
    if not (root / path).is_file():
        result.failures.append(f"missing {field}: {text}")


def _check_command(scenario_id: str, command: object, make_targets: set[str], result: AuditResult) -> None:
    text = str(command or "").strip()
    if not text:
        result.failures.append(f"{scenario_id} command must not be empty")
        return
    parts = text.split()
    if parts[0] == "make" and len(parts) >= 2 and parts[1] not in make_targets:
        result.failures.append(f"{scenario_id} command references missing make target: {parts[1]}")


def _check_scenario(root: Path, scenario: dict[str, Any], make_targets: set[str], result: AuditResult) -> None:
    scenario_id = str(scenario.get("id") or "")
    if not scenario_id:
        result.failures.append("evidence gallery scenario missing id")
        return
    for field in ["title", "category", "assay_type", "dataset", "condition", "comparator_semantics", "validation"]:
        if not str(scenario.get(field) or "").strip():
            result.failures.append(f"{scenario_id} missing {field}")
    if scenario.get("status") not in VALID_STATUSES:
        result.failures.append(f"{scenario_id} has invalid status: {scenario.get('status')}")
    if not _as_list(scenario.get("gallery_roles")):
        result.failures.append(f"{scenario_id} must list gallery_roles")
    for field in ["proves", "limits", "commands", "raw_artifacts", "report_examples"]:
        if not _as_list(scenario.get(field)):
            result.failures.append(f"{scenario_id} must list {field}")

    _check_repo_file(root, "primary_report", scenario.get("primary_report"), result)
    for raw in _as_list(scenario.get("raw_artifacts")):
        _check_repo_file(root, "raw_artifact", raw, result)
    for command in _as_list(scenario.get("commands")):
        _check_command(scenario_id, command, make_targets, result)
    for example in _as_list(scenario.get("report_examples")):
        if not isinstance(example, dict):
            result.failures.append(f"{scenario_id} report_examples entries must be objects")
            continue
        if not str(example.get("label") or "").strip():
            result.failures.append(f"{scenario_id} report example missing label")
        _check_repo_file(root, "report_example", example.get("path"), result)


def _check_generated_files(root: Path, result: AuditResult) -> None:
    try:
        expected = generate_evidence_gallery.render_files(root)
    except Exception as exc:
        result.failures.append(f"could not render evidence gallery: {exc}")
        return
    for rel_path, text in expected.items():
        path = root / rel_path
        if not path.is_file():
            result.failures.append(f"missing generated evidence gallery file: {rel_path}")
            continue
        current = path.read_text(encoding="utf-8")
        if current != text.rstrip() + "\n":
            result.failures.append(f"generated evidence gallery file is stale: {rel_path}")


def audit(root: Path) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    manifest = _load_manifest(root, result)
    if not manifest:
        return result
    if manifest.get("schema_version") != 1:
        result.failures.append("evidence gallery manifest must declare schema_version 1")
    scenarios = manifest.get("scenarios")
    if not isinstance(scenarios, list) or not all(isinstance(item, dict) for item in scenarios):
        result.failures.append("evidence gallery manifest must contain a scenarios object list")
        return result

    seen_ids: set[str] = set()
    seen_roles: set[str] = set()
    make_targets = _make_targets(root)
    for scenario in scenarios:
        scenario_id = str(scenario.get("id") or "")
        if scenario_id in seen_ids:
            result.failures.append(f"duplicate evidence gallery scenario: {scenario_id}")
        seen_ids.add(scenario_id)
        seen_roles.update(str(role) for role in _as_list(scenario.get("gallery_roles")))
        _check_scenario(root, scenario, make_targets, result)

    missing_scenarios = sorted(REQUIRED_SCENARIOS - seen_ids)
    if missing_scenarios:
        result.failures.append("missing required evidence gallery scenarios: " + ", ".join(missing_scenarios))
    missing_roles = sorted(REQUIRED_ROLES - seen_roles)
    if missing_roles:
        result.failures.append("missing required evidence gallery roles: " + ", ".join(missing_roles))

    _check_generated_files(root, result)
    if result.ok:
        result.passed.append("evidence gallery manifest and generated pages are ready")
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
        print("EVIDENCE GALLERY: PASS")
        return 0
    print("EVIDENCE GALLERY: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
