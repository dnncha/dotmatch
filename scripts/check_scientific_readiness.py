#!/usr/bin/env python3

import argparse
import json
import re
from pathlib import Path


REQUIRED_CONTROLS = {
    "input_integrity",
    "memory_safety",
    "oracle_validation",
    "public_assay_evidence",
    "distribution_reproducibility",
    "release_governance",
}


class AuditResult:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def _make_targets(root: Path) -> set[str]:
    text = (root / "Makefile").read_text(encoding="utf-8")
    targets: set[str] = set()
    for line in text.splitlines():
        match = re.match(r"^([A-Za-z0-9_.-]+):(?:\s|$)", line)
        if match:
            targets.add(match.group(1))
    return targets


def _check_path(root: Path, field: str, value: str, result: AuditResult) -> None:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        result.failures.append(f"{field} must be repository-relative: {value}")
        return
    if not (root / path).exists():
        result.failures.append(f"missing {field}: {value}")


def _check_gate(gate: str, make_targets: set[str], result: AuditResult) -> None:
    parts = gate.split()
    if len(parts) != 2 or parts[0] != "make":
        result.failures.append(f"gate must be a simple make target command: {gate}")
        return
    if parts[1] not in make_targets:
        result.failures.append(f"missing make target for scientific readiness gate: {parts[1]}")


def audit(root: Path) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    path = root / "docs" / "scientific-readiness.json"
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        result.failures.append(f"docs/scientific-readiness.json could not be read: {exc}")
        return result

    if manifest.get("schema_version") != 1:
        result.failures.append("docs/scientific-readiness.json must declare schema_version 1")
    if manifest.get("status") != "evidence_bounded":
        result.failures.append("scientific readiness status must be evidence_bounded")
    if not manifest.get("scope"):
        result.failures.append("scientific readiness manifest must declare scope")
    not_validated = manifest.get("not_validated_for")
    if not isinstance(not_validated, list) or not not_validated:
        result.failures.append("scientific readiness manifest must declare not_validated_for boundaries")

    controls = manifest.get("controls")
    if not isinstance(controls, list):
        result.failures.append("scientific readiness manifest must contain a controls list")
        return result

    make_targets = _make_targets(root)
    seen: set[str] = set()
    for control in controls:
        if not isinstance(control, dict):
            result.failures.append("scientific readiness controls must be objects")
            continue
        control_id = str(control.get("id") or "")
        if not control_id:
            result.failures.append("scientific readiness control missing id")
            continue
        if control_id in seen:
            result.failures.append(f"duplicate scientific readiness control: {control_id}")
        seen.add(control_id)
        if control.get("status") != "required":
            result.failures.append(f"{control_id} status must be required")
        if not control.get("acceptance"):
            result.failures.append(f"{control_id} must declare acceptance")
        evidence = control.get("evidence")
        if not isinstance(evidence, list) or not evidence:
            result.failures.append(f"{control_id} must list evidence")
        else:
            for value in evidence:
                _check_path(root, f"{control_id} evidence", str(value), result)
        gates = control.get("gates")
        if not isinstance(gates, list) or not gates:
            result.failures.append(f"{control_id} must list gates")
        else:
            for gate in gates:
                _check_gate(str(gate), make_targets, result)

    missing = REQUIRED_CONTROLS - seen
    for control_id in sorted(missing):
        result.failures.append(f"missing required scientific readiness control: {control_id}")
    if not missing and not result.failures:
        result.passed.append("scientific readiness controls valid")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit DotMatch scientific readiness controls.")
    parser.add_argument("--root", default=".", help="repository root")
    args = parser.parse_args()

    result = audit(Path(args.root))
    for item in result.passed:
        print(f"PASS: {item}")
    for item in result.failures:
        print(f"FAIL: {item}")
    if result.ok:
        print("SCIENTIFIC READINESS: PASS")
        return 0
    print("SCIENTIFIC READINESS: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
