#!/usr/bin/env python3

import argparse
import csv
import json
import re
from pathlib import Path


REQUIRED_ASSAYS = [
    "crispr_guide_counting",
    "inline_barcode",
    "perturb_seq",
    "feature_barcode",
    "amplicon_panel",
    "oligo_adapter",
]
VALID_STATUSES = {"supported", "gated", "smoke", "planned"}


class AuditResult:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def _load_manifest(root: Path) -> dict:
    return json.loads((root / "docs" / "assay-evidence.json").read_text(encoding="utf-8"))


def _make_targets(root: Path) -> set[str]:
    text = (root / "Makefile").read_text(encoding="utf-8")
    targets: set[str] = set()
    for line in text.splitlines():
        match = re.match(r"^([A-Za-z0-9_.-]+):(?:\s|$)", line)
        if match:
            targets.add(match.group(1))
    return targets


def _check_relative_path(root: Path, field: str, value: str, result: AuditResult) -> None:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        result.failures.append(f"{field} must be a repository-relative path: {value}")
        return
    if not (root / path).is_file():
        kind = "raw artifact" if field == "raw_artifacts" else "report"
        result.failures.append(f"missing {kind}: {value}")


def _parse_int(value: str):
    text = value.strip()
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _check_raw_csv(root: Path, value: str, result: AuditResult) -> None:
    path = root / value
    if path.suffix.lower() != ".csv" or not path.is_file():
        return
    try:
        with path.open(encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh)
            rows = list(reader)
    except Exception as exc:
        result.failures.append(f"raw artifact {value} could not be read as CSV: {exc}")
        return

    if not reader.fieldnames:
        result.failures.append(f"raw artifact {value} must contain a CSV header")
        return
    if not rows:
        result.failures.append(f"raw artifact {value} must contain at least one data row")
        return

    if "command" in reader.fieldnames:
        for row_number, row in enumerate(rows, start=2):
            if not str(row.get("command") or "").strip():
                result.failures.append(f"raw artifact {value}:{row_number} must record command provenance")
            if "exit_code" in reader.fieldnames and not str(row.get("exit_code") or "").strip():
                result.failures.append(f"raw artifact {value}:{row_number} must record exit_code provenance")

    for column in ["validation_mismatches", "mismatches"]:
        if column not in reader.fieldnames:
            continue
        for row_number, row in enumerate(rows, start=2):
            raw = str(row.get(column) or "").strip()
            if not raw:
                continue
            value_int = _parse_int(raw)
            if value_int is None:
                result.failures.append(f"raw artifact {value}:{row_number} has nonnumeric {column}: {raw}")
            elif value_int != 0:
                result.failures.append(f"raw artifact {value}:{row_number} must have zero {column}")

    if "checked_reads" in reader.fieldnames:
        for row_number, row in enumerate(rows, start=2):
            raw = str(row.get("checked_reads") or "").strip()
            if not raw:
                continue
            checked = _parse_int(raw)
            if checked is None or checked <= 0:
                result.failures.append(f"raw artifact {value}:{row_number} must record positive checked_reads")


def _check_gate(gate: str, make_targets: set[str], result: AuditResult) -> None:
    parts = gate.split()
    if len(parts) != 2 or parts[0] != "make":
        result.failures.append(f"gate must be a simple make target command: {gate}")
        return
    target = parts[1]
    if target not in make_targets:
        result.failures.append(f"missing make target for assay evidence gate: {target}")


def _check_command(assay_id: str, command: str, make_targets: set[str], result: AuditResult) -> None:
    if not command.strip():
        result.failures.append(f"{assay_id} commands must not contain empty entries")
        return
    parts = command.split()
    if parts and parts[0] == "make" and len(parts) >= 2:
        target = parts[1]
        if target not in make_targets:
            result.failures.append(f"{assay_id} command references missing make target: {target}")


def _assay_id(assay: object) -> str:
    if isinstance(assay, dict):
        return str(assay.get("id") or "")
    return ""


def check_required_assays(assays: list[dict], result: AuditResult) -> None:
    present = {_assay_id(assay) for assay in assays}
    for assay_id in REQUIRED_ASSAYS:
        if assay_id not in present:
            result.failures.append(f"missing required assay lane: {assay_id}")
    if all(assay_id in present for assay_id in REQUIRED_ASSAYS):
        result.passed.append("required assay lanes present")


def check_assay_entries(root: Path, assays: list[dict], result: AuditResult) -> None:
    make_targets = _make_targets(root)
    seen: set[str] = set()
    for assay in assays:
        assay_id = _assay_id(assay)
        if not assay_id:
            result.failures.append("assay entry missing id")
            continue
        if assay_id in seen:
            result.failures.append(f"duplicate assay lane: {assay_id}")
        seen.add(assay_id)

        status = str(assay.get("status") or "")
        if status not in VALID_STATUSES:
            result.failures.append(f"{assay_id} has invalid status: {status}")

        if not assay.get("label"):
            result.failures.append(f"{assay_id} must declare label")
        if not assay.get("claim_boundary"):
            result.failures.append(f"{assay_id} must declare claim_boundary")

        raw_artifacts = assay.get("raw_artifacts") or []
        reports = assay.get("reports") or []
        gates = assay.get("gates") or []
        commands = assay.get("commands") or []

        if status == "supported":
            if not raw_artifacts:
                result.failures.append(f"{assay_id} supported lane must list raw_artifacts")
            if not reports:
                result.failures.append(f"{assay_id} supported lane must list reports")
            if not gates:
                result.failures.append(f"{assay_id} supported lane must list gates")

        if status in {"supported", "gated", "smoke"}:
            if not commands:
                result.failures.append(f"{assay_id} must list exact commands for non-planned evidence lanes")
            if not assay.get("comparator_semantics"):
                result.failures.append(f"{assay_id} must declare comparator_semantics for non-planned evidence lanes")
            if not assay.get("validation"):
                result.failures.append(f"{assay_id} must declare validation for non-planned evidence lanes")

        if status in {"gated", "smoke", "planned"} and not assay.get("next_public_evidence"):
            result.failures.append(f"{assay_id} must declare next_public_evidence until public evidence is complete")

        for command in commands:
            _check_command(assay_id, str(command), make_targets, result)
        for value in raw_artifacts:
            _check_relative_path(root, "raw_artifacts", str(value), result)
            _check_raw_csv(root, str(value), result)
        for value in reports:
            _check_relative_path(root, "reports", str(value), result)
        for gate in gates:
            _check_gate(str(gate), make_targets, result)

    if not any(
        marker in failure
        for failure in result.failures
        for marker in [
            "assay",
            "commands",
            "comparator_semantics",
            "validation",
            "raw artifact",
            "report",
            "make target",
            "next_public_evidence",
        ]
    ):
        result.passed.append("assay evidence entries valid")


def audit(root: Path) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    try:
        manifest = _load_manifest(root)
    except Exception as exc:
        result.failures.append(f"docs/assay-evidence.json could not be read: {exc}")
        return result

    if manifest.get("schema_version") != 1:
        result.failures.append("docs/assay-evidence.json must declare schema_version 1")
    assays = manifest.get("assays")
    if not isinstance(assays, list):
        result.failures.append("docs/assay-evidence.json must contain an assays list")
        return result
    if not all(isinstance(assay, dict) for assay in assays):
        result.failures.append("docs/assay-evidence.json assays must be objects")
        return result

    check_required_assays(assays, result)
    check_assay_entries(root, assays, result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit DotMatch assay-evidence coverage and claim boundaries.")
    parser.add_argument("--root", default=".", help="repository root")
    args = parser.parse_args()

    result = audit(Path(args.root))
    for item in result.passed:
        print(f"PASS: {item}")
    for item in result.failures:
        print(f"FAIL: {item}")
    if result.ok:
        print("ASSAY EVIDENCE: PASS")
        return 0
    print("ASSAY EVIDENCE: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
