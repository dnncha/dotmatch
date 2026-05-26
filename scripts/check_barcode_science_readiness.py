#!/usr/bin/env python3
"""Audit the validation bundle behind DotMatch's fixed-window barcode work."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Iterable


REQUIRED_MIN_PUBLIC_DATASETS = 5
EXACT_BASELINE_TOOLS = {"exact_slice_hash", "exact_prefix_hash", "hash_splitter_exact"}


class AuditResult:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _repo_path(root: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"path must be repository-relative: {value}")
    return root / path


def _read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        return list(reader.fieldnames or []), list(reader)


def _make_targets(root: Path) -> set[str]:
    makefile = root / "Makefile"
    if not makefile.exists():
        return set()
    targets: set[str] = set()
    for line in makefile.read_text(encoding="utf-8").splitlines():
        match = re.match(r"^([A-Za-z0-9_.-]+):(?:\s|$)", line)
        if match:
            targets.add(match.group(1))
    return targets


def _as_int(value: object, default: int = 0) -> int:
    text = str(value or "").strip()
    if not text:
        return default
    try:
        return int(float(text))
    except ValueError:
        return default


def _require(condition: bool, message: str, result: AuditResult) -> None:
    if not condition:
        result.failures.append(message)


def _public_rows(rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    public: list[dict[str, str]] = []
    for row in rows:
        workflow = str(row.get("workflow") or "").lower()
        metadata = str(row.get("metadata") or "").lower()
        status = str(row.get("status") or "").lower()
        if "synthetic" in workflow or "fixture" in workflow or status == "smoke":
            continue
        if "public" in workflow or "real" in workflow or metadata:
            public.append(row)
    return public


def _has_successful_dotmatch(rows: list[dict[str, str]]) -> bool:
    return any(
        str(row.get("tool") or "").startswith("dotmatch")
        and str(row.get("exit_code") or "") == "0"
        and _assigned_reads(row) > 0
        for row in rows
    )


def _assigned_reads(row: dict[str, str]) -> int:
    for key in ["assigned_unique", "assigned_reads", "assigned_pairs"]:
        value = _as_int(row.get(key), default=-1)
        if value >= 0:
            return value
    return 0


def _has_exact_baseline(rows: list[dict[str, str]]) -> bool:
    return any(
        str(row.get("tool") or "") in EXACT_BASELINE_TOOLS
        and str(row.get("exit_code") or "") == "0"
        and _assigned_reads(row) > 0
        for row in rows
    )


def _has_inline_barcode_comparators(rows: list[dict[str, str]]) -> bool:
    tools = {str(row.get("tool") or "") for row in rows if str(row.get("exit_code") or "") == "0"}
    return "cutadapt_demux" in tools and "hash_splitter_exact" in tools


def _check_zero_validation_mismatches(dataset_id: str, raw_artifact: str, rows: list[dict[str, str]], result: AuditResult) -> None:
    for index, row in enumerate(rows, start=2):
        raw = str(row.get("validation_mismatches") or "").strip()
        if not raw:
            continue
        if _as_int(raw, default=-1) != 0:
            result.failures.append(f"{dataset_id}: {raw_artifact}:{index} must have zero validation_mismatches")


def _check_dataset(root: Path, dataset: dict, make_targets: set[str], result: AuditResult) -> None:
    dataset_id = str(dataset.get("id") or "")
    if not dataset_id:
        result.failures.append("barcode validation dataset is missing id")
        return
    for key in ["label", "raw_artifact", "metadata", "gate", "comparator_semantics", "claim_boundary"]:
        if not str(dataset.get(key) or "").strip():
            result.failures.append(f"{dataset_id}: missing {key}")

    try:
        metadata_path = _repo_path(root, str(dataset.get("metadata") or ""))
        raw_path = _repo_path(root, str(dataset.get("raw_artifact") or ""))
    except ValueError as exc:
        result.failures.append(f"{dataset_id}: {exc}")
        return

    _require(metadata_path.is_file(), f"{dataset_id}: missing metadata {dataset.get('metadata')}", result)
    if metadata_path.is_file():
        metadata = _load_json(metadata_path)
        _require(bool(metadata.get("evidence_ready")), f"{dataset_id}: metadata is not evidence_ready", result)
    _require(raw_path.is_file(), f"{dataset_id}: missing raw artifact {dataset.get('raw_artifact')}", result)
    if not raw_path.is_file():
        return

    fieldnames, rows = _read_rows(raw_path)
    _require(bool(fieldnames), f"{dataset_id}: raw artifact must contain a CSV header", result)
    public_rows = _public_rows(rows)
    _require(bool(public_rows), f"{dataset_id}: raw artifact must include public or real-data rows", result)
    _require(_has_successful_dotmatch(public_rows), f"{dataset_id}: public rows must include successful DotMatch assignments", result)
    _check_zero_validation_mismatches(dataset_id, str(dataset.get("raw_artifact")), public_rows, result)

    if dataset_id == "inline_barcode_srp009896":
        _require(
            _has_inline_barcode_comparators(public_rows),
            f"{dataset_id}: inline barcode evidence must include Cutadapt and hash-splitter comparator rows",
            result,
        )
    else:
        _require(
            _has_exact_baseline(public_rows),
            f"{dataset_id}: public evidence must include a transparent exact baseline comparator",
            result,
        )

    gate = str(dataset.get("gate") or "")
    parts = gate.split()
    _require(len(parts) == 2 and parts[0] == "make", f"{dataset_id}: gate must be a simple make target", result)
    if len(parts) == 2 and parts[0] == "make":
        _require(parts[1] in make_targets, f"{dataset_id}: missing make target for gate {parts[1]}", result)


def audit(root: Path) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    manifest_path = root / "docs" / "barcode-science-readiness.json"
    if not manifest_path.is_file():
        result.failures.append("missing docs/barcode-science-readiness.json")
        return result
    try:
        manifest = _load_json(manifest_path)
    except Exception as exc:
        result.failures.append(f"docs/barcode-science-readiness.json could not be read: {exc}")
        return result

    _require(manifest.get("schema_version") == 1, "barcode validation manifest must declare schema_version 1", result)
    datasets = manifest.get("datasets")
    if not isinstance(datasets, list):
        result.failures.append("barcode validation manifest must contain datasets list")
        return result
    _require(
        len(datasets) >= REQUIRED_MIN_PUBLIC_DATASETS,
        f"barcode validation requires at least {REQUIRED_MIN_PUBLIC_DATASETS} public fixed-window datasets",
        result,
    )

    make_targets = _make_targets(root)
    seen: set[str] = set()
    for dataset in datasets:
        if not isinstance(dataset, dict):
            result.failures.append("barcode validation datasets must be objects")
            continue
        dataset_id = str(dataset.get("id") or "")
        if dataset_id in seen:
            result.failures.append(f"duplicate barcode validation dataset: {dataset_id}")
        seen.add(dataset_id)
        _check_dataset(root, dataset, make_targets, result)

    if result.ok:
        result.passed.append(f"{len(datasets)} public fixed-window datasets are comparator-backed")
        result.passed.append("barcode validation manifest references raw artifacts and gates")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="repository root")
    args = parser.parse_args()
    result = audit(Path(args.root))
    for item in result.passed:
        print(f"PASS: {item}")
    for failure in result.failures:
        print(f"FAIL: {failure}")
    if result.ok:
        print("BARCODE VALIDATION: PASS")
        return 0
    print("BARCODE VALIDATION: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
