#!/usr/bin/env python3
"""Check the frozen CRISPR and Perturb-seq records used by agent task pages."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

try:
    from jsonschema import Draft202012Validator
except ModuleNotFoundError:
    Draft202012Validator = None  # type: ignore[assignment,misc]


ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    failures: list[str] = []
    fixture = ROOT / "examples" / "workflows" / "fixtures"
    expected_hashes = {
        "crispr_assay.toml": "871a58c33550d1b31e71eb03605bc7e9d2c4b3b1b3250982b5153054ea10bc9d",
        "expected_counts.mageck.tsv": "82c61e9c2eb2d7701fd7505f4fa466e50fdbabbaebc16c266f304e06f4c424d9",
        "expected_sample_qc.tsv": "2d8d905ee296adfc2a3711f5ae50dd18782468c7ac704b07dee015954f77651d",
    }
    for name, expected in expected_hashes.items():
        path = fixture / name
        if not path.is_file() or sha256(path) != expected:
            failures.append(f"CRISPR workflow fixture hash differs: {name}")

    reference = json.loads((ROOT / "agent-reference-crispr.json").read_text(encoding="utf-8"))
    envelope_contract = json.loads((ROOT / "agent-tools.schema.json").read_text(encoding="utf-8"))
    if Draft202012Validator is not None:
        envelope_schema = {"$ref": "#/$defs/envelope", "$defs": envelope_contract["$defs"]}
        try:
            Draft202012Validator(envelope_schema).validate(reference)
        except Exception as exc:
            failures.append(f"CRISPR fixture envelope does not validate: {exc}")
    elif not set(envelope_contract["$defs"]["envelope"]["required"]) <= set(reference):
        failures.append("CRISPR fixture envelope is missing required contract fields")
    if reference.get("status") != "failed" or reference.get("exit_code") != 2:
        failures.append("CRISPR fixture envelope must retain its failed/2 reliability verdict")
    if reference.get("spec", {}).get("sha256") != expected_hashes["crispr_assay.toml"]:
        failures.append("CRISPR fixture envelope spec hash differs")
    reference_artifacts = {item.get("role"): item for item in reference.get("artifacts", [])}
    for role, name in {
        "counts": "expected_counts.mageck.tsv",
        "sample_qc": "expected_sample_qc.tsv",
    }.items():
        item = reference_artifacts.get(role, {})
        path = ROOT / str(item.get("path", ""))
        if not path.is_file() or item.get("sha256") != sha256(path) or item.get("bytes") != path.stat().st_size:
            failures.append(f"CRISPR fixture envelope {role} artifact is stale")
        if item.get("safe_to_share") is not True:
            failures.append(f"public synthetic CRISPR fixture {role} artifact must be marked safe to share")
    outcomes = set(reference.get("result", {}).get("outcomes_preserved", []))
    if outcomes != {"unique", "ambiguous", "none", "invalid"}:
        failures.append("CRISPR fixture envelope must preserve all four assignment outcomes")

    agreement_path = ROOT / "benchmarks" / "raw" / "count_agreement_summary.csv"
    with agreement_path.open(encoding="utf-8", newline="") as handle:
        agreement = {row["comparison"]: row for row in csv.DictReader(handle)}
    exact = agreement.get("dotmatch_exact_vs_mageck_exact", {})
    for key, value in {
        "status": "ok",
        "n_guides": "87437",
        "total_delta": "0",
        "differing_guides": "0",
        "max_abs_delta": "0",
        "pearson": "1.00000000",
        "spearman": "1.00000000",
    }.items():
        if exact.get(key) != value:
            failures.append(f"Yusa exact-count comparison {key} differs: {exact.get(key)!r}")

    perturb = json.loads((ROOT / "examples" / "perturb_seq_gse146194" / "expected-results.json").read_text(encoding="utf-8"))
    if perturb.get("guide_library", {}).get("target_count") != 32:
        failures.append("Perturb-seq reference must contain 32 guides")
    if perturb.get("input_slice", {}).get("evaluation_records") != 48000:
        failures.append("Perturb-seq reference must contain 48,000 held-out records")
    for run_id in ("k0", "k1"):
        run = perturb.get("runs", {}).get(run_id, {})
        if run.get("validation_mismatches") != 0 or run.get("agreement_rate") != 1.0:
            failures.append(f"Perturb-seq {run_id} must retain zero oracle differences")

    docs = "\n".join(
        (ROOT / name).read_text(encoding="utf-8")
        for name in (
            "docs/agent-crispr.md",
            "docs/agent-perturb-seq.md",
            "app/page.tsx",
            "agent-reference-crispr.json",
        )
    )
    for phrase in ["no downstream screen statistics", "guide-per-cell", "agent-reference-crispr.json", "82c61e9c2eb2d7701fd7505f4fa466e50fdbabbaebc16c266f304e06f4c424d9"]:
        if phrase.lower() not in docs.lower():
            failures.append(f"agent experience is missing reference boundary/evidence: {phrase}")

    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1
    print("CRISPR fixture hashes, Yusa exact-count agreement, and 32-guide Perturb-seq oracle records passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
