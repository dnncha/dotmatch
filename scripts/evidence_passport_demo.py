#!/usr/bin/env python3
"""Run the real stdio MCP gate with a scripted client and fictional reviewer.

No LLM, network, credentials, human data, paid compute or workflow execution.
An existing output directory is never overwritten.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "python/assaycode/evidence_passport.py"
spec = importlib.util.spec_from_file_location("passport_demo_core", MODULE)
assert spec and spec.loader
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)


def write(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def fixture(root: Path) -> bytes:
    root.mkdir(parents=True, exist_ok=False)
    reads = b"@synthetic-read-1\nACGT\n+\nIIII\n"
    reference = b"id\tsequence\ng1\tACGT\n"
    (root / "reads.fastq").write_bytes(reads)
    (root / "reference.tsv").write_bytes(b"id\tsequence\nwrong-library\tTTTT\n")
    project = {
        "schema": "dotmatch.project/v1", "project_id": "catalyst-synthetic-demo",
        "question": "Assign this synthetic library to its declared known targets.",
        "study_design": "One synthetic sample; no biological replication or inferential statistics.",
        "analysis_method": "Exact known-target matching; an external workflow adapter is required.",
        "reporting_limits": "Research-use demonstration only; no biological or clinical conclusion.",
        "samples": [{"id": "sample-1", "biological_unit": "synthetic-library-1", "condition": "demo"}],
        "files": [
            {"id": "reads-1", "path": "reads.fastq", "role": "reads", "sample_id": "sample-1",
             "expected_sha256": hashlib.sha256(reads).hexdigest()},
            {"id": "reference-1", "path": "reference.tsv", "role": "reference", "version": "synthetic-v1",
             "expected_sha256": hashlib.sha256(reference).hexdigest()},
        ],
    }
    write(root / "project.json", project)
    return reference


def run(out: Path) -> dict[str, object]:
    require_new = not out.exists()
    if not require_new:
        raise ValueError("Output directory already exists; choose a new path")
    reference = fixture(out / "project")
    process = subprocess.Popen([sys.executable, str(MODULE), "serve", "--root", str(out / "project")],
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, encoding="utf-8")
    transcript = []
    serial = 0
    assert process.stdin and process.stdout

    def request(method: str, params: dict | None = None) -> dict:
        nonlocal serial
        serial += 1
        envelope = {"jsonrpc": "2.0", "id": serial, "method": method, "params": params or {}}
        process.stdin.write(json.dumps(envelope) + "\n")
        process.stdin.flush()
        response = json.loads(process.stdout.readline())
        transcript.append({"request": envelope, "response": response})
        if "error" in response:
            raise AssertionError(response)
        return response["result"]

    def tool(name: str) -> dict:
        return request("tools/call", {"name": name, "arguments": {}})

    try:
        initialization = request("initialize", {"protocolVersion": "2025-11-25", "capabilities": {},
                                               "clientInfo": {"name": "scripted-grant-demo", "version": "1.0.0"}})
        assert initialization["protocolVersion"] == "2025-11-25"
        process.stdin.write('{"jsonrpc":"2.0","method":"notifications/initialized"}\n')
        process.stdin.flush()
        discovered = request("tools/list")["tools"]
        assert len(discovered) == 5 and all("approve" not in t["name"] for t in discovered)
        refused = tool("inspect_project")["structuredContent"]
        assert refused["status"] == "refused"
        assert any(f["code"] == "FILE_IDENTITY_MISMATCH" for f in refused["findings"])
        assert tool("compile_workflow")["isError"] is True
        refused_passport = tool("export_evidence_passport")["structuredContent"]
        write(out / "refused-evidence-passport.json", refused_passport)
        # The operator restores the independently expected synthetic reference bytes.
        (out / "project/reference.tsv").write_bytes(reference)
        corrected = tool("inspect_project")["structuredContent"]
        assert corrected["status"] == "needs_human_review"
        contract = tool("draft_analysis_contract")["structuredContent"]
        write(out / "project/analysis-contract.json", contract)
        assert tool("compile_workflow")["isError"] is True
        # Explicitly fictional human-approval ceremony: a separate CLI, NOT an MCP tool.
        approval = subprocess.run([
            sys.executable, str(MODULE), "approve", "--root", str(out / "project"),
            "--reviewer", "Dr Casey Example (fictional demo reviewer)",
            "--rationale", "Synthetic fixture only: reviewed all six decisions after restoring the expected reference.",
            "--confirm-all-reviewed"], check=True, capture_output=True, text=True)
        write(out / "project/analysis-contract.reviewed.json", json.loads(approval.stdout))
        os.replace(out / "project/analysis-contract.reviewed.json", out / "project/analysis-contract.json")
        compiled = tool("compile_workflow")["structuredContent"]
        assert compiled["execution_status"] == "not_executed"
        write(out / "workflow-compilation-input.json", compiled)
        passport = tool("export_evidence_passport")["structuredContent"]
        assert passport["payload"]["status"] == "ready_for_workflow_compilation"
        passport_path = out / "research-evidence-passport.json"
        write(passport_path, passport)
        anchor = passport["integrity"]["payload_sha256"]
        (out / "research-evidence-passport.payload-sha256.txt").write_text(anchor + "\n")
        verified = subprocess.run([sys.executable, str(MODULE), "verify", "--passport", str(passport_path),
                                   "--expected-sha256", anchor], check=True, capture_output=True, text=True)
        verification = json.loads(verified.stdout)
        # Recompute the embedded checksum after tampering: independent pin still rejects it.
        altered = copy.deepcopy(passport["payload"])
        altered["agent"]["name"] = "changed-after-export"
        write(out / "tampered-and-rehashed.json", core.seal(altered))
        tampered = subprocess.run([sys.executable, str(MODULE), "verify", "--passport", str(out / "tampered-and-rehashed.json"),
                                  "--expected-sha256", anchor], capture_output=True, text=True)
        assert tampered.returncode == 2
        summary = {
            "scenario": "deliberately_wrong_synthetic_reference", "client": "real_stdio_scripted_client_not_an_LLM",
            "reviewer": "fictional_demonstration_not_a_real_human_attestation",
            "wrong_reference_refused": True, "unapproved_compilation_refused": True,
            "corrected_and_approved_compilation_input_returned": True,
            "historical_refusal_retained": any(e["kind"] == "refusal" for e in passport["payload"]["events"]),
            "rehashed_tampering_rejected_against_pin": True,
            "biological_workflow_executed": False, "external_client_certification": False,
            "payload_sha256": anchor, "verification": verification,
            "note": "The adjacent digest file is convenient, not independent storage. Retain the digest separately.",
        }
        write(out / "demo-result.json", summary)
        write(out / "mcp-transcript.json", transcript)
        return summary
    finally:
        process.stdin.close()
        process.wait(timeout=10)
        if process.returncode:
            raise RuntimeError("MCP process failed: " + (process.stderr.read() if process.stderr else ""))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    print(json.dumps(run(parser.parse_args().out), indent=2))
