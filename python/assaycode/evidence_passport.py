"""Sequence Doctor research-evidence-passport pilot (Python 3.10+, stdlib only).

Run this file directly for an engine-independent source-checkout CLI. No shell,
network, raw-read export or workflow execution. Human identity is self-asserted,
not authenticated. Pin the passport digest OUTSIDE the passport to detect a
complete rewrite. See docs/grants/catalyst-2026/TECHNICAL.md for the trust boundary.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import sys
import zlib
from datetime import datetime, timezone
from typing import Any

VERSION = "0.1.0"
FORMAT = "dotmatch.research-evidence-passport/v1"
CANON = "dotmatch-json-v1"
MAX_JSON = 2 * 1024 * 1024
MAX_FILE = 1024 * 1024 * 1024
MAX_TOTAL = 2 * MAX_FILE
MAX_LINE = 1024 * 1024
MAX_RECORDS = 1000
MAX_PREFIX_BYTES = 8 * 1024 * 1024
HEX = re.compile(r"^[0-9a-f]{64}$")
IDENT = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
DECISIONS = ("biological_question", "sample_identity", "reference_choice",
             "study_design", "analysis_method", "reporting_limits")
LIMITS = [
    "Pre-execution evidence only; no biological result or claim is validated.",
    "FASTQ syntax inspection is bounded to the first 1000 records; hashes cover complete files.",
    "Human and client identities are self-asserted, not cryptographically authenticated.",
    "Events cover this server session, not all historical actions or external tools.",
    "An independently retained digest is required to detect wholesale passport replacement.",
]


class GateError(ValueError):
    """A refused operation, never an instruction to silently continue."""


def require(ok: bool, message: str) -> None:
    if not ok:
        raise GateError(message)


def text(value: Any, label: str, limit: int = 4096) -> str:
    require(isinstance(value, str) and 0 < len(value.strip()) <= limit,
            f"Invalid {label}")
    require(not any(ord(c) < 32 and c not in "\n\t" for c in value), f"Invalid {label}")
    return value


def keys(value: Any, required: set[str], optional: set[str] | None = None) -> None:
    require(isinstance(value, dict), "Expected an object")
    require(required <= value.keys() <= required | (optional or set()),
            "Missing or unsupported fields")


def canonical(value: Any) -> bytes:
    """Specified restricted JSON: no floats, safe integers, UTF-8, sorted keys.

    NOT RFC 8785. Unicode is preserved (not normalized). Empty/duplicate keys
    are handled by the surrounding schema/strict parser, not silently merged.
    """
    def check(v: Any, depth: int = 0) -> None:
        require(depth <= 64, "JSON nesting exceeds limit")
        if v is None or isinstance(v, bool):
            return
        if type(v) is int:
            require(abs(v) <= 9007199254740991, "Integer outside portable range")
        elif isinstance(v, str):
            v.encode("utf-8", errors="strict")
        elif isinstance(v, list):
            for item in v:
                check(item, depth + 1)
        elif isinstance(v, dict):
            require(all(isinstance(k, str) for k in v), "JSON keys must be strings")
            for k, item in v.items():
                check(k, depth + 1)
                check(item, depth + 1)
        else:
            raise GateError("Only restricted JSON values are supported; no floating-point numbers")
    check(value)
    return json.dumps(value, ensure_ascii=False, sort_keys=True,
                      separators=(",", ":"), allow_nan=False).encode("utf-8")


def digest(value: Any) -> str:
    return hashlib.sha256(canonical(value)).hexdigest()


def loads(raw: str | bytes) -> Any:
    require(len(raw) <= MAX_JSON, "JSON exceeds size limit")
    def pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in items:
            require(key not in result, "Duplicate JSON key")
            result[key] = value
        return result
    try:
        result = json.loads(raw, object_pairs_hook=pairs)
        canonical(result)
        return result
    except (UnicodeError, json.JSONDecodeError, RecursionError) as exc:
        raise GateError("Invalid JSON") from exc


def relative_path(value: Any) -> str:
    text(value, "relative path", 512)
    path = PurePosixPath(value)
    require(not path.is_absolute() and "\\" not in value and ":" not in value
            and all(p not in ("", ".", "..") for p in value.split("/")),
            "Paths must be project-relative without traversal")
    return value


def open_local(root: Path, relative: str) -> io.BufferedReader:
    """Open regular files below root using no-follow directory descriptors.

    POSIX pilot: no path-resolution/open race through intermediate symlinks.
    This is NOT a sandbox against another process with the same OS privileges.
    """
    relative_path(relative)
    require(os.name == "posix" and hasattr(os, "O_NOFOLLOW"),
            "This pilot requires POSIX no-follow file access")
    descriptor = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        parts = relative.split("/")
        for part in parts[:-1]:
            child = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                            dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
        fd = os.open(parts[-1], os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK,
                     dir_fd=descriptor)
        try:
            require(stat.S_ISREG(os.fstat(fd).st_mode), "Input is not a regular file")
            return os.fdopen(fd, "rb")
        except BaseException:
            os.close(fd)
            raise
    finally:
        os.close(descriptor)


def read_json(root: Path, name: str) -> Any:
    with open_local(root, name) as handle:
        return loads(handle.read(MAX_JSON + 1))


def file_evidence(root: Path, item: dict[str, Any]) -> dict[str, Any]:
    with open_local(root, item["path"]) as handle:
        before = os.fstat(handle.fileno())
        require(before.st_size <= MAX_FILE, "File exceeds pilot size limit")
        sha, size = hashlib.sha256(), 0
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            size += len(block)
            require(size <= MAX_FILE, "File exceeds pilot size limit")
            sha.update(block)
        qc: dict[str, Any] = {"kind": "not_applicable"}
        if item["role"] == "reads":
            handle.seek(0)
            source = gzip.GzipFile(fileobj=handle) if item["path"].endswith(".gz") else handle
            qc = {"kind": "fastq_prefix", "records_checked": 0,
                  "record_limit": MAX_RECORDS, "whole_fastq_validated": False}
            try:
                prefix_bytes = 0
                for _ in range(MAX_RECORDS):
                    header = source.readline(MAX_LINE + 1)
                    if not header:
                        break
                    seq, plus, quality = [source.readline(MAX_LINE + 1) for _ in range(3)]
                    prefix_bytes += sum(len(line) for line in (header, seq, plus, quality))
                    require(prefix_bytes <= MAX_PREFIX_BYTES, "FASTQ inspection exceeds prefix budget")
                    require(all(len(line) <= MAX_LINE for line in (header, seq, plus, quality)),
                            "FASTQ line exceeds pilot limit")
                    seq, quality = seq.rstrip(b"\r\n"), quality.rstrip(b"\r\n")
                    require(header.startswith(b"@") and plus.startswith(b"+")
                            and bool(seq) and len(seq) == len(quality)
                            and all(base in b"ACGTNacgtn" for base in seq)
                            and all(33 <= q <= 126 for q in quality), "Malformed FASTQ prefix")
                    qc["records_checked"] += 1
                require(qc["records_checked"] > 0, "Empty FASTQ")
            finally:
                if source is not handle:
                    source.close()
        after = os.fstat(handle.fileno())
        require((before.st_size, before.st_mtime_ns, before.st_ctime_ns) ==
                (after.st_size, after.st_mtime_ns, after.st_ctime_ns)
                and size == before.st_size, "File changed during inspection")
    return {**item, "observed_sha256": sha.hexdigest(), "size_bytes": size, "qc": qc}


def validate_project(project: Any) -> None:
    keys(project, {"schema", "project_id", "question", "study_design", "analysis_method",
                   "reporting_limits", "samples", "files"})
    require(project["schema"] == "dotmatch.project/v1", "Unsupported project schema")
    require(isinstance(project["project_id"], str) and bool(IDENT.fullmatch(project["project_id"])),
            "Invalid project identity")
    for field in ("question", "study_design", "analysis_method", "reporting_limits"):
        text(project[field], field)
    require(isinstance(project["samples"], list) and 1 <= len(project["samples"]) <= 128,
            "Declare 1-128 samples")
    sample_ids = set()
    for sample in project["samples"]:
        keys(sample, {"id", "biological_unit", "condition"})
        for field in sample:
            text(sample[field], field, 128)
        require(sample["id"] not in sample_ids, "Duplicate sample identifier")
        sample_ids.add(sample["id"])
    require(isinstance(project["files"], list) and 2 <= len(project["files"]) <= 128,
            "Declare 2-128 files including reads and reference")
    ids, paths, observed_samples, references = set(), set(), set(), 0
    for item in project["files"]:
        keys(item, {"id", "path", "role", "expected_sha256"}, {"sample_id", "version"})
        text(item["id"], "file identity", 128)
        relative_path(item["path"])
        require(item["id"] not in ids and item["path"] not in paths, "Duplicate file identity or path")
        ids.add(item["id"])
        paths.add(item["path"])
        require(isinstance(item["expected_sha256"], str)
                and bool(HEX.fullmatch(item["expected_sha256"])), "Expected SHA-256 required")
        require(item["role"] in ("reads", "reference", "sample_metadata"), "Unsupported file role")
        if item["role"] == "reads":
            require(item.get("sample_id") in sample_ids, "Reads require a declared sample")
            require(item["path"].endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")),
                    "Only FASTQ reads are supported by this pilot")
            observed_samples.add(item["sample_id"])
        if item["role"] == "reference":
            text(item.get("version"), "reference version", 128)
            references += 1
    require(references > 0 and observed_samples == sample_ids,
            "Every declared sample needs reads and the project needs a versioned reference")


def inspect_project(root: Path, project_name: str = "project.json") -> dict[str, Any]:
    project = read_json(root, project_name)
    validate_project(project)
    observed, findings, total = [], [], 0
    for item in project["files"]:
        try:
            evidence = file_evidence(root, item)
            total += evidence["size_bytes"]
            require(total <= MAX_TOTAL, "Project exceeds pilot size limit")
            observed.append(evidence)
            if evidence["observed_sha256"] != item["expected_sha256"]:
                findings.append({"code": "FILE_IDENTITY_MISMATCH", "file_id": item["id"],
                                 "severity": "block", "message": "File differs from the declared expected SHA-256"})
        except (OSError, ValueError, EOFError, zlib.error):
            # Do not copy raw reads, OS paths, or uncontrolled exception strings into evidence.
            findings.append({"code": "INPUT_UNVERIFIABLE", "file_id": item["id"],
                             "severity": "block", "message": "Input unavailable, unsafe, malformed or over budget"})
    body = {"schema": "dotmatch.sequence-doctor/v1", "project": project,
            "files": observed, "findings": findings,
            "status": "refused" if findings else "needs_human_review"}
    return {**body, "report_sha256": digest(body)}


def validate_report(report: Any) -> None:
    keys(report, {"schema", "project", "files", "findings", "status", "report_sha256"})
    require(report["schema"] == "dotmatch.sequence-doctor/v1", "Unsupported report schema")
    validate_project(report["project"])
    require(report["report_sha256"] == digest({k: v for k, v in report.items() if k != "report_sha256"}),
            "Report digest mismatch")
    require(isinstance(report["files"], list) and isinstance(report["findings"], list), "Invalid report lists")
    declarations = {f["id"]: f for f in report["project"]["files"]}
    seen = set()
    for item in report["files"]:
        require(isinstance(item, dict) and item.get("id") in declarations and item["id"] not in seen,
                "Invalid file evidence link")
        seen.add(item["id"])
        expected = declarations[item["id"]]
        keys(item, set(expected) | {"observed_sha256", "size_bytes", "qc"})
        require(all(item[k] == v for k, v in expected.items()), "File declaration mismatch")
        require(isinstance(item["observed_sha256"], str) and bool(HEX.fullmatch(item["observed_sha256"]))
                and type(item["size_bytes"]) is int and 0 <= item["size_bytes"] <= MAX_FILE,
                "Invalid file evidence")
        require(isinstance(item["qc"], dict), "Invalid QC record")
    for finding in report["findings"]:
        keys(finding, {"code", "file_id", "severity", "message"})
        require(finding["file_id"] in declarations and finding["severity"] == "block", "Invalid finding link")
        require(finding["code"] in ("FILE_IDENTITY_MISMATCH", "INPUT_UNVERIFIABLE"), "Unknown finding")
        text(finding["message"], "finding message")
    require(report["status"] == ("refused" if report["findings"] else "needs_human_review"), "Invalid report state")
    if not report["findings"]:
        require(seen == declarations.keys(), "Incomplete file evidence")
        require(all(f["observed_sha256"] == f["expected_sha256"] for f in report["files"]),
                "Unreported identity mismatch")
        require(sum(f["size_bytes"] for f in report["files"]) <= MAX_TOTAL, "Project exceeds pilot limit")
        for item in report["files"]:
            if item["role"] == "reads":
                qc = item["qc"]
                keys(qc, {"kind", "records_checked", "record_limit", "whole_fastq_validated"})
                require(qc["kind"] == "fastq_prefix" and type(qc["records_checked"]) is int
                        and 0 < qc["records_checked"] <= MAX_RECORDS
                        and qc["record_limit"] == MAX_RECORDS and qc["whole_fastq_validated"] is False,
                        "Invalid bounded FASTQ evidence")


def draft_contract(report: dict[str, Any]) -> dict[str, Any]:
    validate_report(report)
    project = report["project"]
    values = [project["question"], project["samples"],
              [f for f in project["files"] if f["role"] == "reference"],
              project["study_design"], project["analysis_method"], project["reporting_limits"]]
    return {"schema": "dotmatch.analysis-contract/v1", "project_id": project["project_id"],
            "report_sha256": report["report_sha256"],
            "decisions": {name: {"value": value, "confirmed": False, "rationale": ""}
                          for name, value in zip(DECISIONS, values)}, "approval": None}


def contract_body(contract: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in contract.items() if k != "approval"}


def validate_contract(contract: Any, report: dict[str, Any], approved: bool = True) -> None:
    validate_report(report)
    keys(contract, {"schema", "project_id", "report_sha256", "decisions", "approval"})
    expected = draft_contract(report)
    require(all(contract[k] == expected[k] for k in ("schema", "project_id", "report_sha256")),
            "Contract is stale or belongs to another inspection")
    keys(contract["decisions"], set(DECISIONS))
    for name in DECISIONS:
        decision = contract["decisions"][name]
        keys(decision, {"value", "confirmed", "rationale"})
        require(decision["value"] == expected["decisions"][name]["value"],
                "Decision differs from inspected project; edit project and inspect again")
        require(type(decision["confirmed"]) is bool and isinstance(decision["rationale"], str),
                "Invalid decision confirmation")
        if approved:
            require(decision["confirmed"], f"Human review required: {name}")
            text(decision["rationale"], "decision rationale")
    if approved:
        require(not report["findings"], "Blocking findings cannot be overridden by approval")
        approval = contract["approval"]
        keys(approval, {"reviewer", "role", "identity_assurance", "approved_at", "contract_sha256"})
        text(approval["reviewer"], "named reviewer", 256)
        require(approval["role"] == "human_reviewer"
                and approval["identity_assurance"] == "self_asserted_local_operator",
                "Unsupported approval identity assurance")
        timestamp(approval["approved_at"])
        require(approval["contract_sha256"] == digest(contract_body(contract)), "Approval no longer matches contract")


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def timestamp(value: Any) -> None:
    require(isinstance(value, str), "Invalid timestamp")
    try:
        datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ")
    except ValueError as exc:
        raise GateError("Invalid UTC timestamp") from exc


def approve_contract(contract: dict[str, Any], report: dict[str, Any], reviewer: str,
                     rationale: str) -> dict[str, Any]:
    """Human CLI surface, intentionally NOT exposed as an MCP tool."""
    validate_contract(contract, report, approved=False)
    require(not report["findings"], "Resolve blocking findings before approval")
    text(reviewer, "named reviewer", 256)
    text(rationale, "review rationale")
    result = copy.deepcopy(contract)
    for decision in result["decisions"].values():
        decision.update(confirmed=True, rationale=rationale)
    result["approval"] = {"reviewer": reviewer, "role": "human_reviewer",
                          "identity_assurance": "self_asserted_local_operator",
                          "approved_at": now(), "contract_sha256": digest(contract_body(result))}
    validate_contract(result, report)
    return result


def compilation_plan(contract: dict[str, Any], report: dict[str, Any]) -> dict[str, Any]:
    validate_contract(contract, report)
    return {"schema": "dotmatch.workflow-compilation-input/v1", "execution_status": "not_executed",
            "report_sha256": report["report_sha256"], "contract_sha256": digest(contract),
            "project_id": contract["project_id"], "engine": "external_adapter_required",
            "inputs": report["files"], "decisions": contract["decisions"]}


class Session:
    """Read-only project access and in-memory, hash-linked session evidence."""
    def __init__(self, root: Path, project_name: str = "project.json", client: dict[str, str] | None = None):
        self.root = root.resolve(strict=True)
        self.project_name = relative_path(project_name)
        self.client = client or {"name": "local-cli", "version": VERSION}
        self.reports: list[dict[str, Any]] = []
        self.events: list[dict[str, Any]] = []
        self.contract: dict[str, Any] | None = None
        self.plan: dict[str, Any] | None = None

    def event(self, kind: str, subject: str, details: dict[str, Any]) -> None:
        require(len(self.events) < 4096, "Session event limit reached; export and start a new session")
        body = {"index": len(self.events), "at": now(), "kind": kind,
                "subject_sha256": subject, "details": details,
                "previous_sha256": self.events[-1]["event_sha256"] if self.events else None}
        self.events.append({**body, "event_sha256": digest(body)})

    def inspect(self) -> dict[str, Any]:
        report = inspect_project(self.root, self.project_name)
        if not self.reports or report != self.reports[-1]:
            require(len(self.reports) < 256, "Session inspection limit reached")
            self.reports.append(report)
            self.contract = self.plan = None
        self.event("inspection", report["report_sha256"], {"status": report["status"]})
        if report["findings"]:
            self.event("refusal", report["report_sha256"], {"reason": "blocking_input_findings"})
        return report

    def check(self) -> dict[str, Any]:
        report = self.inspect()
        try:
            contract = read_json(self.root, "analysis-contract.json")
            validate_contract(contract, report)
        except (OSError, ValueError) as exc:
            self.contract = self.plan = None
            self.event("escalation", report["report_sha256"], {"reason": "human_review_required_or_stale"})
            raise GateError("Valid named-human approval of this exact inspection is required") from exc
        if self.contract is not None and self.contract != contract:
            self.plan = None
        self.contract = contract
        self.event("approval_observed", digest(contract), {"contract": copy.deepcopy(contract)})
        return {"status": "ready_for_workflow_compilation", "report_sha256": report["report_sha256"],
                "contract_sha256": digest(contract)}

    def compile(self) -> dict[str, Any]:
        self.check()  # Re-read bytes and approval immediately before authorizing compilation.
        assert self.contract is not None
        self.plan = compilation_plan(self.contract, self.reports[-1])
        self.event("compilation_input", digest(self.plan), {"plan": copy.deepcopy(self.plan)})
        return self.plan

    def export(self) -> dict[str, Any]:
        # Never label a previously successful session ready after source changes.
        report = self.inspect()
        if self.contract is not None:
            try:
                self.check()  # preserve the history even after approval is revoked
            except GateError:
                pass
        module_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        payload = {"schema": FORMAT, "scope": "pre_execution", "created_at": now(),
                   "producer": {"name": "Sequence Doctor", "version": VERSION, "source_sha256": module_hash},
                   "agent": {**self.client, "identity_assurance": "self_reported"},
                   "status": ("ready_for_workflow_compilation" if self.plan is not None
                              else report["status"]),
                   "reports": self.reports, "analysis_contract": self.contract,
                   "workflow_compilation_input": self.plan, "events": self.events,
                   "limitations": LIMITS}
        passport = seal(payload)
        verify_passport(passport)
        require(len(canonical(passport)) <= MAX_JSON, "Passport exceeds pilot export size limit")
        return copy.deepcopy(passport)


def seal(payload: dict[str, Any]) -> dict[str, Any]:
    return {"payload": copy.deepcopy(payload), "integrity": {"algorithm": "sha256",
            "canonicalization": CANON, "payload_sha256": digest(payload)}}


def verify_passport(passport: Any, expected_sha256: str | None = None) -> dict[str, Any]:
    keys(passport, {"payload", "integrity"})
    integrity, p = passport["integrity"], passport["payload"]
    keys(integrity, {"algorithm", "canonicalization", "payload_sha256"})
    require(integrity["algorithm"] == "sha256" and integrity["canonicalization"] == CANON,
            "Unsupported integrity format")
    actual = digest(p)
    require(integrity["payload_sha256"] == actual, "Passport payload digest mismatch")
    if expected_sha256 is not None:
        require(isinstance(expected_sha256, str) and bool(HEX.fullmatch(expected_sha256))
                and expected_sha256 == actual, "Independent expected digest mismatch")
    keys(p, {"schema", "scope", "created_at", "producer", "agent", "status", "reports",
             "analysis_contract", "workflow_compilation_input", "events", "limitations"})
    require(p["schema"] == FORMAT and p["scope"] == "pre_execution", "Unsupported passport schema or scope")
    timestamp(p["created_at"])
    keys(p["producer"], {"name", "version", "source_sha256"})
    require(isinstance(p["producer"]["source_sha256"], str)
            and bool(HEX.fullmatch(p["producer"]["source_sha256"])), "Invalid tool identity")
    for name in ("name", "version"):
        text(p["producer"][name], "producer identity", 128)
    keys(p["agent"], {"name", "version", "identity_assurance"})
    for name in ("name", "version"):
        text(p["agent"][name], "agent identity", 128)
    require(p["agent"]["identity_assurance"] == "self_reported" and p["limitations"] == LIMITS,
            "Unsupported assurance claim or omitted limitations")
    require(isinstance(p["reports"], list) and 1 <= len(p["reports"]) <= 256, "Invalid inspection history")
    for report in p["reports"]:
        validate_report(report)
    reports = {r["report_sha256"]: r for r in p["reports"]}
    project_ids = {r["project"]["project_id"] for r in p["reports"]}
    require(len(project_ids) == 1, "Mixed project history")
    current = p["reports"][-1]
    contract, plan = p["analysis_contract"], p["workflow_compilation_input"]
    if contract is not None:
        validate_contract(contract, current)
    if plan is not None:
        require(contract is not None and plan == compilation_plan(contract, current), "Compilation binding mismatch")
    require(p["status"] == ("ready_for_workflow_compilation" if plan is not None else current["status"]),
            "False passport readiness state")
    require(isinstance(p["events"], list) and 1 <= len(p["events"]) <= 4096, "Invalid event history")
    previous, seen_reports, compiled = None, set(), False
    approved_contracts: dict[str, Any] = {}
    inspection_order: list[str] = []
    for index, event in enumerate(p["events"]):
        keys(event, {"index", "at", "kind", "subject_sha256", "details", "previous_sha256", "event_sha256"})
        require(type(event["index"]) is int and event["index"] == index
                and event["previous_sha256"] == previous, "Broken event order or chain")
        timestamp(event["at"])
        require(event["event_sha256"] == digest({k: v for k, v in event.items() if k != "event_sha256"}),
                "Event digest mismatch")
        previous = event["event_sha256"]
        kind, subject, details = event["kind"], event["subject_sha256"], event["details"]
        require(isinstance(subject, str) and bool(HEX.fullmatch(subject)) and isinstance(details, dict),
                "Invalid event subject or details")
        require(kind in ("inspection", "refusal", "escalation", "approval_observed", "compilation_input"),
                "Unknown event kind")
        if kind == "inspection":
            require(subject in reports and details == {"status": reports[subject]["status"]}, "Unknown inspection")
            seen_reports.add(subject)
            if not inspection_order or inspection_order[-1] != subject:
                inspection_order.append(subject)
        elif kind in ("refusal", "escalation"):
            require(subject in seen_reports, "Decision precedes its inspection")
            if kind == "refusal":
                require(bool(reports[subject]["findings"]), "Refusal has no blocking evidence")
        elif kind == "approval_observed":
            keys(details, {"contract"})
            observed_contract = details["contract"]
            require(isinstance(observed_contract, dict), "Invalid approval snapshot")
            report_hash = observed_contract.get("report_sha256")
            require(report_hash in seen_reports and report_hash == inspection_order[-1],
                    "Approval is not for the latest inspection")
            validate_contract(observed_contract, reports[report_hash])
            require(subject == digest(observed_contract), "Approval event digest mismatch")
            approved_contracts[subject] = observed_contract
        elif kind == "compilation_input":
            keys(details, {"plan"})
            observed_plan = details["plan"]
            require(isinstance(observed_plan, dict), "Invalid compilation snapshot")
            contract_hash = observed_plan.get("contract_sha256")
            require(contract_hash in approved_contracts, "Compilation precedes approval")
            observed_contract = approved_contracts[contract_hash]
            report_hash = observed_contract["report_sha256"]
            require(report_hash == inspection_order[-1], "Compilation uses a superseded inspection")
            require(observed_plan == compilation_plan(observed_contract, reports[report_hash])
                    and subject == digest(observed_plan), "Compilation event binding mismatch")
            if plan is not None and subject == digest(plan):
                compiled = True
    require(inspection_order == [r["report_sha256"] for r in p["reports"]], "Inspection history order mismatch")
    require(seen_reports == reports.keys(), "Unlogged inspection")
    require(plan is None or compiled, "Missing approved compilation event")
    return {"integrity": "verified_against_independent_digest" if expected_sha256 else "internally_consistent_only",
            "payload_sha256": actual, "status": p["status"], "biological_validity": "not_established",
            "identity_authenticated": False, "external_file_bytes_rechecked": False}


TOOL_NAMES = ("inspect_project", "draft_analysis_contract", "check_analysis_contract",
              "compile_workflow", "export_evidence_passport")


def tools() -> list[dict[str, Any]]:
    descriptions = ["Inspect declared project files and refuse identity/syntax failures.",
                    "Draft six human-review decisions; this does NOT approve them.",
                    "Verify externally saved named-human approval against freshly inspected bytes.",
                    "Return approved engine-neutral compilation input; never execute a workflow.",
                    "Return a portable passport and session audit history; never write project files."]
    return [{"name": name, "description": description,
             "inputSchema": {"type": "object", "properties": {}, "additionalProperties": False},
             "annotations": {"readOnlyHint": True, "destructiveHint": False, "openWorldHint": False}}
            for name, description in zip(TOOL_NAMES, descriptions)]


def serve(root: Path, project_name: str = "project.json") -> None:
    session = Session(root, project_name)
    initialized, ready = False, False
    while True:
        raw = sys.stdin.buffer.readline(MAX_JSON + 1)
        if not raw:
            break
        if len(raw) > MAX_JSON:
            # Drain only the oversized frame; do not reinterpret its suffix as requests.
            while raw and not raw.endswith(b"\n"):
                raw = sys.stdin.buffer.readline(MAX_JSON + 1)
            sys.stdout.write('{"jsonrpc":"2.0","id":null,"error":{"code":-32700,"message":"Message exceeds size limit"}}\n')
            sys.stdout.flush()
            continue
        request_id: Any = None
        notification = False
        try:
            request = loads(raw)
            require(isinstance(request, dict), "JSON-RPC batches are not supported")
            request_id = request.get("id")
            notification = "id" not in request
            require(request_id is None or type(request_id) is int or isinstance(request_id, str), "Invalid request ID")
            require(request.get("jsonrpc") == "2.0" and isinstance(request.get("method"), str), "Invalid request")
            method, params = request["method"], request.get("params", {})
            require(isinstance(params, dict), "Invalid request parameters")
            if notification:
                if method == "notifications/initialized" and initialized:
                    ready = True
                continue
            if method == "initialize":
                require(not initialized, "Already initialized")
                client = params.get("clientInfo")
                require(isinstance(client, dict), "Client identity required")
                session.client = {k: text(client.get(k), "client identity", 128) for k in ("name", "version")}
                result: Any = {"protocolVersion": "2025-11-25", "capabilities": {"tools": {}},
                               "serverInfo": {"name": "dotmatch-sequence-doctor", "version": VERSION}}
                initialized = True
            elif method == "ping":
                result = {}
            elif method == "tools/list":
                require(ready, "Initialize the MCP session first")
                result = {"tools": tools()}
            elif method == "tools/call":
                require(ready, "Initialize the MCP session first")
                name = params.get("name")
                require(name in TOOL_NAMES and params.get("arguments", {}) == {}, "Unknown tool or arguments")
                try:
                    actions = {"inspect_project": session.inspect,
                               "draft_analysis_contract": lambda: draft_contract(session.inspect()),
                               "check_analysis_contract": session.check, "compile_workflow": session.compile,
                               "export_evidence_passport": session.export}
                    output = actions[name]()
                    result = {"content": [{"type": "text", "text": canonical(output).decode()}],
                              "structuredContent": output, "isError": False}
                except (OSError, ValueError) as exc:
                    result = {"content": [{"type": "text", "text": "Refused: input or approval validation failed"}],
                              "isError": True}
            else:
                raise GateError("Unsupported method")
            response = {"jsonrpc": "2.0", "id": request_id, "result": result}
        except (OSError, ValueError, TypeError, UnicodeError, RecursionError):
            if notification:
                continue
            response = {"jsonrpc": "2.0", "id": request_id, "error": {"code": -32600, "message": "Invalid or unsupported request"}}
        sys.stdout.write(canonical(response).decode() + "\n")
        sys.stdout.flush()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["inspect", "draft", "approve", "compile", "export", "verify", "serve"])
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--project", default="project.json")
    parser.add_argument("--reviewer")
    parser.add_argument("--rationale")
    parser.add_argument("--confirm-all-reviewed", action="store_true")
    parser.add_argument("--passport", type=Path)
    parser.add_argument("--expected-sha256")
    args = parser.parse_args(argv)
    try:
        if args.command == "serve":
            serve(args.root, args.project)
            return 0
        if args.command == "verify":
            require(args.passport is not None, "--passport required")
            with args.passport.open("rb") as handle:
                result = verify_passport(loads(handle.read(MAX_JSON + 1)), args.expected_sha256)
        else:
            session = Session(args.root, args.project)
            if args.command == "inspect":
                result = session.inspect()
            elif args.command == "draft":
                result = draft_contract(session.inspect())
            elif args.command == "approve":
                require(args.confirm_all_reviewed, "Explicit --confirm-all-reviewed required")
                result = approve_contract(read_json(session.root, "analysis-contract.json"), session.inspect(),
                                          args.reviewer, args.rationale)
            elif args.command == "compile":
                result = session.compile()
            else:
                try:
                    session.compile()
                except GateError:
                    pass
                result = session.export()
        print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
        return 2 if isinstance(result, dict) and result.get("status") == "refused" else 0
    except (OSError, ValueError, TypeError, UnicodeError, RecursionError) as exc:
        print(json.dumps({"status": "refused", "error": str(exc) if isinstance(exc, GateError)
                          else "Input could not be safely read or validated"}), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
