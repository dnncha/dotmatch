"""Versioned, host-neutral tools for local research agents.

The public functions in this module deliberately accept and return JSON-shaped
objects.  They never accept command strings and never make network requests.
Scientific execution is delegated to the existing AssaySpec engine so its
ambiguity, target-safety, reliability, and CPU-authority rules remain binding.
"""

from __future__ import annotations

import argparse
import contextlib
import copy
import hashlib
import io
import json
import os
import re
import resource
import shutil
import sys
import time
from importlib import resources
from pathlib import Path
from typing import Any, Mapping, Sequence

from . import __version__
from .assayspec import (
    AssaySpecError,
    command_assay,
    compile_assay_plan,
    format_plan,
    load_assay_spec,
    scaffold_assay_project,
)


TOOL_CONTRACT_VERSION = "1.0"
ENVELOPE_SCHEMA_VERSION = "1.0"
DEFAULT_MINIMUM_FREE_BYTES = 1024 * 1024 * 1024
SUPPORTED_INTENTS = {
    "crispr-guide-counting": "crispr",
    "perturb-seq-guide-capture": "feature-barcode",
}
STATUS_EXIT_CODES = {
    "passed": 0,
    "needs_review": 1,
    "failed": 2,
    "blocked": 2,
    "invalid_input": 2,
    "interrupted": 130,
}
ALLOWED_REMEDIATIONS = {
    ("extract", "start"),
    ("extract", "length"),
    ("extract", "orientation"),
    ("assignment", "metric"),
    ("assignment", "k"),
    ("backend", "mode"),
}
FORBIDDEN_REMEDIATION_SECTIONS = {"targets", "reliability", "outputs"}


class AgentToolError(ValueError):
    """A structured-input or local safety failure."""


class AgentResourceError(AgentToolError):
    """A local resource gate that must block execution."""


def _package_json(name: str) -> dict[str, Any]:
    path = resources.files("dotmatch").joinpath("data", name)
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise AgentToolError(f"installed {name} must contain a JSON object")
    return value


def list_tools() -> dict[str, Any]:
    """Return a detached copy of the installed agent-tool contract."""

    return copy.deepcopy(_package_json("agent-tools.json"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _artifact(path: Path, *, role: str, safe_to_share: bool = True) -> dict[str, Any]:
    resolved = path.resolve(strict=False)
    result: dict[str, Any] = {
        "role": role,
        "path": str(resolved),
        "exists": path.is_file(),
        "safe_to_share": bool(safe_to_share),
    }
    if path.is_file():
        result.update({"sha256": _sha256(path), "bytes": path.stat().st_size})
    return result


def _resource_snapshot() -> tuple[float, float, int]:
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return time.monotonic(), float(usage.ru_utime + usage.ru_stime), int(usage.ru_maxrss)


def _resource_use(start: tuple[float, float, int], anchor: Path | None = None) -> dict[str, Any]:
    now = _resource_snapshot()
    target = anchor if anchor and anchor.exists() else Path.cwd()
    disk = shutil.disk_usage(target)
    return {
        "wall_seconds": round(now[0] - start[0], 6),
        "cpu_seconds": round(now[1] - start[1], 6),
        "max_rss_platform_units": max(now[2], start[2]),
        "disk_free_bytes": int(disk.free),
    }


def _base_envelope(tool: str) -> dict[str, Any]:
    return {
        "schema_version": ENVELOPE_SCHEMA_VERSION,
        "tool_contract_version": TOOL_CONTRACT_VERSION,
        "dotmatch_version": __version__,
        "tool": tool,
        "status": "passed",
        "exit_code": 0,
        "spec": {"revision": 0, "path": "", "sha256": ""},
        "artifacts": [],
        "findings": [],
        "resource_use": {},
        "next_actions": [],
    }


def _set_status(envelope: dict[str, Any], status: str) -> None:
    envelope["status"] = status
    envelope["exit_code"] = STATUS_EXIT_CODES[status]


def _finding(
    finding_id: str,
    severity: str,
    message: str,
    *,
    source_artifact: str = "",
) -> dict[str, str]:
    return {
        "finding_id": finding_id,
        "severity": severity,
        "message": message,
        "source_artifact": source_artifact,
    }


def _as_object(value: Mapping[str, Any] | object) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise AgentToolError("tool input must be a JSON object")
    return dict(value)


def _validate_keys(data: Mapping[str, Any], *, required: set[str], allowed: set[str]) -> None:
    missing = sorted(required - set(data))
    unknown = sorted(set(data) - allowed)
    if missing:
        raise AgentToolError("missing required input field(s): " + ", ".join(missing))
    if unknown:
        raise AgentToolError("unknown input field(s): " + ", ".join(unknown))


def _input_path(value: object, label: str, *, directory: bool = False) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise AgentToolError(f"{label} must be a non-empty path string")
    path = Path(value).expanduser().resolve(strict=True)
    if directory and not path.is_dir():
        raise AgentToolError(f"{label} must be a directory: {path}")
    if not directory and not path.is_file():
        raise AgentToolError(f"{label} must be a file: {path}")
    return path


def _empty_output_path(value: object, label: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise AgentToolError(f"{label} must be a non-empty path string")
    requested = Path(value).expanduser()
    if requested.is_symlink():
        raise AgentToolError(f"{label} must not be a symlink: {requested}")
    path = requested.resolve(strict=False)
    if path.exists() and not path.is_dir():
        raise AgentToolError(f"{label} must be a directory: {path}")
    if path.exists() and any(path.iterdir()):
        raise AgentToolError(f"{label} must be empty: {path}")
    parent = path.parent.resolve(strict=True)
    if not parent.is_dir():
        raise AgentToolError(f"{label} parent must be a directory: {parent}")
    return path


def _integer_parameter(
    value: object,
    label: str,
    *,
    minimum: int,
    maximum: int | None = None,
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise AgentToolError(f"{label} must be an integer")
    if value < minimum or (maximum is not None and value > maximum):
        if maximum is None:
            raise AgentToolError(f"{label} must be at least {minimum}")
        raise AgentToolError(f"{label} must be between {minimum} and {maximum}")
    return value


def _boolean_parameter(value: object, label: str) -> bool:
    if not isinstance(value, bool):
        raise AgentToolError(f"{label} must be a boolean")
    return value


def _require_resources(anchor: Path, minimum_free_bytes: object) -> None:
    minimum = _integer_parameter(
        minimum_free_bytes,
        "minimum_free_bytes",
        minimum=0,
    )
    free = shutil.disk_usage(anchor).free
    if free < minimum:
        raise AgentResourceError(
            f"insufficient disk: {free} bytes free; {minimum} bytes required before local assay execution"
        )


def _spec_record(path: Path, revision: int) -> dict[str, Any]:
    return {"revision": revision, "path": str(path.resolve()), "sha256": _sha256(path)}


def _load_reliability(spec_path: Path) -> tuple[dict[str, Any], dict[str, Path]]:
    plan = compile_assay_plan(load_assay_spec(spec_path))
    summary_path = plan.artifacts["reliability_summary"]
    if not summary_path.is_file():
        return {}, plan.artifacts
    value = json.loads(summary_path.read_text(encoding="utf-8"))
    return (value if isinstance(value, dict) else {}), plan.artifacts


def _normalized_findings(summary: Mapping[str, Any]) -> list[dict[str, Any]]:
    fields = (
        "finding_id",
        "severity",
        "stage",
        "sample_id",
        "metric",
        "observed",
        "threshold",
        "message",
        "recommended_action",
        "source_artifact",
    )
    result: list[dict[str, Any]] = []
    for item in summary.get("findings", []) or []:
        if isinstance(item, Mapping):
            result.append({key: item.get(key, "") for key in fields})
    return sorted(
        result,
        key=lambda item: (str(item["severity"]), str(item["finding_id"]), str(item["sample_id"])),
    )


def _artifact_set(paths: Mapping[str, Path]) -> list[dict[str, Any]]:
    independently_shareable_roles = {"citation_bib", "software_versions"}
    return [
        _artifact(path, role=role, safe_to_share=role in independently_shareable_roles)
        for role, path in sorted(paths.items())
        if path.is_file()
    ]


def _captured_assay(command: Sequence[str]) -> tuple[int, str, str]:
    stdout = io.StringIO()
    stderr = io.StringIO()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
        rc = command_assay(command)
    diagnostics = stderr.getvalue()
    if diagnostics:
        print(diagnostics, end="", file=sys.stderr)
    return rc, stdout.getvalue(), diagnostics


def _tool_discover(data: dict[str, Any], envelope: dict[str, Any]) -> None:
    _validate_keys(data, required=set(), allowed={"intent"})
    intent = data.get("intent")
    if intent is not None and intent not in SUPPORTED_INTENTS:
        raise AgentToolError("intent must be one of: " + ", ".join(sorted(SUPPORTED_INTENTS)))
    capabilities = _package_json("agent-capabilities.json")
    contract = list_tools()
    envelope["result"] = {
        "local_only": True,
        "network_requests": "none",
        "supported_intents": sorted(SUPPORTED_INTENTS),
        "selected_intent": intent or "",
        "capabilities": capabilities,
        "tools": contract["tools"],
    }
    envelope["next_actions"] = [
        {"tool": "prepare_assay", "reason": "Scaffold a reviewed local assay from targets and FASTQ files."}
    ]


def _tool_prepare(data: dict[str, Any], envelope: dict[str, Any]) -> None:
    allowed = {"intent", "targets", "reads_dir", "output_dir", "link_reads", "threads", "max_reads", "max_start", "minimum_free_bytes"}
    _validate_keys(data, required={"intent", "targets", "reads_dir", "output_dir"}, allowed=allowed)
    intent = str(data["intent"])
    if intent not in SUPPORTED_INTENTS:
        raise AgentToolError("intent must be one of: " + ", ".join(sorted(SUPPORTED_INTENTS)))
    targets = _input_path(data["targets"], "targets")
    reads_dir = _input_path(data["reads_dir"], "reads_dir", directory=True)
    output_dir = _empty_output_path(data["output_dir"], "output_dir")
    _require_resources(output_dir.parent, data.get("minimum_free_bytes", DEFAULT_MINIMUM_FREE_BYTES))
    threads = _integer_parameter(data.get("threads", 1), "threads", minimum=1)
    max_reads = _integer_parameter(data.get("max_reads", 50000), "max_reads", minimum=1)
    max_start = _integer_parameter(data.get("max_start", 32), "max_start", minimum=0)
    link_reads = _boolean_parameter(data.get("link_reads", False), "link_reads")
    result = scaffold_assay_project(
        template=SUPPORTED_INTENTS[intent],
        project_dir=output_dir,
        reads_dir=reads_dir,
        targets=targets,
        link_reads=link_reads,
        threads=threads,
        max_reads=max_reads,
        max_start=max_start,
    )
    spec_path = result["spec"]
    assay = load_assay_spec(spec_path)
    plan = compile_assay_plan(assay)
    envelope["spec"] = _spec_record(spec_path, 0)
    envelope["artifacts"] = [
        _artifact(path, role=role, safe_to_share=False)
        for role, path in sorted(result.items())
        if path.is_file()
    ]
    envelope["result"] = {
        "intent": intent,
        "assay_status": assay.status,
        "plan": format_plan(plan).splitlines(),
        "output_directory": str(output_dir),
    }
    if threads > 1:
        envelope["findings"] = [
            _finding(
                "aggregate_diagnostics_only",
                "info",
                "Row-level diagnostic files are disabled for multi-threaded counting; aggregate ambiguity and assignment counts remain available.",
                source_artifact=str(result["report"]),
            )
        ]
    if assay.status == "draft":
        _set_status(envelope, "needs_review")
        envelope["findings"].append(
            _finding(
                "draft_assayspec",
                "warning",
                "Inference was not decisive; review the candidate window before execution.",
                source_artifact=str(result["report"]),
            )
        )
        envelope["next_actions"] = [
            {"tool": "inspect_assay", "input": {"spec": str(spec_path)}, "reason": "Review preflight evidence."}
        ]
    else:
        envelope["next_actions"] = [
            {"tool": "inspect_assay", "input": {"spec": str(spec_path)}, "reason": "Run target and resource preflight."}
        ]


def _tool_inspect(data: dict[str, Any], envelope: dict[str, Any]) -> None:
    _validate_keys(data, required={"spec"}, allowed={"spec", "minimum_free_bytes"})
    spec_path = _input_path(data["spec"], "spec")
    _require_resources(spec_path.parent, data.get("minimum_free_bytes", DEFAULT_MINIMUM_FREE_BYTES))
    assay = load_assay_spec(spec_path)
    plan = compile_assay_plan(assay)
    rc, _stdout, _stderr = _captured_assay(["check", str(spec_path)])
    summary, artifacts = _load_reliability(spec_path)
    status = str(summary.get("overall_status", "passed" if rc == 0 else "blocked"))
    _set_status(envelope, status if status in STATUS_EXIT_CODES else "blocked")
    envelope["spec"] = _spec_record(spec_path, _revision_from_path(spec_path))
    envelope["artifacts"] = _artifact_set(artifacts)
    envelope["findings"] = _normalized_findings(summary)
    envelope["result"] = {"plan": format_plan(plan).splitlines(), "reliability": summary}
    envelope["next_actions"] = (
        [{"tool": "run_assay", "input": {"spec": str(spec_path)}, "reason": "Preflight permits local execution."}]
        if envelope["status"] in {"passed", "needs_review"}
        else [{"tool": "inspect_assay", "reason": "Resolve blocking findings without changing targets or QC thresholds."}]
    )


def _revision_from_path(path: Path) -> int:
    match = re.search(r"\.agent-r([0-9]+)\.toml$", path.name)
    return int(match.group(1)) if match else 0


def _finding_state(summary: Mapping[str, Any]) -> str:
    rows = [
        (str(item.get("finding_id", "")), str(item.get("severity", "")), str(item.get("observed", "")))
        for item in summary.get("findings", []) or []
        if isinstance(item, Mapping)
    ]
    return hashlib.sha256(json.dumps(sorted(rows)).encode("utf-8")).hexdigest()


def _scientific_spec_state(path: Path) -> str:
    data = copy.deepcopy(load_assay_spec(path).data)
    run = data.get("run")
    if isinstance(run, dict):
        run.pop("out_dir", None)
    return hashlib.sha256(json.dumps(data, sort_keys=True, default=str).encode("utf-8")).hexdigest()


def _toml_scalar(raw: str) -> str:
    value = raw.strip()
    if value in {"true", "false"} or re.fullmatch(r"-?[0-9]+(?:\.[0-9]+)?", value):
        return value
    if value in {"hamming", "levenshtein", "cpu", "forward", "reverse_complement"}:
        return json.dumps(value)
    raise AgentToolError(f"unsupported remediation value: {raw}")


def _replace_toml_value(text: str, section: str, key: str, suggested: str) -> str:
    section_re = re.compile(rf"(?ms)(^\[{re.escape(section)}\]\s*$)(.*?)(?=^\[|\Z)")
    match = section_re.search(text)
    if not match:
        raise AgentToolError(f"cannot apply remediation; missing [{section}] section")
    body = match.group(2)
    key_re = re.compile(rf"(?m)^(\s*{re.escape(key)}\s*=\s*).*$")
    scalar = _toml_scalar(suggested)
    if key_re.search(body):
        body = key_re.sub(rf"\g<1>{scalar}", body, count=1)
    else:
        body = body.rstrip() + f"\n{key} = {scalar}\n"
    return text[: match.start(2)] + body + text[match.end(2) :]


def _replace_out_dir(text: str, out_dir: Path) -> str:
    section_re = re.compile(r"(?ms)(^\[run\]\s*$)(.*?)(?=^\[|\Z)")
    match = section_re.search(text)
    if not match:
        raise AgentToolError("cannot revise spec; missing [run] section")
    body = match.group(2)
    key_re = re.compile(r"(?m)^(\s*out_dir\s*=\s*).*$")
    quoted = json.dumps(str(out_dir))
    if key_re.search(body):
        body = key_re.sub(rf"\g<1>{quoted}", body, count=1)
    else:
        body = body.rstrip() + f"\nout_dir = {quoted}\n"
    return text[: match.start(2)] + body + text[match.end(2) :]


def _allowed_fixes(summary: Mapping[str, Any]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    allowed: list[dict[str, str]] = []
    forbidden: list[dict[str, str]] = []
    for raw in summary.get("assay_fixes", []) or []:
        if not isinstance(raw, Mapping):
            continue
        fix = {str(key): str(value) for key, value in raw.items()}
        section_key = (fix.get("section", ""), fix.get("key", ""))
        if section_key not in ALLOWED_REMEDIATIONS or fix.get("section") in FORBIDDEN_REMEDIATION_SECTIONS:
            forbidden.append(fix)
            continue
        if section_key == ("backend", "mode") and fix.get("suggested_value") != "cpu":
            forbidden.append(fix)
            continue
        allowed.append(fix)
    return allowed, forbidden


def _write_candidate_spec(source: Path, fixes: Sequence[Mapping[str, str]], revision: int) -> Path:
    text = source.read_text(encoding="utf-8")
    for fix in fixes:
        text = _replace_toml_value(text, str(fix["section"]), str(fix["key"]), str(fix["suggested_value"]))
    candidate = source.with_name(re.sub(r"(?:\.agent-r[0-9]+)?\.toml$", "", source.name) + f".agent-r{revision}.toml")
    candidate_out = candidate.parent / f"assay_out.agent-r{revision}"
    text = _replace_out_dir(text, candidate_out)
    try:
        with candidate.open("x", encoding="utf-8") as handle:
            handle.write(text)
    except FileExistsError as exc:
        raise AgentToolError(f"refusing to overwrite existing candidate spec: {candidate}") from exc
    load_assay_spec(candidate)
    return candidate


def _tool_run(data: dict[str, Any], envelope: dict[str, Any]) -> None:
    _validate_keys(data, required={"spec"}, allowed={"spec", "max_revisions", "minimum_free_bytes"})
    current = _input_path(data["spec"], "spec")
    _require_resources(current.parent, data.get("minimum_free_bytes", DEFAULT_MINIMUM_FREE_BYTES))
    max_revisions = _integer_parameter(
        data.get("max_revisions", 3),
        "max_revisions",
        minimum=0,
        maximum=3,
    )
    revision = _revision_from_path(current)
    seen: set[tuple[str, str]] = set()
    history: list[dict[str, Any]] = []

    initial_plan = compile_assay_plan(load_assay_spec(current))
    completed_markers = [
        initial_plan.artifacts.get(name)
        for name in ("manifest", "counts", "summary", "sample_qc", "assignments")
    ]
    existing_completed = [path for path in completed_markers if path is not None and path.exists()]
    if existing_completed:
        raise AgentToolError(
            "refusing to overwrite completed assay outputs: "
            + ", ".join(str(path) for path in existing_completed)
        )

    while True:
        rc, _stdout, _stderr = _captured_assay(["start", str(current)])
        summary, artifacts = _load_reliability(current)
        status = str(summary.get("overall_status", "passed" if rc == 0 else "blocked"))
        file_hash = _sha256(current)
        state = (_scientific_spec_state(current), _finding_state(summary))
        history.append(
            {
                "revision": revision,
                "spec_sha256": file_hash,
                "scientific_state_sha256": state[0],
                "finding_state_sha256": state[1],
                "status": status,
                "exit_code": rc,
            }
        )
        if status == "passed" and rc == 0:
            break
        if state in seen:
            status = "blocked"
            summary.setdefault("findings", []).append(
                _finding("repeated_spec_finding_state", "blocked", "Revision stopped because the same spec and finding state repeated.")
            )
            break
        seen.add(state)
        allowed, forbidden = _allowed_fixes(summary)
        if forbidden:
            summary.setdefault("findings", []).append(
                _finding(
                    "forbidden_remediation",
                    "blocked",
                    "DotMatch refused a suggested change to targets, output policy, scientific scope, or reliability thresholds.",
                )
            )
        if not allowed or revision >= max_revisions:
            if revision >= max_revisions:
                summary.setdefault("findings", []).append(
                    _finding("revision_limit_exhausted", "blocked", "The assay reached the configured revision limit.")
                )
                status = "blocked"
            break
        try:
            _require_resources(current.parent, data.get("minimum_free_bytes", DEFAULT_MINIMUM_FREE_BYTES))
        except AgentResourceError as exc:
            summary.setdefault("findings", []).append(
                _finding("insufficient_resources", "blocked", str(exc))
            )
            status = "blocked"
            break
        revision += 1
        current = _write_candidate_spec(current, allowed, revision)

    if status not in STATUS_EXIT_CODES:
        status = "blocked"
    _set_status(envelope, status)
    envelope["spec"] = _spec_record(current, revision)
    envelope["artifacts"] = _artifact_set(artifacts)
    envelope["findings"] = _normalized_findings(summary)
    envelope["result"] = {"revision_history": history, "reliability": summary}
    envelope["next_actions"] = (
        [{"tool": "review_assay", "input": {"spec": str(current)}, "reason": "Normalize the completed evidence record."}]
        if status == "passed"
        else [{"tool": "inspect_assay", "input": {"spec": str(current)}, "reason": "Resolve the remaining reliability block explicitly."}]
    )


def _tool_review(data: dict[str, Any], envelope: dict[str, Any]) -> None:
    _validate_keys(data, required={"spec"}, allowed={"spec", "minimum_free_bytes"})
    spec_path = _input_path(data["spec"], "spec")
    summary, artifacts = _load_reliability(spec_path)
    if not summary:
        raise AgentToolError("review requires reliability_summary.json; run or inspect the assay first")
    status = str(summary.get("overall_status", "blocked"))
    _set_status(envelope, status if status in STATUS_EXIT_CODES else "blocked")
    envelope["spec"] = _spec_record(spec_path, _revision_from_path(spec_path))
    envelope["artifacts"] = _artifact_set(artifacts)
    envelope["findings"] = _normalized_findings(summary)
    envelope["result"] = {
        "reliability": summary,
        "evidence_boundary": summary.get("evidence_boundary", {}),
        "ambiguity_policy": load_assay_spec(spec_path).data.get("assignment", {}).get("ambiguous", "discard"),
    }
    envelope["next_actions"] = (
        [{"tool": "handoff_assay", "reason": "Create a hashed, raw-data-free review bundle."}]
        if status == "passed"
        else [{"tool": "inspect_assay", "reason": "Do not hand off failed results as production-ready evidence."}]
    )


def _tool_handoff(data: dict[str, Any], envelope: dict[str, Any]) -> None:
    _validate_keys(
        data,
        required={"spec", "output_dir"},
        allowed={"spec", "output_dir", "minimum_free_bytes"},
    )
    spec_path = _input_path(data["spec"], "spec")
    out_dir = _empty_output_path(data["output_dir"], "output_dir")
    _require_resources(out_dir.parent, data.get("minimum_free_bytes", DEFAULT_MINIMUM_FREE_BYTES))
    rc, _stdout, _stderr = _captured_assay(["handoff", str(spec_path), "--out-dir", str(out_dir)])
    if rc:
        raise AgentToolError("handoff failed; inspect diagnostics on stderr")
    manifest = out_dir / "handoff_manifest.json"
    sums = out_dir / "SHA256SUMS"
    reliability, _artifacts = _load_reliability(spec_path)
    reliability_status = str(reliability.get("overall_status", "blocked"))
    _set_status(envelope, reliability_status if reliability_status in STATUS_EXIT_CODES else "blocked")
    envelope["spec"] = _spec_record(spec_path, _revision_from_path(spec_path))
    envelope["artifacts"] = [
        _artifact(manifest, role="handoff_manifest", safe_to_share=False),
        _artifact(sums, role="handoff_checksums", safe_to_share=False),
    ]
    envelope["result"] = {
        "output_directory": str(out_dir),
        "raw_data_included": False,
        "manifest": json.loads(manifest.read_text(encoding="utf-8")),
        "reliability_status": reliability_status,
    }
    envelope["next_actions"] = [
        {
            "action": "review_handoff",
            "reason": (
                "Verify hashes and reliability boundaries before downstream interpretation."
                if reliability_status == "passed"
                else "The bundle was created for review but is not production-ready; resolve its reliability findings."
            ),
        }
    ]


_TOOL_HANDLERS = {
    "discover": _tool_discover,
    "prepare_assay": _tool_prepare,
    "inspect_assay": _tool_inspect,
    "run_assay": _tool_run,
    "review_assay": _tool_review,
    "handoff_assay": _tool_handoff,
}


def invoke_tool(tool: str, input_data: Mapping[str, Any]) -> dict[str, Any]:
    """Invoke one local tool and always return the stable JSON envelope."""

    start = _resource_snapshot()
    envelope = _base_envelope(tool if isinstance(tool, str) else "")
    try:
        if not isinstance(tool, str):
            raise AgentToolError("tool must be a string")
        if tool not in _TOOL_HANDLERS:
            raise AgentToolError("unknown tool; choose one of: " + ", ".join(sorted(_TOOL_HANDLERS)))
        data = _as_object(input_data)
        _TOOL_HANDLERS[tool](data, envelope)
    except KeyboardInterrupt:
        _set_status(envelope, "interrupted")
        envelope["findings"] = [_finding("interrupted_execution", "blocked", "Execution was interrupted; partial outputs are not trusted.")]
        envelope["next_actions"] = [{"tool": "inspect_assay", "reason": "Inspect partial artifacts before a clean rerun."}]
    except (AgentToolError, AssaySpecError, FileNotFoundError, json.JSONDecodeError, OSError, ValueError) as exc:
        if isinstance(exc, AgentResourceError):
            error_status = "blocked"
        elif isinstance(exc, (AgentToolError, json.JSONDecodeError, ValueError)):
            error_status = "invalid_input"
        else:
            error_status = "blocked"
        _set_status(envelope, error_status)
        envelope["findings"] = [_finding("tool_error", "blocked", str(exc))]
        envelope["next_actions"] = [{"tool": "discover", "reason": "Inspect the installed structured contract and required inputs."}]
    envelope["resource_use"] = _resource_use(start)
    return envelope


def _read_input(source: str) -> dict[str, Any]:
    if source == "-":
        raw = sys.stdin.read()
    else:
        raw = _input_path(source, "input").read_text(encoding="utf-8")
    value = json.loads(raw)
    if not isinstance(value, dict):
        raise AgentToolError("input JSON must contain an object")
    return value


def _skill_source(host: str) -> Any:
    directory = {
        "codex": "codex-skill",
        "claude-code": "claude-code-skill",
    }.get(host)
    if directory is None:
        raise AgentToolError(f"unsupported skill host: {host}")
    return resources.files("dotmatch").joinpath("data", directory)


def _export_skill(target: str, host: str = "codex") -> dict[str, Any]:
    out_dir = _empty_output_path(target, "target")
    source = _skill_source(host)
    out_dir.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        destination = out_dir / item.name
        if item.is_dir():
            shutil.copytree(str(item), destination)
        else:
            destination.write_bytes(item.read_bytes())
    files = [_artifact(path, role="skill_file") for path in sorted(out_dir.rglob("*")) if path.is_file()]
    envelope = _base_envelope("export-skill")
    envelope["artifacts"] = files
    envelope["result"] = {"target": str(out_dir), "host": host, "files": len(files)}
    restart_action = "restart_claude_code" if host == "claude-code" else "restart_codex"
    envelope["next_actions"] = [{"action": restart_action, "reason": "Reload local skill discovery after installation."}]
    return envelope


class _AgentArgumentParser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        raise AgentToolError(message)


def command_agent(argv: Sequence[str]) -> int:
    """CLI namespace for ``dotmatch agent``."""

    parser = _AgentArgumentParser(prog="dotmatch agent", description="Invoke versioned local research-agent tools.")
    sub = parser.add_subparsers(dest="command", required=True, parser_class=_AgentArgumentParser)
    tools_parser = sub.add_parser("tools", help="print the installed tool contract")
    tools_parser.add_argument("--json", action="store_true", required=True)
    invoke = sub.add_parser("invoke", help="invoke one structured local tool")
    invoke.add_argument("tool")
    invoke.add_argument("--input", required=True)
    export = sub.add_parser("export-skill", help="copy a bundled agent skill into an empty directory")
    export.add_argument("--target", required=True)
    export.add_argument("--host", choices=("codex", "claude-code"), default="codex")
    command = "agent"
    try:
        args = parser.parse_args(list(argv))
        command = str(args.command)
        if args.command == "tools":
            result = list_tools()
            print(json.dumps(result, indent=2, sort_keys=True))
            return 0
        if args.command == "export-skill":
            result = _export_skill(args.target, host=args.host)
        else:
            result = invoke_tool(args.tool, _read_input(args.input))
    except (AgentToolError, OSError, json.JSONDecodeError, ValueError) as exc:
        result = _base_envelope(command)
        _set_status(result, "invalid_input")
        result["findings"] = [_finding("tool_error", "blocked", str(exc))]
        result["next_actions"] = [{"tool": "discover", "reason": "Inspect the structured contract."}]
    if result.get("status") != "passed":
        for finding in result.get("findings", []):
            if isinstance(finding, Mapping) and finding.get("message"):
                print(f"dotmatch agent: {finding['message']}", file=sys.stderr)
    print(json.dumps(result, indent=2, sort_keys=True))
    return int(result.get("exit_code", 2))


__all__ = ["invoke_tool", "list_tools"]
