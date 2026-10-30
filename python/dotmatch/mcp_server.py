from __future__ import annotations

import contextlib
import json
import os
import subprocess
import sys
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import IO, Any

from . import __version__
from .agent_tools import _package_json, invoke_tool, list_tools

SERVER_NAME = "dotmatch"
SUPPORTED_PROTOCOL_VERSIONS = ("2025-11-25", "2025-06-18", "2025-03-26", "2024-11-05")
FAILED_ENVELOPE_STATUSES = {"invalid_input", "blocked", "interrupted"}
OUTPUT_TAIL_CHARS = 4000

SERVER_INSTRUCTIONS = (
    "DotMatch assigns FASTQ reads to a finite list of known short DNA targets "
    "(CRISPR guides, inline barcodes, feature barcodes, primers, panel targets) using one fixed read window "
    "with Hamming k<=3 or Levenshtein k<=2 correction. Start with route_task when unsure which tool fits. "
    "For a known window use audit_targets before any k>=2 run, then count_targets or demux_barcodes. "
    "For CRISPR screens or Perturb-seq direct guide capture where the window is unknown, use the reviewed "
    "assay flow: discover, prepare_assay, inspect_assay, run_assay, review_assay, handoff_assay, following each "
    "result's next_actions. Coordinates are zero-based. Count only unique assignments; ambiguous, none and invalid "
    "reads stay separate. Everything runs locally; no data is uploaded. DotMatch is not a genome aligner, "
    "UMI or cell quantifier, variant caller, or screen-statistics package."
)

ASSAY_TOOL_DESCRIPTIONS = {
    "discover": (
        "Return the installed DotMatch assay contract: supported intents (CRISPR guide counting, Perturb-seq "
        "direct guide capture), tool schemas, safety policy and capability manifest. Read-only."
    ),
    "prepare_assay": (
        "Infer the guide window from a target table and a directory of FASTQ files and scaffold an immutable, "
        "reviewable AssaySpec project in an empty output directory. Use for CRISPR screen guide counting or "
        "Perturb-seq guide capture when the read window is not known in advance."
    ),
    "inspect_assay": "Preflight an AssaySpec: validate inputs, disk space and target-library risks without counting.",
    "run_assay": (
        "Run an AssaySpec end to end, producing raw counts (MAGeCK-compatible for CRISPR) and reliability evidence. "
        "Applies at most the allow-listed, evidence-backed spec revisions; never edits targets or counts ambiguity."
    ),
    "review_assay": "Summarise reliability findings, artifact hashes and next steps for a completed AssaySpec run.",
    "handoff_assay": "Package a raw-data-free, hash-verified review bundle for a completed AssaySpec run.",
}

ROUTE_TASK_SCHEMA: dict[str, Any] = {
    "type": "object",
    "additionalProperties": False,
    "properties": {
        "query": {
            "type": "string",
            "description": "Optional free-text task, e.g. 'count sgRNAs in FASTQ' or 'demultiplex inline barcodes'.",
        }
    },
}

COUNT_SCHEMA: dict[str, Any] = {
    "type": "object",
    "additionalProperties": False,
    "required": ["targets", "reads", "target_start", "target_length", "out"],
    "properties": {
        "targets": {"type": "string", "description": "CSV or TSV target table with id and sequence columns."},
        "reads": {
            "type": "array",
            "items": {"type": "string"},
            "minItems": 1,
            "description": "FASTQ or FASTQ.GZ paths, one per sample.",
        },
        "sample_labels": {
            "type": "array",
            "items": {"type": "string"},
            "description": "Sample names in the same order as reads; defaults to file stems.",
        },
        "target_start": {"type": "integer", "minimum": 0, "description": "Zero-based start of the target window."},
        "target_length": {"type": "integer", "minimum": 1, "description": "Length of the target window."},
        "k": {"type": "integer", "minimum": 0, "maximum": 3, "default": 0, "description": "Correction radius."},
        "metric": {"enum": ["hamming", "levenshtein"], "default": "hamming"},
        "ambiguity_policy": {"enum": ["radius", "best"], "default": "radius"},
        "format": {"enum": ["dotmatch", "mageck"], "default": "dotmatch"},
        "out": {"type": "string", "description": "New count matrix path; must not exist."},
    },
}

DEMUX_SCHEMA: dict[str, Any] = {
    "type": "object",
    "additionalProperties": False,
    "required": ["barcodes", "reads", "barcode_start", "barcode_length", "out_dir"],
    "properties": {
        "barcodes": {"type": "string", "description": "CSV or TSV barcode table with id and sequence columns."},
        "reads": {"type": "string", "description": "Pooled FASTQ or FASTQ.GZ path."},
        "barcode_start": {"type": "integer", "minimum": 0, "description": "Zero-based start of the inline barcode."},
        "barcode_length": {
            "oneOf": [{"type": "integer", "minimum": 1}, {"const": "auto"}],
            "description": "Barcode length, or 'auto' to infer from the table.",
        },
        "k": {"type": "integer", "minimum": 0, "maximum": 2, "default": 0},
        "metric": {"enum": ["hamming", "levenshtein"], "default": "hamming"},
        "ambiguity_policy": {"enum": ["radius", "best"], "default": "radius"},
        "out_dir": {"type": "string", "description": "New or empty directory for split FASTQs."},
    },
}

AUDIT_SCHEMA: dict[str, Any] = {
    "type": "object",
    "additionalProperties": False,
    "required": ["targets", "out_dir"],
    "properties": {
        "targets": {"type": "string", "description": "CSV or TSV target table with id and sequence columns."},
        "k": {"type": "integer", "minimum": 0, "maximum": 3, "default": 3},
        "audit_mode": {"enum": ["auto", "exact", "fast"], "default": "exact"},
        "out_dir": {"type": "string", "description": "New or empty directory for audit tables."},
    },
}

DISTANCE_SCHEMA: dict[str, Any] = {
    "type": "object",
    "additionalProperties": False,
    "required": ["a", "b"],
    "properties": {
        "a": {"type": "string", "minLength": 1, "maxLength": 4096},
        "b": {"type": "string", "minLength": 1, "maxLength": 4096},
    },
}

RESULT_SCHEMA: dict[str, Any] = {
    "type": "object",
    "required": ["status", "exit_code"],
    "properties": {"status": {"type": "string"}, "exit_code": {"type": "integer"}},
}


class McpToolInputError(ValueError):
    pass


def _require_keys(arguments: Mapping[str, Any], schema: Mapping[str, Any]) -> None:
    allowed = set(schema["properties"])
    unknown = sorted(set(arguments) - allowed)
    if unknown:
        raise McpToolInputError("unknown arguments: " + ", ".join(unknown))
    missing = [key for key in schema.get("required", []) if key not in arguments]
    if missing:
        raise McpToolInputError("missing required arguments: " + ", ".join(missing))


def _existing_file(value: object, label: str) -> str:
    path = Path(str(value)).expanduser()
    if not path.is_file():
        raise McpToolInputError(f"{label} must be an existing file: {path}")
    return str(path)


def _fresh_file(value: object, label: str) -> str:
    path = Path(str(value)).expanduser()
    if path.exists():
        raise McpToolInputError(f"{label} already exists; choose a new path: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return str(path)


def _empty_dir(value: object, label: str) -> str:
    path = Path(str(value)).expanduser()
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise McpToolInputError(f"{label} must be a new or empty directory: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return str(path)


def _integer(arguments: Mapping[str, Any], key: str, default: int, low: int, high: int | None = None) -> int:
    value = arguments.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int) or value < low or (high is not None and value > high):
        bound = f"{low}..{high}" if high is not None else f">= {low}"
        raise McpToolInputError(f"{key} must be an integer {bound}")
    return value


def _choice(arguments: Mapping[str, Any], key: str, schema: Mapping[str, Any]) -> str:
    spec = schema["properties"][key]
    value = arguments.get(key, spec.get("default"))
    if value not in spec["enum"]:
        raise McpToolInputError(f"{key} must be one of: " + ", ".join(spec["enum"]))
    return str(value)


def _tail(text: str) -> str:
    return text if len(text) <= OUTPUT_TAIL_CHARS else "..." + text[-OUTPUT_TAIL_CHARS:]


def _read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _run_cli(argv: Sequence[str]) -> dict[str, Any]:
    command = [sys.executable, "-m", "dotmatch", *argv]
    completed = subprocess.run(
        command, check=False, text=True, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    return {
        "status": "passed" if completed.returncode == 0 else "failed",
        "exit_code": completed.returncode,
        "command": ["dotmatch", *argv],
        "stdout": _tail(completed.stdout),
        "stderr": _tail(completed.stderr),
    }


def _route_task(arguments: Mapping[str, Any]) -> dict[str, Any]:
    _require_keys(arguments, ROUTE_TASK_SCHEMA)
    manifest = _package_json("agent-capabilities.json")
    query = str(arguments.get("query", "")).strip().lower()
    intents = manifest.get("intents", [])
    if query and isinstance(intents, list):
        words = [word for word in query.split() if len(word) > 2]
        scored = []
        for intent in intents:
            haystack = json.dumps(intent).lower()
            score = sum(haystack.count(word) for word in words)
            if score:
                scored.append((score, intent))
        scored.sort(key=lambda item: -item[0])
        matched = [intent for _, intent in scored[:3]]
    else:
        matched = intents
    return {
        "status": "passed",
        "exit_code": 0,
        "query": query,
        "matched_intents": matched,
        "scope": manifest.get("scope"),
        "error_recovery": manifest.get("error_recovery"),
        "out_of_scope_note": (
            "If no intent matches, DotMatch is probably the wrong tool: it only assigns reads to a known target list."
            if query and not matched
            else ""
        ),
    }


def _count_targets(arguments: Mapping[str, Any]) -> dict[str, Any]:
    _require_keys(arguments, COUNT_SCHEMA)
    reads_value = arguments["reads"]
    if not isinstance(reads_value, list) or not reads_value:
        raise McpToolInputError("reads must be a non-empty array of FASTQ paths")
    reads = [_existing_file(item, "reads") for item in reads_value]
    labels_value = arguments.get("sample_labels") or [Path(item).name.split(".")[0] for item in reads]
    if not isinstance(labels_value, list) or len(labels_value) != len(reads):
        raise McpToolInputError("sample_labels must have one entry per reads path")
    labels = [str(label) for label in labels_value]
    if any("," in label or not label for label in labels):
        raise McpToolInputError("sample_labels must be non-empty and contain no commas")
    metric = _choice(arguments, "metric", COUNT_SCHEMA)
    k = _integer(arguments, "k", 0, 0, 2 if metric == "levenshtein" else 3)
    out = _fresh_file(arguments["out"], "out")
    summary = out + ".summary.json"
    argv = ["count", "--targets", _existing_file(arguments["targets"], "targets")]
    for path in reads:
        argv += ["--reads", path]
    argv += [
        "--sample-label", ",".join(labels),
        "--target-start", str(_integer(arguments, "target_start", 0, 0)),
        "--target-length", str(_integer(arguments, "target_length", 0, 1)),
        "--k", str(k),
        "--metric", metric,
        "--ambiguity-policy", _choice(arguments, "ambiguity_policy", COUNT_SCHEMA),
        "--format", _choice(arguments, "format", COUNT_SCHEMA),
        "--no-progress",
        "--out", out,
        "--summary", summary,
    ]
    result = _run_cli(argv)
    result["outputs"] = {"counts": out, "summary": summary}
    result["summary"] = _read_json(Path(summary))
    if metric == "hamming" and k >= 2:
        result["warning"] = f"Confirm audit_targets reports safe_at_hamming_k{k}=true before trusting these counts."
    return result


def _demux_barcodes(arguments: Mapping[str, Any]) -> dict[str, Any]:
    _require_keys(arguments, DEMUX_SCHEMA)
    length = arguments["barcode_length"]
    if length != "auto":
        length = _integer(arguments, "barcode_length", 0, 1)
    out_dir = _empty_dir(arguments["out_dir"], "out_dir")
    summary = str(Path(out_dir) / "demux_summary.json")
    argv = [
        "demux",
        "--barcodes", _existing_file(arguments["barcodes"], "barcodes"),
        "--reads", _existing_file(arguments["reads"], "reads"),
        "--barcode-start", str(_integer(arguments, "barcode_start", 0, 0)),
        "--barcode-length", str(length),
        "--k", str(_integer(arguments, "k", 0, 0, 2)),
        "--metric", _choice(arguments, "metric", DEMUX_SCHEMA),
        "--ambiguity-policy", _choice(arguments, "ambiguity_policy", DEMUX_SCHEMA),
        "--out-dir", out_dir,
        "--summary", summary,
    ]
    result = _run_cli(argv)
    result["outputs"] = {"out_dir": out_dir, "files": sorted(p.name for p in Path(out_dir).iterdir())}
    result["summary"] = _read_json(Path(summary))
    return result


def _audit_targets(arguments: Mapping[str, Any]) -> dict[str, Any]:
    _require_keys(arguments, AUDIT_SCHEMA)
    out_dir = _empty_dir(arguments["out_dir"], "out_dir")
    argv = [
        "audit",
        "--targets", _existing_file(arguments["targets"], "targets"),
        "--k", str(_integer(arguments, "k", 3, 0, 3)),
        "--audit-mode", _choice(arguments, "audit_mode", AUDIT_SCHEMA),
        "--out-dir", out_dir,
    ]
    result = _run_cli(argv)
    result["outputs"] = {"out_dir": out_dir, "files": sorted(p.name for p in Path(out_dir).iterdir())}
    result["summary"] = _read_json(Path(out_dir) / "audit_summary.json")
    return result


def _edit_distance(arguments: Mapping[str, Any]) -> dict[str, Any]:
    _require_keys(arguments, DISTANCE_SCHEMA)
    from .core import distance

    a, b = str(arguments["a"]), str(arguments["b"])
    return {"status": "passed", "exit_code": 0, "a": a, "b": b, "levenshtein_distance": int(distance(a, b))}


DIRECT_TOOLS: dict[str, tuple[str, dict[str, Any], Callable[[Mapping[str, Any]], dict[str, Any]], bool]] = {
    "route_task": (
        "Decide whether DotMatch fits a sequencing task and which tool to call. Returns matching intents with "
        "commands, inputs, outputs and limitations from the installed capability manifest. Read-only.",
        ROUTE_TASK_SCHEMA,
        _route_task,
        True,
    ),
    "count_targets": (
        "Count FASTQ reads per known target (CRISPR sgRNAs, barcodes, primers, feature tags, panel targets) "
        "from one fixed zero-based read window with Hamming k<=3 or Levenshtein k<=2 correction. Writes a count "
        "matrix (optionally MAGeCK format) plus a JSON summary with assignment rates. Run audit_targets first for k>=2.",
        COUNT_SCHEMA,
        _count_targets,
        False,
    ),
    "demux_barcodes": (
        "Demultiplex a pooled FASTQ by inline barcodes at a fixed read position into per-barcode FASTQ files, "
        "with unassigned and ambiguous reads kept separate. Writes a QC summary JSON.",
        DEMUX_SCHEMA,
        _demux_barcodes,
        False,
    ),
    "audit_targets": (
        "Audit a guide, barcode or target library for duplicates and near-neighbour pairs that make correction "
        "ambiguous. Reports min Hamming distance and safe_at_hamming_k2/k3. Run before choosing a correction radius.",
        AUDIT_SCHEMA,
        _audit_targets,
        False,
    ),
    "edit_distance": (
        "Compute the global Levenshtein edit distance between two short DNA sequences. Read-only.",
        DISTANCE_SCHEMA,
        _edit_distance,
        True,
    ),
}


def _annotations(title: str, read_only: bool) -> dict[str, Any]:
    return {
        "title": title,
        "readOnlyHint": read_only,
        "destructiveHint": False,
        "idempotentHint": read_only,
        "openWorldHint": False,
    }


def tool_definitions() -> list[dict[str, Any]]:
    tools = []
    for name, (description, schema, _handler, read_only) in DIRECT_TOOLS.items():
        tools.append(
            {
                "name": name,
                "title": name.replace("_", " ").capitalize(),
                "description": description,
                "inputSchema": schema,
                "outputSchema": RESULT_SCHEMA,
                "annotations": _annotations(name.replace("_", " ").capitalize(), read_only),
            }
        )
    for item in list_tools()["tools"]:
        name = item["name"]
        read_only = name in {"discover", "inspect_assay", "review_assay"}
        tools.append(
            {
                "name": name,
                "title": name.replace("_", " ").capitalize(),
                "description": ASSAY_TOOL_DESCRIPTIONS.get(name, item.get("purpose", "")),
                "inputSchema": item["input_schema"],
                "outputSchema": RESULT_SCHEMA,
                "annotations": _annotations(name.replace("_", " ").capitalize(), read_only),
            }
        )
    return tools


def call_tool(name: str, arguments: Mapping[str, Any]) -> dict[str, Any]:
    if name in DIRECT_TOOLS:
        try:
            result = DIRECT_TOOLS[name][2](arguments)
        except (McpToolInputError, OSError, ValueError) as exc:
            result = {"status": "invalid_input", "exit_code": 2, "error": str(exc)}
    elif name in ASSAY_TOOL_DESCRIPTIONS:
        result = invoke_tool(name, arguments)
    else:
        raise KeyError(name)
    is_error = result.get("status") in FAILED_ENVELOPE_STATUSES or result.get("status") == "failed"
    return {
        "content": [{"type": "text", "text": json.dumps(result, indent=2, sort_keys=True)}],
        "structuredContent": result,
        "isError": is_error,
    }


RESOURCES = {
    "dotmatch://capabilities": ("agent-capabilities.json", "DotMatch capability manifest"),
    "dotmatch://tool-contract": ("agent-tools.json", "DotMatch assay tool contract"),
}


def _resources_list() -> dict[str, Any]:
    return {
        "resources": [
            {"uri": uri, "name": filename, "title": title, "mimeType": "application/json"}
            for uri, (filename, title) in RESOURCES.items()
        ]
    }


def _resources_read(params: Mapping[str, Any]) -> dict[str, Any]:
    uri = str(params.get("uri", ""))
    if uri not in RESOURCES:
        raise KeyError(uri)
    text = json.dumps(_package_json(RESOURCES[uri][0]), indent=2, sort_keys=True)
    return {"contents": [{"uri": uri, "mimeType": "application/json", "text": text}]}


def _initialize(params: Mapping[str, Any]) -> dict[str, Any]:
    requested = params.get("protocolVersion")
    version = requested if requested in SUPPORTED_PROTOCOL_VERSIONS else SUPPORTED_PROTOCOL_VERSIONS[0]
    return {
        "protocolVersion": version,
        "capabilities": {"tools": {"listChanged": False}, "resources": {"listChanged": False}},
        "serverInfo": {"name": SERVER_NAME, "title": "DotMatch", "version": __version__},
        "instructions": SERVER_INSTRUCTIONS,
    }


def _error(message_id: Any, code: int, message: str) -> dict[str, Any]:
    return {"jsonrpc": "2.0", "id": message_id, "error": {"code": code, "message": message}}


def handle_message(message: Any) -> dict[str, Any] | None:
    if not isinstance(message, dict) or message.get("jsonrpc") != "2.0":
        return _error(None, -32600, "invalid JSON-RPC 2.0 request")
    method = message.get("method")
    message_id = message.get("id")
    if "id" not in message:
        return None
    params = message.get("params") or {}
    if not isinstance(params, dict):
        return _error(message_id, -32602, "params must be an object")
    try:
        if method == "initialize":
            result = _initialize(params)
        elif method == "ping":
            result = {}
        elif method == "tools/list":
            result = {"tools": tool_definitions()}
        elif method == "tools/call":
            arguments = params.get("arguments") or {}
            if not isinstance(arguments, dict):
                return _error(message_id, -32602, "arguments must be an object")
            try:
                result = call_tool(str(params.get("name", "")), arguments)
            except KeyError:
                return _error(message_id, -32602, f"unknown tool: {params.get('name')}")
        elif method == "resources/list":
            result = _resources_list()
        elif method == "resources/read":
            try:
                result = _resources_read(params)
            except KeyError:
                return _error(message_id, -32002, f"resource not found: {params.get('uri')}")
        elif method in {"prompts/list", "resources/templates/list"}:
            result = {"prompts": []} if method == "prompts/list" else {"resourceTemplates": []}
        else:
            return _error(message_id, -32601, f"method not found: {method}")
    except Exception as exc:
        return _error(message_id, -32603, f"internal error: {exc}")
    return {"jsonrpc": "2.0", "id": message_id, "result": result}


def serve(stdin: IO[str], protocol_out: IO[str]) -> None:
    for line in stdin:
        if not line.strip():
            continue
        try:
            message = json.loads(line)
        except ValueError:
            response: dict[str, Any] | None = _error(None, -32700, "parse error")
        else:
            with contextlib.redirect_stdout(sys.stderr):
                response = handle_message(message)
        if response is not None:
            protocol_out.write(json.dumps(response, separators=(",", ":")) + "\n")
            protocol_out.flush()


def _protocol_stream() -> IO[str]:
    sys.stdout.flush()
    protocol_fd = os.dup(1)
    os.dup2(2, 1)
    return os.fdopen(protocol_fd, "w", encoding="utf-8", buffering=1)


def main(argv: Sequence[str] = ()) -> int:
    if any(arg in {"-h", "--help"} for arg in argv):
        print(
            "usage: dotmatch mcp\n\n"
            "Serve DotMatch as a local Model Context Protocol server over stdio.\n"
            "Add it to an MCP client, for example:\n"
            "  claude mcp add dotmatch -- dotmatch mcp\n"
            "  codex mcp add dotmatch -- dotmatch mcp\n"
            '  {"mcpServers": {"dotmatch": {"command": "dotmatch", "args": ["mcp"]}}}\n'
            "Or without installing: uvx --from git+https://github.com/dnncha/dotmatch dotmatch mcp"
        )
        return 0
    if argv:
        print(f"dotmatch mcp: unexpected arguments: {' '.join(argv)}", file=sys.stderr)
        return 2
    serve(sys.stdin, _protocol_stream())
    return 0


__all__ = ["call_tool", "handle_message", "main", "serve", "tool_definitions"]
