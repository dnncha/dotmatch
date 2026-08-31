#!/usr/bin/env python3
"""Validate the installed/source agent contract without executing assay data."""

from __future__ import annotations

import json
import os
import socket
import subprocess
import sys
from pathlib import Path

try:
    from jsonschema import Draft202012Validator
except ModuleNotFoundError:  # clean source check; release dev environments install the declared extra
    Draft202012Validator = None  # type: ignore[assignment,misc]


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from dotmatch.agent_tools import invoke_tool, list_tools  # noqa: E402


def main() -> int:
    failures: list[str] = []
    contract = json.loads((ROOT / "agent-tools.json").read_text(encoding="utf-8"))
    schema = json.loads((ROOT / "agent-tools.schema.json").read_text(encoding="utf-8"))
    if Draft202012Validator is not None:
        try:
            Draft202012Validator.check_schema(schema)
            Draft202012Validator(schema).validate(contract)
        except Exception as exc:
            failures.append(f"agent tool contract/schema validation failed: {exc}")
    else:
        if schema.get("$schema") != "https://json-schema.org/draft/2020-12/schema":
            failures.append("agent tool schema must declare JSON Schema draft 2020-12")
        if not isinstance(schema.get("$defs", {}).get("envelope"), dict):
            failures.append("agent tool schema must define the stable envelope")
    if list_tools() != contract:
        failures.append("installed source-package agent contract differs from the canonical contract")

    names = {item.get("name") for item in contract.get("tools", []) if isinstance(item, dict)}
    expected = {"discover", "prepare_assay", "inspect_assay", "run_assay", "review_assay", "handoff_assay"}
    if names != expected:
        failures.append(f"agent tool names differ: {sorted(names)}")

    original_connect = socket.socket.connect
    try:
        def refuse_connect(*_args: object, **_kwargs: object) -> None:
            raise AssertionError("outbound network connection attempted")

        socket.socket.connect = refuse_connect  # type: ignore[method-assign]
        envelope = invoke_tool("discover", {})
    finally:
        socket.socket.connect = original_connect  # type: ignore[method-assign]
    if Draft202012Validator is not None:
        envelope_schema = {"$ref": "#/$defs/envelope", "$defs": schema["$defs"]}
        envelope_validator = Draft202012Validator(envelope_schema)
        try:
            envelope_validator.validate(envelope)
        except Exception as exc:
            failures.append(f"discover envelope does not validate: {exc}")
        invalid_mapping = dict(envelope)
        invalid_mapping["exit_code"] = 2
        if not list(envelope_validator.iter_errors(invalid_mapping)):
            failures.append("agent tool schema must reject a passed status paired with exit code 2")
    else:
        required_envelope = set(schema["$defs"]["envelope"]["required"])
        missing_envelope = sorted(required_envelope - set(envelope))
        if missing_envelope:
            failures.append("discover envelope is missing schema fields: " + ", ".join(missing_envelope))
        if len(schema["$defs"]["envelope"].get("allOf", [])) < 4:
            failures.append("agent tool schema must bind statuses to exit codes")
    if envelope.get("status") != "passed" or envelope.get("exit_code") != 0:
        failures.append("discover must map passed to exit code 0")

    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "python")
    cli = subprocess.run(
        [sys.executable, "-m", "dotmatch.cli", "agent", "invoke", "discover", "--input", "-"],
        input='{"intent":"crispr-guide-counting"}\n',
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=ROOT,
        env=env,
        check=False,
    )
    try:
        cli_envelope = json.loads(cli.stdout)
    except json.JSONDecodeError as exc:
        failures.append(f"agent CLI stdout is not one JSON document: {exc}")
        cli_envelope = {}
    if cli.returncode != 0 or cli_envelope.get("status") != "passed":
        failures.append("agent CLI discover forward test failed")

    invalid = subprocess.run(
        [sys.executable, "-m", "dotmatch.cli", "agent", "invoke", "discover", "--input", "-"],
        input='{"command":"touch must-not-run"}\n',
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=ROOT,
        env=env,
        check=False,
    )
    try:
        invalid_envelope = json.loads(invalid.stdout)
    except json.JSONDecodeError:
        invalid_envelope = {}
    if invalid.returncode != 2 or invalid_envelope.get("status") != "invalid_input":
        failures.append("invalid structured input must map to status invalid_input and exit code 2")
    if "unknown input field(s): command" not in invalid.stderr:
        failures.append("agent CLI diagnostics must be written to stderr")
    if (ROOT / "must-not-run").exists():
        failures.append("shell-shaped structured input was executed")

    missing_args = subprocess.run(
        [sys.executable, "-m", "dotmatch.cli", "agent", "invoke", "discover"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=ROOT,
        env=env,
        check=False,
    )
    try:
        missing_args_envelope = json.loads(missing_args.stdout)
    except json.JSONDecodeError:
        missing_args_envelope = {}
    if missing_args.returncode != 2 or missing_args_envelope.get("status") != "invalid_input":
        failures.append("agent CLI argument errors must retain the JSON envelope and exit code 2")

    copies = {
        "agent-tools.json": ["docs/agent-tools.json", "public/agent-tools.json", "python/dotmatch/data/agent-tools.json"],
        "agent-tools.schema.json": ["docs/agent-tools.schema.json", "public/agent-tools.schema.json", "python/dotmatch/data/agent-tools.schema.json"],
    }
    for source, destinations in copies.items():
        expected_bytes = (ROOT / source).read_bytes()
        for destination in destinations:
            path = ROOT / destination
            if not path.is_file() or path.read_bytes() != expected_bytes:
                failures.append(f"stale agent tool copy: {destination}")

    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1
    print("Agent tool contract, envelope, CLI routing, no-network discover, and drift checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
