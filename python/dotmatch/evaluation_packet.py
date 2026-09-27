"""Create and verify prospectively locked independent-evaluation packets."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ALLOWED_RANK_POLICIES = {
    "not_applicable",
    "compare_exact_unless_identical_statistical_tie",
}
ALLOWED_TOP_N_POLICIES = {"not_primary", "report_boundary_ties"}
HEX = set("0123456789abcdef")


class PacketError(ValueError):
    """Raised when an evaluation packet is incomplete or inconsistent."""


def _object(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise PacketError(f"{label} must be an object")
    return value


def _list(value: Any, label: str) -> list[Any]:
    if not isinstance(value, list) or not value:
        raise PacketError(f"{label} must be a non-empty array")
    return value


def _text(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise PacketError(f"{label} must be a non-empty string")
    if "REPLACE_ME" in value.upper():
        raise PacketError(f"{label} still contains a template placeholder")
    return value.strip()


def _sha256(value: Any, label: str) -> str:
    digest = _text(value, label).lower()
    if len(digest) != 64 or any(character not in HEX for character in digest):
        raise PacketError(f"{label} must be a 64-character lowercase SHA-256")
    return digest


def _exact_keys(document: dict[str, Any], expected: set[str], label: str) -> None:
    missing = expected - set(document)
    extra = set(document) - expected
    if missing:
        raise PacketError(f"{label} missing fields: {', '.join(sorted(missing))}")
    if extra:
        raise PacketError(f"{label} has unsupported fields: {', '.join(sorted(extra))}")


def canonical_protocol_payload(document: dict[str, Any]) -> bytes:
    payload = {
        "schema_version": document["schema_version"],
        "evaluation_id": document["evaluation_id"],
        "protocol": document["protocol"],
    }
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode()


def protocol_sha256(document: dict[str, Any]) -> str:
    return hashlib.sha256(canonical_protocol_payload(document)).hexdigest()


def validate_protocol_document(document: dict[str, Any]) -> None:
    _exact_keys(
        document, {"schema_version", "evaluation_id", "protocol"}, "protocol document"
    )
    if document["schema_version"] != 1:
        raise PacketError("schema_version must be 1")
    _text(document["evaluation_id"], "evaluation_id")
    protocol = _object(document["protocol"], "protocol")
    _exact_keys(
        protocol,
        {
            "protocol_locked_before_results",
            "protocol_owner",
            "workflow",
            "inputs",
            "dotmatch",
            "comparator",
            "primary_endpoints",
            "tie_policy",
            "downstream_endpoint",
            "failure_criteria",
            "scope_limits",
        },
        "protocol",
    )
    if protocol["protocol_locked_before_results"] is not True:
        raise PacketError("protocol_locked_before_results must be true")
    _text(protocol["protocol_owner"], "protocol.protocol_owner")

    workflow = _object(protocol["workflow"], "protocol.workflow")
    _exact_keys(
        workflow,
        {"class", "dataset_scope", "private_data_handling"},
        "protocol.workflow",
    )
    for field in workflow:
        _text(workflow[field], f"protocol.workflow.{field}")

    inputs = _object(protocol["inputs"], "protocol.inputs")
    _exact_keys(inputs, {"library", "samples"}, "protocol.inputs")
    library = _object(inputs["library"], "protocol.inputs.library")
    _exact_keys(library, {"name", "revision", "sha256"}, "protocol.inputs.library")
    _text(library["name"], "protocol.inputs.library.name")
    _text(library["revision"], "protocol.inputs.library.revision")
    _sha256(library["sha256"], "protocol.inputs.library.sha256")
    samples = _list(inputs["samples"], "protocol.inputs.samples")
    sample_ids: set[str] = set()
    for index, raw_sample in enumerate(samples):
        sample = _object(raw_sample, f"protocol.inputs.samples[{index}]")
        _exact_keys(
            sample,
            {"id", "role", "read_mate", "sha256"},
            f"protocol.inputs.samples[{index}]",
        )
        sample_id = _text(sample["id"], f"protocol.inputs.samples[{index}].id")
        if sample_id in sample_ids:
            raise PacketError(f"duplicate sample id: {sample_id}")
        sample_ids.add(sample_id)
        _text(sample["role"], f"protocol.inputs.samples[{index}].role")
        _text(sample["read_mate"], f"protocol.inputs.samples[{index}].read_mate")
        _sha256(sample["sha256"], f"protocol.inputs.samples[{index}].sha256")

    for tool_name in ("dotmatch", "comparator"):
        tool = _object(protocol[tool_name], f"protocol.{tool_name}")
        required = {"version", "install_route", "command", "settings"}
        if tool_name == "comparator":
            required.add("name")
        _exact_keys(tool, required, f"protocol.{tool_name}")
        for field in required - {"settings"}:
            _text(tool[field], f"protocol.{tool_name}.{field}")
        if not isinstance(tool["settings"], dict):
            raise PacketError(f"protocol.{tool_name}.settings must be an object")

    endpoints = _list(protocol["primary_endpoints"], "protocol.primary_endpoints")
    endpoint_ids: set[str] = set()
    for index, raw_endpoint in enumerate(endpoints):
        endpoint = _object(raw_endpoint, f"protocol.primary_endpoints[{index}]")
        _exact_keys(
            endpoint,
            {"id", "metric", "expected", "failure_criterion"},
            f"protocol.primary_endpoints[{index}]",
        )
        endpoint_id = _text(endpoint["id"], f"protocol.primary_endpoints[{index}].id")
        if endpoint_id in endpoint_ids:
            raise PacketError(f"duplicate primary endpoint id: {endpoint_id}")
        endpoint_ids.add(endpoint_id)
        for field in ("metric", "expected", "failure_criterion"):
            _text(endpoint[field], f"protocol.primary_endpoints[{index}].{field}")

    tie_policy = _object(protocol["tie_policy"], "protocol.tie_policy")
    _exact_keys(tie_policy, {"ordinal_ranks", "top_n"}, "protocol.tie_policy")
    if tie_policy["ordinal_ranks"] not in ALLOWED_RANK_POLICIES:
        raise PacketError("protocol.tie_policy.ordinal_ranks is unsupported")
    if tie_policy["top_n"] not in ALLOWED_TOP_N_POLICIES:
        raise PacketError("protocol.tie_policy.top_n is unsupported")

    downstream = _object(
        protocol["downstream_endpoint"], "protocol.downstream_endpoint"
    )
    _exact_keys(
        downstream,
        {"planned", "tool", "version", "command", "endpoint"},
        "protocol.downstream_endpoint",
    )
    if not isinstance(downstream["planned"], bool):
        raise PacketError("protocol.downstream_endpoint.planned must be boolean")
    if downstream["planned"]:
        for field in ("tool", "version", "command", "endpoint"):
            _text(downstream[field], f"protocol.downstream_endpoint.{field}")
    elif any(
        downstream[field] not in (None, "")
        for field in ("tool", "version", "command", "endpoint")
    ):
        raise PacketError("unplanned downstream fields must be null or empty")

    for label in ("failure_criteria", "scope_limits"):
        values = _list(protocol[label], f"protocol.{label}")
        for index, value in enumerate(values):
            _text(value, f"protocol.{label}[{index}]")


def lock_protocol(document: dict[str, Any], locked_at: str) -> dict[str, Any]:
    validate_protocol_document(document)
    _text(locked_at, "locked_at")
    return {
        **document,
        "status": "protocol_locked",
        "locked_at": locked_at,
        "protocol_sha256": protocol_sha256(document),
    }


def verify_locked_packet(packet: dict[str, Any]) -> None:
    required = {
        "schema_version",
        "evaluation_id",
        "protocol",
        "status",
        "locked_at",
        "protocol_sha256",
    }
    if packet.get("status") == "completed":
        required.add("results")
    _exact_keys(packet, required, "evaluation packet")
    validate_protocol_document(
        {key: packet[key] for key in ("schema_version", "evaluation_id", "protocol")}
    )
    if packet["status"] not in {"protocol_locked", "completed"}:
        raise PacketError("status must be protocol_locked or completed")
    _text(packet["locked_at"], "locked_at")
    recorded_hash = _sha256(packet["protocol_sha256"], "protocol_sha256")
    actual_hash = protocol_sha256(packet)
    if recorded_hash != actual_hash:
        raise PacketError(
            f"protocol hash mismatch: expected {recorded_hash}, calculated {actual_hash}"
        )
    if packet["status"] == "completed":
        validate_results(packet["results"], packet)


def validate_results(results: Any, locked: dict[str, Any]) -> None:
    result = _object(results, "results")
    _exact_keys(
        result,
        {
            "schema_version",
            "evaluation_id",
            "protocol_sha256",
            "completed_at",
            "protocol_record_reference",
            "observations",
            "operational_metrics",
            "conclusion",
            "limitations",
            "public_use_permission",
        },
        "results",
    )
    if result["schema_version"] != 1:
        raise PacketError("results.schema_version must be 1")
    if result["evaluation_id"] != locked["evaluation_id"]:
        raise PacketError("results evaluation_id does not match the locked protocol")
    if result["protocol_sha256"] != locked["protocol_sha256"]:
        raise PacketError("results protocol_sha256 does not match the locked protocol")
    for field in (
        "completed_at",
        "protocol_record_reference",
        "conclusion",
        "public_use_permission",
    ):
        _text(result[field], f"results.{field}")
    observations = _list(result["observations"], "results.observations")
    endpoint_ids = {item["id"] for item in locked["protocol"]["primary_endpoints"]}
    observed_ids: set[str] = set()
    for index, raw_observation in enumerate(observations):
        observation = _object(raw_observation, f"results.observations[{index}]")
        _exact_keys(
            observation,
            {"endpoint_id", "observed", "status"},
            f"results.observations[{index}]",
        )
        endpoint_id = _text(
            observation["endpoint_id"], f"results.observations[{index}].endpoint_id"
        )
        if endpoint_id not in endpoint_ids:
            raise PacketError(f"result references unknown endpoint: {endpoint_id}")
        if endpoint_id in observed_ids:
            raise PacketError(f"duplicate result endpoint: {endpoint_id}")
        observed_ids.add(endpoint_id)
        _text(observation["observed"], f"results.observations[{index}].observed")
        if observation["status"] not in {"pass", "fail", "inconclusive", "not_run"}:
            raise PacketError(
                f"unsupported observation status: {observation['status']}"
            )
    missing = endpoint_ids - observed_ids
    if missing:
        raise PacketError(
            f"results missing primary endpoints: {', '.join(sorted(missing))}"
        )
    if not isinstance(result["operational_metrics"], list):
        raise PacketError("results.operational_metrics must be an array")
    limitations = _list(result["limitations"], "results.limitations")
    for index, limitation in enumerate(limitations):
        _text(limitation, f"results.limitations[{index}]")


def complete_packet(locked: dict[str, Any], results: dict[str, Any]) -> dict[str, Any]:
    verify_locked_packet(locked)
    if locked["status"] != "protocol_locked":
        raise PacketError("only a protocol_locked packet can be completed")
    validate_results(results, locked)
    return {**locked, "status": "completed", "results": results}


def template() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "evaluation_id": "REPLACE_ME",
        "protocol": {
            "protocol_locked_before_results": True,
            "protocol_owner": "REPLACE_ME",
            "workflow": {
                "class": "REPLACE_ME",
                "dataset_scope": "REPLACE_ME",
                "private_data_handling": "REPLACE_ME",
            },
            "inputs": {
                "library": {
                    "name": "REPLACE_ME",
                    "revision": "REPLACE_ME",
                    "sha256": "REPLACE_ME",
                },
                "samples": [
                    {
                        "id": "REPLACE_ME",
                        "role": "control",
                        "read_mate": "R1",
                        "sha256": "REPLACE_ME",
                    }
                ],
            },
            "dotmatch": {
                "version": "REPLACE_ME",
                "install_route": "REPLACE_ME",
                "command": "REPLACE_ME",
                "settings": {},
            },
            "comparator": {
                "name": "REPLACE_ME",
                "version": "REPLACE_ME",
                "install_route": "REPLACE_ME",
                "command": "REPLACE_ME",
                "settings": {},
            },
            "primary_endpoints": [
                {
                    "id": "count_agreement",
                    "metric": "REPLACE_ME",
                    "expected": "REPLACE_ME",
                    "failure_criterion": "REPLACE_ME",
                }
            ],
            "tie_policy": {
                "ordinal_ranks": "not_applicable",
                "top_n": "not_primary",
            },
            "downstream_endpoint": {
                "planned": False,
                "tool": None,
                "version": None,
                "command": None,
                "endpoint": None,
            },
            "failure_criteria": ["REPLACE_ME"],
            "scope_limits": ["REPLACE_ME"],
        },
    }


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PacketError(f"could not read {path}: {exc}") from exc
    return _object(value, str(path))


def write_json(path: Path, value: dict[str, Any]) -> None:
    if path.exists():
        raise PacketError(f"refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def now_utc() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    template_parser = subparsers.add_parser(
        "template", help="write an incomplete protocol template"
    )
    template_parser.add_argument("output", type=Path)
    lock_parser = subparsers.add_parser(
        "lock", help="validate and lock a protocol before outcomes are known"
    )
    lock_parser.add_argument("protocol", type=Path)
    lock_parser.add_argument("output", type=Path)
    lock_parser.add_argument("--locked-at", default=None)
    verify_parser = subparsers.add_parser(
        "verify", help="verify a locked or completed packet"
    )
    verify_parser.add_argument("packet", type=Path)
    complete_parser = subparsers.add_parser(
        "complete", help="attach outcomes to an unchanged locked protocol"
    )
    complete_parser.add_argument("locked_packet", type=Path)
    complete_parser.add_argument("results", type=Path)
    complete_parser.add_argument("output", type=Path)
    args = parser.parse_args(argv)

    try:
        if args.command == "template":
            write_json(args.output, template())
            print(args.output)
        elif args.command == "lock":
            locked = lock_protocol(
                read_json(args.protocol), args.locked_at or now_utc()
            )
            write_json(args.output, locked)
            print(locked["protocol_sha256"])
        elif args.command == "verify":
            packet = read_json(args.packet)
            verify_locked_packet(packet)
            print(f"{packet['status']}: {packet['protocol_sha256']}")
        else:
            completed = complete_packet(
                read_json(args.locked_packet), read_json(args.results)
            )
            write_json(args.output, completed)
            print(completed["protocol_sha256"])
    except PacketError as exc:
        print(f"dotmatch evaluation-packet: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
