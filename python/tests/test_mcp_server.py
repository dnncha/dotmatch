from __future__ import annotations

import io
import json
import subprocess
import sys
from pathlib import Path

from dotmatch.mcp_server import call_tool, handle_message, serve, tool_definitions

ROOT = Path(__file__).resolve().parents[2]
TARGETS = ROOT / "demo-data" / "crispr_guides.tsv"
READS = ROOT / "demo-data" / "reads.fastq"


def _rpc(message_id: int, method: str, params: dict[str, object] | None = None) -> dict[str, object]:
    message: dict[str, object] = {"jsonrpc": "2.0", "id": message_id, "method": method}
    if params is not None:
        message["params"] = params
    return message


def test_lists_direct_and_assay_tools_with_schemas() -> None:
    names = [tool["name"] for tool in tool_definitions()]
    assert names[:5] == ["route_task", "count_targets", "demux_barcodes", "audit_targets", "edit_distance"]
    assert set(names[5:]) == {"discover", "prepare_assay", "inspect_assay", "run_assay", "review_assay", "handoff_assay"}
    for tool in tool_definitions():
        assert tool["description"]
        assert tool["inputSchema"]["type"] == "object"


def test_initialize_negotiates_protocol_and_ignores_notifications() -> None:
    response = handle_message(_rpc(1, "initialize", {"protocolVersion": "2025-06-18"}))
    assert response is not None
    result = response["result"]
    assert result["protocolVersion"] == "2025-06-18"
    assert result["serverInfo"]["name"] == "dotmatch"
    assert "known short DNA targets" in result["instructions"]
    assert handle_message({"jsonrpc": "2.0", "method": "notifications/initialized"}) is None
    unknown = handle_message(_rpc(2, "initialize", {"protocolVersion": "1999-01-01"}))
    assert unknown is not None and unknown["result"]["protocolVersion"] == "2025-11-25"


def test_unknown_tool_and_method_return_jsonrpc_errors() -> None:
    tool = handle_message(_rpc(1, "tools/call", {"name": "nope", "arguments": {}}))
    method = handle_message(_rpc(2, "sampling/createMessage"))
    assert tool is not None and tool["error"]["code"] == -32602
    assert method is not None and method["error"]["code"] == -32601


def test_count_audit_and_refusals(tmp_path: Path) -> None:
    audit = call_tool("audit_targets", {"targets": str(TARGETS), "out_dir": str(tmp_path / "audit")})
    assert audit["isError"] is False
    assert audit["structuredContent"]["summary"]["min_hamming_distance"] == 3

    out = tmp_path / "counts.tsv"
    arguments = {
        "targets": str(TARGETS),
        "reads": [str(READS)],
        "target_start": 0,
        "target_length": 4,
        "k": 1,
        "out": str(out),
    }
    counted = call_tool("count_targets", arguments)
    assert counted["isError"] is False, counted["structuredContent"]
    rows = {line.split("\t")[0]: line.split("\t")[-1] for line in out.read_text().splitlines()[1:]}
    assert rows["guide_1"] == "2"

    again = call_tool("count_targets", arguments)
    assert again["isError"] is True
    assert "already exists" in again["structuredContent"]["error"]

    unknown = call_tool("count_targets", {**arguments, "out": str(tmp_path / "x.tsv"), "shell": "rm -rf /"})
    assert unknown["structuredContent"]["status"] == "invalid_input"


def test_assay_tools_pass_through_contract_envelope() -> None:
    result = call_tool("discover", {"intent": "crispr-guide-counting"})
    assert result["isError"] is False
    assert result["structuredContent"]["result"]["local_only"] is True


def test_serve_writes_only_jsonrpc_lines() -> None:
    requests = [
        _rpc(1, "initialize", {"protocolVersion": "2025-06-18"}),
        {"jsonrpc": "2.0", "method": "notifications/initialized"},
        _rpc(2, "tools/call", {"name": "edit_distance", "arguments": {"a": "ACGT", "b": "AGGT"}}),
        _rpc(3, "resources/read", {"uri": "dotmatch://capabilities"}),
    ]
    stdin = io.StringIO("".join(json.dumps(item) + "\n" for item in requests) + "not json\n")
    stdout = io.StringIO()
    serve(stdin, stdout)
    responses = [json.loads(line) for line in stdout.getvalue().splitlines()]
    assert [item.get("id") for item in responses] == [1, 2, 3, None]
    assert responses[1]["result"]["structuredContent"]["levenshtein_distance"] == 1
    assert json.loads(responses[2]["result"]["contents"][0]["text"])["project"]["package"] == "dotmatch"
    assert responses[3]["error"]["code"] == -32700


def test_cli_subcommand_speaks_mcp_over_stdio() -> None:
    request = json.dumps(_rpc(1, "tools/list")) + "\n"
    completed = subprocess.run(
        [sys.executable, "-m", "dotmatch", "mcp"], input=request, capture_output=True, text=True, timeout=60, check=False
    )
    assert completed.returncode == 0, completed.stderr
    lines = completed.stdout.splitlines()
    assert len(lines) == 1
    assert len(json.loads(lines[0])["result"]["tools"]) == 11
