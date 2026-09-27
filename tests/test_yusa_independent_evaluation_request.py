"""Validate the unsent Yusa independent-evaluation draft against evidence."""

from __future__ import annotations

import copy
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "python/dotmatch/evaluation_packet.py"
SPEC = importlib.util.spec_from_file_location("dotmatch_evaluation_packet", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("could not load evaluation packet implementation")
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    draft = load_json(
        ROOT / "examples/evaluation_packets/yusa-independent-protocol.draft.json"
    )
    counts = load_json(
        ROOT / "benchmarks/raw/mageck_yusa_full_implementation_comparison.json"
    )
    downstream = load_json(
        ROOT / "benchmarks/raw/mageck_yusa_full_downstream_comparison.json"
    )

    try:
        MODULE.lock_protocol(draft, "2026-09-27T12:00:00Z")
    except MODULE.PacketError as exc:
        assert "template placeholder" in str(exc)
    else:
        raise AssertionError("unsent draft must refuse locking until evaluator-owned")

    evidence_inputs = counts["inputs"]
    protocol_inputs = draft["protocol"]["inputs"]
    assert (
        protocol_inputs["library"]["sha256"]
        == evidence_inputs["yusa_library.csv"]["sha256"]
    )
    assert {
        sample["id"]: sample["sha256"] for sample in protocol_inputs["samples"]
    } == {
        "ERR376998": evidence_inputs["ERR376998.fastq.gz"]["sha256"],
        "ERR376999": evidence_inputs["ERR376999.fastq.gz"]["sha256"],
    }
    assert (
        draft["protocol"]["dotmatch"]["version"]
        == counts["dotmatch_version"]
        == "0.6.1"
    )
    assert (
        draft["protocol"]["comparator"]["version"]
        == counts["mageck_runtime"]["version"]
    )
    assert (
        counts["mageck_runtime"]["package_sha256"]
        in draft["protocol"]["comparator"]["install_route"]
    )
    assert draft["protocol"]["tie_policy"] == {
        "ordinal_ranks": "compare_exact_unless_identical_statistical_tie",
        "top_n": "report_boundary_ties",
    }
    assert downstream["protocol"]["settings"]["fdr_threshold"] == 0.05

    evaluator_owned = copy.deepcopy(draft)
    evaluator_owned["evaluation_id"] = "independent-yusa-example-001"
    evaluator_owned["protocol"]["protocol_owner"] = "independent example evaluator"
    locked = MODULE.lock_protocol(evaluator_owned, "2026-09-27T12:00:00Z")
    MODULE.verify_locked_packet(locked)
    assert locked["protocol_sha256"] == MODULE.protocol_sha256(evaluator_owned)

    print(f"yusa independent evaluation draft: PASS {locked['protocol_sha256']}")


if __name__ == "__main__":
    main()
