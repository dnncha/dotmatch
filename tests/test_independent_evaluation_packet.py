"""Regression tests for prospectively locked independent-evaluation packets."""

from __future__ import annotations

import copy
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/independent_evaluation_packet.py"
MODULE_SCRIPT = ROOT / "python/dotmatch/evaluation_packet.py"
SPEC = importlib.util.spec_from_file_location(
    "dotmatch_evaluation_packet", MODULE_SCRIPT
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("could not load packaged independent-evaluation implementation")
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def protocol() -> dict[str, object]:
    value = MODULE.template()
    value["evaluation_id"] = "external-lab-crispr-001"
    body = value["protocol"]
    body["protocol_owner"] = "independent evaluator"
    body["workflow"] = {
        "class": "pooled CRISPR screen guide counting",
        "dataset_scope": "one plasmid control and two treatment replicates",
        "private_data_handling": "input files remain in the evaluator-controlled workspace",
    }
    body["inputs"] = {
        "library": {"name": "library-v1", "revision": "2026-09", "sha256": "a" * 64},
        "samples": [
            {"id": "control", "role": "control", "read_mate": "R1", "sha256": "b" * 64},
            {
                "id": "treatment-a",
                "role": "treatment",
                "read_mate": "R1",
                "sha256": "c" * 64,
            },
        ],
    }
    body["dotmatch"] = {
        "version": "0.6.1",
        "install_route": "PyPI wheel with recorded artifact SHA-256",
        "command": "dotmatch crispr-count --library library.csv --samples samples.csv --out result",
        "settings": {
            "metric": "hamming",
            "threshold": 0,
            "read_start": 23,
            "read_length": 20,
        },
    }
    body["comparator"] = {
        "name": "established laboratory workflow",
        "version": "recorded immutable version",
        "install_route": "existing evaluator-controlled environment",
        "command": "lab-counter --manifest samples.csv --output baseline.tsv",
        "settings": {"assignment": "exact", "read_start": 23, "read_length": 20},
    }
    body["primary_endpoints"] = [
        {
            "id": "count_agreement",
            "metric": "identifier-keyed exact guide counts",
            "expected": "all guide and sample cells match",
            "failure_criterion": "any missing guide, annotation mismatch, or count mismatch",
        },
        {
            "id": "downstream_fdr_calls",
            "metric": "identifier-keyed MAGeCK score, p-value, FDR and hit-set agreement",
            "expected": "non-rank statistics and FDR <= 0.05 hit sets match",
            "failure_criterion": "any non-rank statistic or FDR hit-set difference",
        },
    ]
    body["tie_policy"] = {
        "ordinal_ranks": "compare_exact_unless_identical_statistical_tie",
        "top_n": "report_boundary_ties",
    }
    body["downstream_endpoint"] = {
        "planned": True,
        "tool": "MAGeCK",
        "version": "0.5.9.5",
        "command": "mageck test -k counts.tsv -t treatment-a -c control -n result",
        "endpoint": "score, p-value, FDR and FDR <= 0.05 hit sets",
    }
    body["failure_criteria"] = [
        "stop if input identity or sample roles cannot be established",
        "fail on any primary endpoint failure without post hoc threshold changes",
    ]
    body["scope_limits"] = [
        "software-output agreement is not biological correctness",
        "one evaluator-controlled workflow is not general adoption evidence",
    ]
    return value


def results(locked: dict[str, object]) -> dict[str, object]:
    return {
        "schema_version": 1,
        "evaluation_id": locked["evaluation_id"],
        "protocol_sha256": locked["protocol_sha256"],
        "completed_at": "2026-09-27T11:30:00Z",
        "protocol_record_reference": "private-vcs:commit-0123456789abcdef",
        "observations": [
            {
                "endpoint_id": "count_agreement",
                "observed": "all cells exact",
                "status": "pass",
            },
            {
                "endpoint_id": "downstream_fdr_calls",
                "observed": "scores and hit sets exact; tied ranks permuted",
                "status": "pass",
            },
        ],
        "operational_metrics": [
            {"metric": "wall_time_seconds", "dotmatch": 2.0, "comparator": 15.0}
        ],
        "conclusion": "outputs agreed under the locked tie semantics",
        "limitations": ["no biological truth source", "no claim beyond this workflow"],
        "public_use_permission": "no public organization or dataset name approved",
    }


def expect_error(callable_object, phrase: str) -> None:
    try:
        callable_object()
    except MODULE.PacketError as exc:
        assert phrase in str(exc), str(exc)
    else:
        raise AssertionError(f"expected PacketError containing {phrase!r}")


def main() -> None:
    source = protocol()
    locked = MODULE.lock_protocol(source, "2026-09-27T10:30:00Z")
    assert locked["status"] == "protocol_locked"
    assert locked["protocol_sha256"] == MODULE.protocol_sha256(source)
    MODULE.verify_locked_packet(locked)

    repeated = MODULE.lock_protocol(copy.deepcopy(source), "2026-09-28T10:30:00Z")
    assert repeated["protocol_sha256"] == locked["protocol_sha256"]

    completed = MODULE.complete_packet(locked, results(locked))
    assert completed["status"] == "completed"
    MODULE.verify_locked_packet(completed)

    tampered = copy.deepcopy(locked)
    tampered["protocol"]["dotmatch"]["settings"]["threshold"] = 1
    expect_error(
        lambda: MODULE.verify_locked_packet(tampered), "protocol hash mismatch"
    )

    placeholder = protocol()
    placeholder["protocol"]["comparator"]["version"] = "REPLACE_ME"
    expect_error(
        lambda: MODULE.lock_protocol(placeholder, "2026-09-27T10:30:00Z"), "placeholder"
    )

    missing_identity = protocol()
    missing_identity["protocol"]["inputs"]["samples"][0]["sha256"] = "unknown"
    expect_error(
        lambda: MODULE.lock_protocol(missing_identity, "2026-09-27T10:30:00Z"),
        "SHA-256",
    )

    unsafe_rank = protocol()
    unsafe_rank["protocol"]["tie_policy"]["ordinal_ranks"] = "exact_ordinal_rank"
    expect_error(
        lambda: MODULE.lock_protocol(unsafe_rank, "2026-09-27T10:30:00Z"), "unsupported"
    )

    unsafe_top_n = protocol()
    unsafe_top_n["protocol"]["tie_policy"]["top_n"] = "primary_without_tie_check"
    expect_error(
        lambda: MODULE.lock_protocol(unsafe_top_n, "2026-09-27T10:30:00Z"),
        "unsupported",
    )

    wrong_results = results(locked)
    wrong_results["evaluation_id"] = "other-evaluation"
    expect_error(lambda: MODULE.complete_packet(locked, wrong_results), "evaluation_id")

    incomplete_results = results(locked)
    incomplete_results["observations"] = incomplete_results["observations"][:1]
    expect_error(
        lambda: MODULE.complete_packet(locked, incomplete_results),
        "missing primary endpoints",
    )

    premature = copy.deepcopy(source)
    premature["results"] = {"conclusion": "known before lock"}
    expect_error(
        lambda: MODULE.lock_protocol(premature, "2026-09-27T10:30:00Z"),
        "unsupported fields",
    )

    with tempfile.TemporaryDirectory(prefix="dotmatch-evaluation-packet-") as raw_tmp:
        tmp = Path(raw_tmp)
        template_path = tmp / "template.json"
        subprocess.run(
            [
                sys.executable,
                str(MODULE_SCRIPT),
                "template",
                str(template_path),
            ],
            check=True,
            cwd=tmp,
        )
        assert template_path.exists()

        protocol_path = tmp / "protocol.json"
        protocol_path.write_text(MODULE.json.dumps(source), encoding="utf-8")
        locked_path = tmp / "locked.json"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "lock",
                str(protocol_path),
                str(locked_path),
                "--locked-at",
                "2026-09-27T10:30:00Z",
            ],
            check=True,
        )
        subprocess.run(
            [sys.executable, str(SCRIPT), "verify", str(locked_path)], check=True
        )

        tampered_path = tmp / "tampered.json"
        tampered_path.write_text(
            MODULE.json.dumps(
                {**MODULE.read_json(locked_path), "protocol_sha256": "0" * 64}
            ),
            encoding="utf-8",
        )
        rejected = subprocess.run(
            [
                sys.executable,
                str(MODULE_SCRIPT),
                "verify",
                str(tampered_path),
            ],
            cwd=tmp,
            capture_output=True,
            text=True,
            check=False,
        )
        assert rejected.returncode == 2
        assert "protocol hash mismatch" in rejected.stderr

    print("independent evaluation packet: PASS")


if __name__ == "__main__":
    main()
