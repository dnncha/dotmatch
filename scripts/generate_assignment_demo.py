#!/usr/bin/env python3
"""Generate the website's synthetic example with the actual native matcher."""
from __future__ import annotations
import argparse
import json
import hashlib
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from dotmatch.core import Matcher, iter_fastq, status_name
from dotmatch.target_io import read_target_table
from dotmatch.sensitivity import compare_hamming_policies, MODES, STATES


def generate():
    fixture = ROOT / "examples/assignment_sensitivity"
    targets = read_target_table(fixture / "targets.tsv")
    reads = list(iter_fastq(fixture / "reads.fastq"))
    with Matcher([row.sequence for row in targets]) as matcher:
        results = compare_hamming_policies(
            matcher, [row.seq if len(row.seq) == 20 else None for row in reads]
        )
    outcomes = {mode: {state: 0 for state in STATES} for mode in MODES}
    records = []
    for read, result in zip(reads, results):
        calls = {}
        for mode in MODES:
            call = getattr(result, mode)
            state = status_name(call.status)
            outcomes[mode][state] += 1
            calls[mode] = {
                "status": state,
                "target_id": (
                    targets[call.target_index].target_id
                    if call.target_index >= 0
                    else None
                ),
            }
        candidates = [
            target.target_id
            for target in targets
            if len(read.seq) == len(target.sequence)
            and sum(a != b for a, b in zip(read.seq, target.sequence)) <= 1
        ]
        assert len(candidates) == result.candidates_within_one
        records.append(
            {
                "id": read.read_id,
                "sequence": read.seq,
                "calls": calls,
                "candidate_ids": candidates,
            }
        )
    return {
        "schema_version": "dotmatch.assignment-demo.v1",
        "kind": "synthetic_software_example",
        "metric": "hamming",
        "target_length": 20,
        "target_start": 0,
        "read_count": len(reads),
        "target_count": len(targets),
        "changed_reads": sum(row.changed for row in results),
        "targets": [{"id": row.target_id, "sequence": row.sequence} for row in targets],
        "records": records,
        "outcomes": outcomes,
        "inputs_sha256": {
            name: hashlib.sha256((fixture / name).read_bytes()).hexdigest()
            for name in ["targets.tsv", "reads.fastq"]
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    text = json.dumps(generate(), indent=2) + "\n"
    output = ROOT / "public/assignment-demo.json"
    if args.check:
        if not output.exists() or output.read_text() != text:
            raise SystemExit(
                "Website assignment demo differs from native results. Regenerate and review the fixture."
            )
        print(
            "Website synthetic example agrees with native policy results and exhaustive candidate checks"
        )
    else:
        output.write_text(text)


if __name__ == "__main__":
    main()
