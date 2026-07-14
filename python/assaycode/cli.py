"""Compatibility-safe AssayCode command surface."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Sequence

from dotmatch import __version__
from dotmatch import cli as _engine_cli
from dotmatch.assayscript import AssayScriptError, load_and_compile, write_compiled_plan
from dotmatch.assaywatch import WatchThresholds, watch_jsonl


_SHORTCUTS = {"new", "infer", "check", "plan", "run", "start"}


def command_compile(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="assaycode compile",
        description="Validate AssayScript v2 and write a deterministic portable execution plan.",
    )
    parser.add_argument("spec", help="AssayScript v2 TOML")
    parser.add_argument("--out", required=True, help="compiled plan JSON")
    args = parser.parse_args(list(argv))
    try:
        plan = load_and_compile(args.spec)
        output = write_compiled_plan(plan, args.out)
    except AssayScriptError as exc:
        print(f"assaycode compile: {exc}", file=sys.stderr)
        return 2
    print(json.dumps({
        "status": "review" if plan.findings else "ready",
        "plan": str(output),
        "segments": len(plan.segments),
        "execution_order": plan.execution_order,
        "findings": plan.findings,
    }, indent=2, sort_keys=True))
    return 0


def command_inspect(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="assaycode inspect",
        description="Summarize a compiled AssayScript execution plan.",
    )
    parser.add_argument("plan", help="compiled plan JSON")
    args = parser.parse_args(list(argv))
    try:
        data = json.loads(Path(args.plan).read_text(encoding="utf-8"))
        segments = data["segments"]
        if data.get("compiler_schema_version") != 1 or not isinstance(segments, list):
            raise ValueError("unsupported compiled plan schema")
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"assaycode inspect: invalid compiled plan: {exc}", file=sys.stderr)
        return 2
    summary = {
        "name": data.get("name"),
        "assay_type": data.get("assay_type"),
        "segments": [
            {
                "name": segment.get("name"),
                "read": segment.get("read"),
                "strategy": segment.get("strategy"),
                "target_count": segment.get("target_count"),
                "safety_status": segment.get("safety_status"),
            }
            for segment in segments
        ],
        "execution_order": data.get("execution_order", []),
        "findings": data.get("findings", []),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0



def command_watch(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="assaycode watch",
        description="Stream assignment JSONL into bounded-memory sequential QC snapshots.",
    )
    parser.add_argument("events", help="assignment JSONL, or - for stdin")
    parser.add_argument("--out", default="-", help="snapshot JSONL, or - for stdout")
    parser.add_argument("--every", type=int, default=100000)
    parser.add_argument("--min-reads", type=int, default=1000)
    parser.add_argument("--min-assignment-rate", type=float, default=0.80)
    parser.add_argument("--max-ambiguous-rate", type=float, default=0.05)
    parser.add_argument("--max-unmatched-rate", type=float, default=0.15)
    parser.add_argument("--max-invalid-rate", type=float, default=0.02)
    parser.add_argument("--fail-on-review", action="store_true")
    args = parser.parse_args(list(argv))
    try:
        latest = watch_jsonl(
            args.events,
            args.out,
            every=args.every,
            thresholds=WatchThresholds(
                min_assignment_rate=args.min_assignment_rate,
                max_ambiguous_rate=args.max_ambiguous_rate,
                max_unmatched_rate=args.max_unmatched_rate,
                max_invalid_rate=args.max_invalid_rate,
                min_reads=args.min_reads,
            ),
        )
    except (OSError, ValueError) as exc:
        print(f"assaycode watch: {exc}", file=sys.stderr)
        return 2
    if args.fail_on_review and latest is not None and latest.decision == "review":
        return 1
    return 0

def print_help() -> None:
    print(
        f"""AssayCode {__version__} — powered by the DotMatch engine

Compile, validate, decode, and diagnose sequencing assays built from known
guides, barcodes, primers, feature tags, and panel targets.

Usage:
  assaycode --help
  assaycode --version
  assaycode compile assay-v2.toml --out assay.plan.json
  assaycode inspect assay.plan.json
  assaycode watch assignments.jsonl --out snapshots.jsonl
  assaycode check assay.toml
  assaycode plan assay.toml
  assaycode run assay.toml
  assaycode start assay.toml
  assaycode assay <command> [options]
  assaycode engine <dotmatch-command> [options]

AssayScript v2:
  compile   validate a multi-read specification and select deterministic strategies
  inspect   summarize a compiled plan, safety status, fingerprints, and findings
  watch     stream assignment events into sequential QC decisions

AssaySpec v1 workflow shortcuts:
  new       scaffold an AssayScript/AssaySpec project
  infer     infer a reviewable assay window from reads
  check     validate a specification and its referenced inputs
  plan      print the deterministic execution plan
  run       execute an already-reviewed specification
  start     check, plan, run, and write reliability artifacts

Specialized DotMatch workflows remain available unchanged:
  assaycode crispr ...
  assaycode barcode ...
  assaycode panel ...
  assaycode count ...
  assaycode demux ...

Compatibility:
  The dotmatch executable, Python package, native ABI, output schemas, DOI, and
  citation remain authoritative and supported. AssayCode is an additive
  assay-level identity, not a fork or destructive package rename.
"""
    )


def main(argv: Sequence[str] | None = None) -> int:
    raw_args = list(sys.argv[1:] if argv is None else argv)
    if not raw_args or raw_args[0] in {"-h", "--help", "help"}:
        print_help()
        return 0
    if raw_args == ["--version"]:
        print(f"assaycode {__version__} (DotMatch engine {__version__})")
        return 0
    if raw_args[0] == "compile":
        return command_compile(raw_args[1:])
    if raw_args[0] == "inspect":
        return command_inspect(raw_args[1:])
    if raw_args[0] == "watch":
        return command_watch(raw_args[1:])
    if raw_args[0] == "engine":
        if len(raw_args) == 1:
            print("usage: assaycode engine <dotmatch-command> [options]", file=sys.stderr)
            return 2
        return _engine_cli.main(raw_args[1:])
    if raw_args[0] in _SHORTCUTS:
        return _engine_cli.main(["assay", *raw_args])
    return _engine_cli.main(raw_args)


if __name__ == "__main__":
    raise SystemExit(main())
