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
from dotmatch.assayruntime import run_compiled_plan
from dotmatch.assaysim import simulate_panel
from dotmatch.core import load_targets
from dotmatch.assaywatch import WatchThresholds, watch_jsonl
from dotmatch.calibration_io import decode_tsv, fit_model_tsv, read_model, write_model


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


def command_execute(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="assaycode execute",
        description="Execute a compiled AssayScript v2 plan over synchronized FASTQ inputs.",
    )
    parser.add_argument("plan", help="compiled AssayScript plan JSON")
    parser.add_argument("--r1", help="R1 FASTQ or FASTQ.GZ")
    parser.add_argument("--r2", help="R2 FASTQ or FASTQ.GZ")
    parser.add_argument("--i1", help="I1 FASTQ or FASTQ.GZ")
    parser.add_argument("--i2", help="I2 FASTQ or FASTQ.GZ")
    parser.add_argument("--out", required=True, help="new or empty output directory")
    parser.add_argument("--max-reads", type=int)
    parser.add_argument(
        "--accept-findings",
        action="store_true",
        help="execute after explicitly reviewing findings recorded by the compiler",
    )
    args = parser.parse_args(list(argv))
    read_paths = {
        name: path
        for name, path in {"R1": args.r1, "R2": args.r2, "I1": args.i1, "I2": args.i2}.items()
        if path
    }
    try:
        result = run_compiled_plan(
            args.plan,
            read_paths,
            args.out,
            max_reads=args.max_reads,
            accept_findings=args.accept_findings,
        )
    except (AssayScriptError, OSError, ValueError) as exc:
        print(f"assaycode execute: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "status": "experimental",
                "output_dir": str(result.output_dir),
                "total_reads": result.total_reads,
                "status_counts": result.status_counts,
                "assignments": str(result.assignments),
                "counts": str(result.counts),
                "events": str(result.events),
                "summary": str(result.summary),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0




def command_calibrate(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="assaycode calibrate",
        description="Fit an experimental per-cycle error model from independently trusted TSV pairs.",
    )
    parser.add_argument("training", help="TSV with observed, expected, and quality columns")
    parser.add_argument("--out", required=True, help="model JSON")
    parser.add_argument("--prior-strength", type=float, default=100.0)
    args = parser.parse_args(list(argv))
    try:
        model = fit_model_tsv(args.training, prior_strength=args.prior_strength)
        output = write_model(model, args.out)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"assaycode calibrate: {exc}", file=sys.stderr)
        return 2
    print(json.dumps({
        "status": "experimental",
        "model": str(output),
        "cycles": len(model.cycle_totals),
        "observations": sum(model.cycle_totals),
        "errors": sum(model.cycle_errors),
    }, indent=2, sort_keys=True))
    return 0


def command_decode_quality(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="assaycode decode-quality",
        description="Apply an experimental calibrated selective decoder to short read windows.",
    )
    parser.add_argument("--reads", required=True, help="TSV with read_id, observed, and quality")
    parser.add_argument("--targets", required=True, help="TSV target table")
    parser.add_argument("--model", required=True, help="calibration model JSON")
    parser.add_argument("--out", required=True, help="calls TSV")
    parser.add_argument("--posterior-min", type=float, default=0.99)
    parser.add_argument("--likelihood-ratio-min", type=float, default=10.0)
    args = parser.parse_args(list(argv))
    try:
        model = read_model(args.model)
        summary = decode_tsv(
            args.reads,
            args.targets,
            model,
            args.out,
            posterior_min=args.posterior_min,
            likelihood_ratio_min=args.likelihood_ratio_min,
        )
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"assaycode decode-quality: {exc}", file=sys.stderr)
        return 2
    print(json.dumps({"status": "experimental", **summary}, indent=2, sort_keys=True))
    return 0

def command_simulate(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="assaycode simulate",
        description="Run an experimental deterministic substitution-error simulation for a target panel.",
    )
    parser.add_argument("--targets", required=True, help="TSV target table")
    parser.add_argument("--out", required=True, help="simulation result JSON")
    parser.add_argument("--reads-per-target", type=int, default=1000)
    parser.add_argument("--error-rate", type=float, default=0.01)
    parser.add_argument("-k", type=int, default=1)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args(list(argv))
    try:
        targets = dict(load_targets(args.targets))
        result = simulate_panel(
            targets,
            reads_per_target=args.reads_per_target,
            error_rate=args.error_rate,
            k=args.k,
            seed=args.seed,
        )
        output = Path(args.out)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    except (OSError, ValueError) as exc:
        print(f"assaycode simulate: {exc}", file=sys.stderr)
        return 2
    print(json.dumps({
        "status": "experimental",
        "result": str(output),
        "total_reads": result.total_reads,
        "usable_yield": result.usable_yield,
        "false_discovery_rate": result.false_discovery_rate,
    }, indent=2, sort_keys=True))
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
  assaycode execute assay.plan.json --r1 R1.fastq.gz --i1 I1.fastq.gz --out run
  assaycode calibrate trusted.tsv --out error-model.json
  assaycode decode-quality --reads windows.tsv --targets targets.tsv --model error-model.json --out calls.tsv
  assaycode simulate --targets targets.tsv --out simulation.json
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
  execute   stream synchronized FASTQs through a compiled multi-segment plan
  calibrate fit an experimental cycle-error model from trusted pairs
  decode-quality apply calibrated selective decoding with abstention
  simulate  estimate yield, ambiguity, no-calls, and FDR under an error model
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
  citation remain authoritative and supported.
  AssayCode is an additive assay-level identity, not a fork or destructive package rename.
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
    if raw_args[0] == "execute":
        return command_execute(raw_args[1:])
    if raw_args[0] == "calibrate":
        return command_calibrate(raw_args[1:])
    if raw_args[0] == "decode-quality":
        return command_decode_quality(raw_args[1:])
    if raw_args[0] == "simulate":
        return command_simulate(raw_args[1:])
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
