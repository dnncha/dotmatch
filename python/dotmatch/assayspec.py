from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
from importlib import resources as importlib_resources
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence
from urllib.parse import quote

from . import __version__ as PYTHON_PACKAGE_VERSION
from .core import MATCH_AMBIGUOUS, MATCH_NONE, MATCH_UNIQUE, Matcher
from .native import find_native_cli

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - exercised on Python 3.9/3.10 through dependency
    import tomli as tomllib  # type: ignore[no-redef]


MODES = {"count", "demux", "pair-count"}
ASSAY_TYPES = {"crispr", "feature_barcode", "inline_barcode", "amplicon_panel", "oligo_adapter", "generic"}
METRICS = {"hamming", "levenshtein"}
AMBIGUITY_POLICIES = {"best", "radius"}
AMBIGUOUS_OUTPUT = {"discard", "report"}
SPEC_STATUS = {"ready", "draft"}
AUTOPSY_THRESHOLDS = {
    "assignment_rate_min": 0.80,
    "ambiguous_rate_max": 0.05,
    "no_match_rate_max": 0.15,
    "invalid_rate_max": 0.02,
}
CRISPR_REPRESENTATION_THRESHOLDS = {
    "min_coverage_fraction": 0.90,
    "max_zero_count_fraction": 0.10,
    "max_gini_index": 0.50,
    "max_top_1pct_fraction": 0.30,
    "min_pairwise_sample_pearson": 0.80,
}
RELIABILITY_PROFILES = {"production", "exploratory"}
BACKEND_MODES = {"auto", "cpu", "gpu-metal-experimental"}
RELIABILITY_DEFAULTS: dict[str, Any] = {
    "profile": "production",
    "fail_on_unsafe_targets": True,
    "fail_on_draft_inference": True,
    "min_assignment_rate": AUTOPSY_THRESHOLDS["assignment_rate_min"],
    "max_ambiguous_rate": AUTOPSY_THRESHOLDS["ambiguous_rate_max"],
    "max_unmatched_rate": AUTOPSY_THRESHOLDS["no_match_rate_max"],
    "max_invalid_rate": AUTOPSY_THRESHOLDS["invalid_rate_max"],
    "min_coverage_fraction": CRISPR_REPRESENTATION_THRESHOLDS["min_coverage_fraction"],
    "max_zero_count_fraction": CRISPR_REPRESENTATION_THRESHOLDS["max_zero_count_fraction"],
    "max_gini_index": CRISPR_REPRESENTATION_THRESHOLDS["max_gini_index"],
    "max_top_1pct_fraction": CRISPR_REPRESENTATION_THRESHOLDS["max_top_1pct_fraction"],
    "min_pairwise_sample_pearson": CRISPR_REPRESENTATION_THRESHOLDS["min_pairwise_sample_pearson"],
    "require_public_evidence_boundary": True,
}
BACKEND_DEFAULTS: dict[str, Any] = {
    "mode": "auto",
    "allow_gpu": True,
}
GPU_BENCHMARK_PRIORS = [
    {
        "workload": "public_crispr_yusa_hamming",
        "mode": "count",
        "assay_type": "crispr",
        "target_count": 87437,
        "target_length": 19,
        "total_speedup": 1.92,
    },
    {
        "workload": "synthetic_hamming_737_targets",
        "mode": "count",
        "assay_type": "generic",
        "target_count": 737,
        "target_length": 20,
        "total_speedup": 13.50,
    },
    {
        "workload": "synthetic_hamming_4096_targets",
        "mode": "count",
        "assay_type": "generic",
        "target_count": 4096,
        "target_length": 20,
        "total_speedup": 10.05,
    },
]
RELIABILITY_FINDING_COLUMNS = [
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
]
ASSAY_FIX_COLUMNS = [
    "fix_id",
    "finding_id",
    "section",
    "key",
    "current_value",
    "suggested_value",
    "rationale",
]
ASSAY_EVIDENCE_IDS = {
    "crispr": "crispr_guide_counting",
    "feature_barcode": "feature_barcode",
    "inline_barcode": "inline_barcode",
    "amplicon_panel": "amplicon_panel",
    "oligo_adapter": "oligo_adapter",
    "generic": "paired_combinatorial",
}
TEMPLATES = {
    "crispr",
    "feature-barcode",
    "inline-barcode-count",
    "inline-barcode-demux",
    "amplicon-panel",
    "oligo-adapter",
    "pair-count",
}
SCAFFOLD_TEMPLATES = {
    "crispr": {
        "mode": "count",
        "assay_type": "crispr",
        "input_key": "targets",
        "format": "mageck",
    },
    "feature-barcode": {
        "mode": "count",
        "assay_type": "feature_barcode",
        "input_key": "targets",
        "format": "dotmatch",
    },
    "inline-barcode-count": {
        "mode": "count",
        "assay_type": "inline_barcode",
        "input_key": "targets",
        "format": "dotmatch",
    },
    "amplicon-panel": {
        "mode": "count",
        "assay_type": "amplicon_panel",
        "input_key": "targets",
        "format": "dotmatch",
    },
    "oligo-adapter": {
        "mode": "count",
        "assay_type": "oligo_adapter",
        "input_key": "targets",
        "format": "dotmatch",
    },
    "inline-barcode-demux": {
        "mode": "demux",
        "assay_type": "inline_barcode",
        "input_key": "barcodes",
        "format": "",
    },
}
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


class AssaySpecError(ValueError):
    pass


@dataclass(frozen=True)
class AssaySpec:
    path: Path
    data: Mapping[str, Any]

    @property
    def mode(self) -> str:
        return str(self.data["mode"])

    @property
    def assay_type(self) -> str:
        return str(self.data["assay_type"])

    @property
    def out_dir(self) -> Path:
        return _path_from_spec(
            self.path,
            str(_table(self.data, "run").get("out_dir", "dotmatch_assay_out")),
            allow_absolute=True,
            name="run.out_dir",
        )

    @property
    def k(self) -> int:
        return int(_table(self.data, "assignment").get("k", 1))

    @property
    def status(self) -> str:
        return str(self.data.get("status", "ready"))

    @property
    def reliability(self) -> dict[str, Any]:
        config = dict(RELIABILITY_DEFAULTS)
        config.update(_table(self.data, "reliability"))
        return config

    @property
    def backend(self) -> dict[str, Any]:
        config = dict(BACKEND_DEFAULTS)
        config.update(_table(self.data, "backend"))
        return config


@dataclass(frozen=True)
class PlanStep:
    name: str
    argv: list[str]
    warning_ok: bool = False


@dataclass(frozen=True)
class AssayPlan:
    spec: AssaySpec
    steps: list[PlanStep]
    artifacts: dict[str, Path]
    generated_files: dict[str, Path]


def load_assay_spec(path: str | Path) -> AssaySpec:
    spec_path = Path(path)
    try:
        with spec_path.open("rb") as fh:
            data = tomllib.load(fh)
    except tomllib.TOMLDecodeError as exc:
        raise AssaySpecError(f"{spec_path}: invalid TOML: {exc}") from exc
    if not isinstance(data, dict):
        raise AssaySpecError(f"{spec_path}: top-level TOML document must be a table")
    assay = AssaySpec(path=spec_path, data=data)
    validate_assay_spec(assay)
    return assay


def validate_assay_spec(assay: AssaySpec) -> None:
    data = assay.data
    _require_equal(data.get("schema_version"), 1, "schema_version")
    _require_enum(data.get("status", "ready"), SPEC_STATUS, "status")
    _require_enum(data.get("mode"), MODES, "mode")
    _require_enum(data.get("assay_type"), ASSAY_TYPES, "assay_type")

    assignment = _table(data, "assignment")
    if "k" in assignment:
        _require_int_range(assignment["k"], 0, 2, "assignment.k")
    _require_enum(assignment.get("metric", "levenshtein"), METRICS, "assignment.metric")
    _require_enum(assignment.get("ambiguity_policy", "radius"), AMBIGUITY_POLICIES, "assignment.ambiguity_policy")
    _require_enum(assignment.get("ambiguous", "discard"), AMBIGUOUS_OUTPUT, "assignment.ambiguous")
    if int(assignment.get("k", 1)) == 2 and assignment.get("metric", "levenshtein") == "hamming":
        raise AssaySpecError("assignment.k=2 is only valid with assignment.metric='levenshtein'")

    reliability = _table(data, "reliability")
    if "profile" in reliability:
        _require_enum(reliability["profile"], RELIABILITY_PROFILES, "reliability.profile")
    for key in ["fail_on_unsafe_targets", "fail_on_draft_inference", "require_public_evidence_boundary"]:
        if key in reliability:
            _require_bool(reliability[key], f"reliability.{key}")
    for key in [
        "min_assignment_rate",
        "max_ambiguous_rate",
        "max_unmatched_rate",
        "max_invalid_rate",
        "min_coverage_fraction",
        "max_zero_count_fraction",
        "max_gini_index",
        "max_top_1pct_fraction",
    ]:
        if key in reliability:
            _require_float_range(reliability[key], 0.0, 1.0, f"reliability.{key}")
    if "min_pairwise_sample_pearson" in reliability:
        _require_float_range(reliability["min_pairwise_sample_pearson"], -1.0, 1.0, "reliability.min_pairwise_sample_pearson")

    backend = _table(data, "backend")
    if "mode" in backend:
        _require_enum(backend["mode"], BACKEND_MODES, "backend.mode")
    if "allow_gpu" in backend:
        _require_bool(backend["allow_gpu"], "backend.allow_gpu")

    mode = str(data["mode"])
    if mode == "count":
        _require_path(assay, "targets")
        samples = data.get("samples")
        if not isinstance(samples, list) or not samples:
            raise AssaySpecError("samples must contain at least one [[samples]] entry")
        for i, sample in enumerate(samples):
            if not isinstance(sample, dict):
                raise AssaySpecError(f"samples[{i}] must be a table")
            if not sample.get("id"):
                raise AssaySpecError(f"samples[{i}].id is required")
            _require_safe_identifier(sample.get("id"), f"samples[{i}].id")
            _require_existing_path(assay, sample.get("fastq"), f"samples[{i}].fastq")
        _require_extract(data, "extract")
    elif mode == "demux":
        _require_path(assay, "barcodes")
        _require_path(assay, "reads")
        _require_extract(data, "extract", allow_auto=True)
    else:
        _require_path(assay, "left_targets")
        _require_path(assay, "right_targets")
        _require_path(assay, "reads")
        _require_extract(data, "left")
        _require_extract(data, "right")


def compile_assay_plan(assay: AssaySpec) -> AssayPlan:
    out_dir = assay.out_dir
    generated: dict[str, Path] = {}
    artifacts: dict[str, Path] = {
        "manifest": out_dir / "assay_manifest.json",
        "manifest_summary": out_dir / "assay_manifest.summary.tsv",
        "assay_report": out_dir / "assay_report.html",
        "normalized_spec": out_dir / "assay.normalized.json",
        "reliability_summary": out_dir / "reliability_summary.json",
        "reliability_findings": out_dir / "reliability_findings.tsv",
        "reliability_report": out_dir / "reliability_report.html",
        "reliability_manifest_summary": out_dir / "reliability_manifest.summary.tsv",
        "assay_fixes": out_dir / "assay_fixes.tsv",
        "backend_optimization": out_dir / "backend_optimization.json",
        "methods": out_dir / "methods.md",
        "citation_bib": out_dir / "CITATION.bib",
        "software_versions": out_dir / "software_versions.yml",
    }
    steps: list[PlanStep] = []

    if assay.mode == "count":
        audit_dir = out_dir / "audit"
        artifacts["audit"] = audit_dir
        steps.append(PlanStep("audit", _audit_cmd(_spec_path(assay, "targets"), audit_dir, assay.k), warning_ok=True))
        samples_path = out_dir / "assay_samples.tsv"
        generated["samples"] = samples_path
        _compile_count(assay, steps, artifacts, samples_path)
    elif assay.mode == "demux":
        audit_dir = out_dir / "audit"
        artifacts["audit"] = audit_dir
        steps.append(PlanStep("audit", _audit_cmd(_spec_path(assay, "barcodes"), audit_dir, assay.k), warning_ok=True))
        _compile_demux(assay, steps, artifacts)
    else:
        left_audit = out_dir / "audit" / "left"
        right_audit = out_dir / "audit" / "right"
        artifacts["left_audit"] = left_audit
        artifacts["right_audit"] = right_audit
        steps.append(PlanStep("audit-left", _audit_cmd(_spec_path(assay, "left_targets"), left_audit, assay.k), warning_ok=True))
        steps.append(PlanStep("audit-right", _audit_cmd(_spec_path(assay, "right_targets"), right_audit, assay.k), warning_ok=True))
        _compile_pair(assay, steps, artifacts)

    return AssayPlan(spec=assay, steps=steps, artifacts=artifacts, generated_files=generated)


def format_plan(plan: AssayPlan) -> str:
    lines = [shlex.join(step.argv) for step in plan.steps]
    lines.extend(
        [
            "# Reliability artifacts",
            f"# reliability_summary: {plan.artifacts['reliability_summary']}",
            f"# reliability_findings: {plan.artifacts['reliability_findings']}",
            f"# reliability_report: {plan.artifacts['reliability_report']}",
            f"# reliability_manifest_summary: {plan.artifacts['reliability_manifest_summary']}",
            f"# assay_fixes: {plan.artifacts['assay_fixes']}",
            f"# methods: {plan.artifacts['methods']}",
            f"# citation_bib: {plan.artifacts['citation_bib']}",
            f"# software_versions: {plan.artifacts['software_versions']}",
        ]
    )
    return "\n".join(lines) + "\n"


def run_assay_plan(plan: AssayPlan, *, skip_audit: bool = False) -> int:
    out_dir = plan.spec.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_generated_files(plan)
    _write_normalized_spec(plan)
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "mode": plan.spec.mode,
        "assay_type": plan.spec.assay_type,
        "status": plan.spec.status,
        "spec_path": str(plan.spec.path),
        "native_cli": "",
        "commands": [],
        "artifacts": {key: str(value) for key, value in plan.artifacts.items()},
        "inference_report": str(plan.spec.data.get("inference_report", "")),
        "autopsy_triggered": False,
        "autopsy_thresholds": _reliability_thresholds(plan.spec),
        "autopsy_artifacts": {},
        "production_warnings": [],
        "warnings": [],
    }
    try:
        native = find_native_cli()
    except FileNotFoundError as exc:
        manifest["native_cli"] = ""
        manifest["production_warnings"].append(str(exc))
        reliability = _build_reliability_summary(plan, stage="preflight", manifest=manifest)
        reliability["findings"].append(
            _finding(
                "native_cli_missing",
                "blocked",
                "preflight",
                "",
                "native_cli",
                "missing",
                "executable",
                str(exc),
                "Build with make dotmatch, install a wheel with the bundled native executable, or set DOTMATCH_NATIVE_CLI.",
                "",
            )
        )
        reliability["finding_counts"] = _finding_counts(reliability["findings"])
        reliability["overall_status"] = _overall_reliability_status(reliability["finding_counts"])
        manifest["reliability"] = {
            "overall_status": reliability["overall_status"],
            "finding_counts": reliability["finding_counts"],
            "summary": str(plan.artifacts["reliability_summary"]),
            "report": str(plan.artifacts["reliability_report"]),
        }
        _write_reliability_artifacts(plan, reliability)
        _write_manifest(plan, manifest)
        print(f"dotmatch assay: {exc}", file=sys.stderr)
        return 2
    manifest["native_cli"] = str(native)
    version = subprocess.run([str(native), "--version"], check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    manifest["native_version"] = version.stdout.strip() if version.returncode == 0 else ""

    for step in plan.steps:
        if skip_audit and step.name.startswith("audit"):
            continue
        if step.name.startswith("audit"):
            Path(step.argv[-1]).parent.mkdir(parents=True, exist_ok=True)
        argv = _resolve_native(step.argv, native)
        result = subprocess.run(argv, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        command_record = {
            "name": step.name,
            "argv": argv,
            "exit_code": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }
        manifest["commands"].append(command_record)
        if step.name.startswith("audit") and result.returncode == 0:
            _append_audit_warnings(plan, step, manifest, warn=not _audit_step_unsafe(plan, step))
            if _audit_step_unsafe(plan, step) and _blocks_on_unsafe_targets(plan.spec):
                message = f"{step.name}: unsafe target audit at k={plan.spec.k}; blocked by production reliability profile"
                manifest["production_warnings"].append(message)
                reliability = _build_reliability_summary(plan, stage="preflight", manifest=manifest)
                manifest["reliability"] = {
                    "overall_status": reliability["overall_status"],
                    "finding_counts": reliability["finding_counts"],
                    "summary": str(plan.artifacts["reliability_summary"]),
                    "report": str(plan.artifacts["reliability_report"]),
                }
                _write_reliability_artifacts(plan, reliability)
                _write_manifest(plan, manifest)
                print(f"dotmatch assay: {message}", file=sys.stderr)
                return 2
        if result.returncode != 0:
            reliability = _build_reliability_summary(plan, stage="postrun", manifest=manifest)
            manifest["reliability"] = {
                "overall_status": reliability["overall_status"],
                "finding_counts": reliability["finding_counts"],
                "summary": str(plan.artifacts["reliability_summary"]),
                "report": str(plan.artifacts["reliability_report"]),
            }
            _write_reliability_artifacts(plan, reliability)
            _write_manifest(plan, manifest)
            sys.stderr.write(result.stderr)
            return int(result.returncode)

    autopsy_reasons = _autopsy_trigger_reasons(plan)
    if autopsy_reasons:
        autopsy_dir = plan.spec.out_dir / "autopsy"
        autopsy_result = run_autopsy(plan.spec, autopsy_dir)
        manifest["autopsy_triggered"] = True
        manifest["autopsy_artifacts"] = {key: str(value) for key, value in autopsy_result.items()}
        manifest["production_warnings"].extend(autopsy_reasons)

    reliability = _build_reliability_summary(plan, stage="postrun", manifest=manifest)
    manifest["reliability"] = {
        "overall_status": reliability["overall_status"],
        "finding_counts": reliability["finding_counts"],
        "summary": str(plan.artifacts["reliability_summary"]),
        "report": str(plan.artifacts["reliability_report"]),
    }
    _write_reliability_artifacts(plan, reliability)
    _write_manifest(plan, manifest)
    return 0


def command_assay(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(prog="dotmatch assay", description="Validate, plan, and run DotMatch AssaySpec TOML workflows.")
    sub = parser.add_subparsers(dest="command", required=True)

    init = sub.add_parser("init", help="write a starter AssaySpec TOML file")
    init.add_argument("--template", required=True, choices=sorted(TEMPLATES))
    init.add_argument("--out", required=True)

    infer = sub.add_parser("infer", help="infer a fixed-window AssaySpec from FASTQ and target tables")
    infer.add_argument("--mode", required=True, choices=sorted(MODES))
    infer.add_argument("--assay-type", required=True, choices=sorted(ASSAY_TYPES))
    infer.add_argument("--targets")
    infer.add_argument("--barcodes")
    infer.add_argument("--left-targets")
    infer.add_argument("--right-targets")
    infer.add_argument("--reads", required=True)
    infer.add_argument("--sample-id", default="sample")
    infer.add_argument("--out", required=True)
    infer.add_argument("--report", required=True)
    infer.add_argument("--candidates")
    infer.add_argument("--max-reads", type=int, default=50000)
    infer.add_argument("--max-start", type=int, default=32)

    new = sub.add_parser(
        "new",
        help="scaffold a reviewable assay project from a target table and a directory of FASTQs",
    )
    new.add_argument("template", choices=sorted(SCAFFOLD_TEMPLATES))
    new.add_argument("--out", required=True, help="empty project directory to create")
    new.add_argument("--targets", help="target or guide table for count scaffolds")
    new.add_argument("--library", help="alias for --targets on CRISPR scaffolds")
    new.add_argument("--barcodes", help="barcode table for demux scaffolds")
    new.add_argument("--reads-dir", required=True, help="directory containing sample FASTQ files")
    new.add_argument("--link-reads", action="store_true", help="symlink reads instead of copying them into the project")
    new.add_argument("--threads", type=int, default=1)
    new.add_argument("--max-reads", type=int, default=50000)
    new.add_argument("--max-start", type=int, default=32)

    autopsy = sub.add_parser("autopsy", help="diagnose suspicious fixed-window assay runs")
    autopsy.add_argument("spec")
    autopsy.add_argument("--out-dir", required=True)

    workflow_help = {
        "check": "validate an AssaySpec and write preflight reliability artifacts",
        "optimize": "write a benchmark-informed CPU/GPU backend recommendation",
        "plan": "print the native commands and outputs for an AssaySpec",
        "run": "run an AssaySpec workflow and write reports, manifests, and outputs",
        "start": "check an AssaySpec, run it, and print the reliability verdict",
    }
    for name, help_text in workflow_help.items():
        child = sub.add_parser(name, help=help_text)
        child.add_argument("spec")
        if name == "start":
            child.add_argument(
                "--check-only",
                action="store_true",
                help="run preflight check only (same as assay check); skip counting",
            )

    args = parser.parse_args(list(argv))
    try:
        if args.command == "init":
            return _command_init(args.template, Path(args.out))
        if args.command == "new":
            return _command_new(args)
        if args.command == "infer":
            return _command_infer(args)
        assay = load_assay_spec(args.spec)
        if args.command == "autopsy":
            run_autopsy(assay, Path(args.out_dir))
            return 0
        if args.command == "check":
            plan = compile_assay_plan(assay)
            return _command_assay_check(plan)
        if args.command == "start" and getattr(args, "check_only", False):
            plan = compile_assay_plan(assay)
            return _command_assay_check(plan)
        if args.command == "optimize":
            plan = compile_assay_plan(assay)
            optimization = optimize_assay_backend(assay)
            _write_backend_optimization(plan, optimization)
            print(_format_backend_optimization(optimization), end="")
            return 0
        if args.command in {"run", "start"} and assay.status == "draft" and _blocks_on_draft_inference(assay):
            plan = compile_assay_plan(assay)
            reliability = _build_reliability_summary(plan, stage="preflight")
            _write_reliability_artifacts(plan, reliability)
            if args.command == "start":
                print(f"{_spec_user_label(assay.path)}: preflight blocked (draft_assayspec)", file=sys.stderr)
                _print_reliability_verdict(plan)
                return 2
            raise AssaySpecError("refusing to run draft AssaySpec; review inference report and promote status to 'ready'")
        plan = compile_assay_plan(assay)
        if args.command == "plan":
            print(format_plan(plan), end="")
            return 0
        if args.command == "start":
            preflight_status = _preflight_for_assay_start(plan)
            spec_label = _spec_user_label(assay.path)
            if preflight_status == "blocked":
                return _finish_assay_preflight(plan, spec_label, preflight_status, verb="preflight")
            if preflight_status == "passed":
                print(f"{spec_label}: preflight passed", file=sys.stderr)
            else:
                summary = _read_reliability_summary(plan) or {}
                reason = _primary_reliability_reason(summary)
                print(
                    f"{spec_label}: preflight {preflight_status.replace('_', ' ')} ({reason}); continuing run",
                    file=sys.stderr,
                )
                _print_reliability_verdict(plan, label="preflight", include_next=False)
            print(f"{spec_label}: running", file=sys.stderr)
            rc = run_assay_plan(plan, skip_audit=True)
            _print_reliability_verdict(plan, label="postrun")
            if rc != 0:
                return rc
            return _reliability_exit_code(_read_reliability_status(plan))
        return run_assay_plan(plan)
    except AssaySpecError as exc:
        print(f"dotmatch assay: {exc}", file=sys.stderr)
        return 2
    except FileNotFoundError as exc:
        print(f"dotmatch assay: {exc}", file=sys.stderr)
        return 2


def _command_new(args: argparse.Namespace) -> int:
    targets = Path(args.targets) if args.targets else None
    if targets is None and args.library:
        targets = Path(args.library)
    try:
        result = scaffold_assay_project(
            template=args.template,
            project_dir=Path(args.out),
            reads_dir=Path(args.reads_dir),
            targets=targets,
            barcodes=Path(args.barcodes) if args.barcodes else None,
            link_reads=bool(args.link_reads),
            threads=args.threads,
            max_reads=args.max_reads,
            max_start=args.max_start,
        )
    except AssaySpecError as exc:
        print(f"dotmatch assay: {exc}", file=sys.stderr)
        return 2
    print(str(result["project"]))
    return 0


def _command_infer(args: argparse.Namespace) -> int:
    try:
        result = infer_assay_spec(
            mode=args.mode,
            assay_type=args.assay_type,
            reads=Path(args.reads),
            out=Path(args.out),
            report=Path(args.report),
            candidates_path=Path(args.candidates) if args.candidates else None,
            targets=Path(args.targets) if args.targets else None,
            barcodes=Path(args.barcodes) if args.barcodes else None,
            left_targets=Path(args.left_targets) if args.left_targets else None,
            right_targets=Path(args.right_targets) if args.right_targets else None,
            sample_id=args.sample_id,
            max_reads=args.max_reads,
            max_start=args.max_start,
        )
    except AssaySpecError as exc:
        print(f"dotmatch assay: {exc}", file=sys.stderr)
        return 2
    print(str(result["spec"]))
    return 0


def infer_assay_spec(
    *,
    mode: str,
    assay_type: str,
    reads: Path,
    out: Path,
    report: Path,
    candidates_path: Path | None = None,
    targets: Path | None = None,
    barcodes: Path | None = None,
    left_targets: Path | None = None,
    right_targets: Path | None = None,
    sample_id: str = "sample",
    max_reads: int = 50000,
    max_start: int = 32,
) -> dict[str, Path]:
    if mode not in MODES:
        raise AssaySpecError(f"mode must be one of: {', '.join(sorted(MODES))}")
    if assay_type not in ASSAY_TYPES:
        raise AssaySpecError(f"assay_type must be one of: {', '.join(sorted(ASSAY_TYPES))}")
    if not reads.exists():
        raise AssaySpecError(f"reads does not exist: {reads}")
    read_seqs = _read_fastq_sequences(reads, max_reads=max_reads)
    if not read_seqs:
        raise AssaySpecError(f"reads contains no FASTQ records: {reads}")

    if mode == "count":
        if targets is None:
            raise AssaySpecError("--targets is required for count inference")
        target_set = _read_target_sequences(targets)
        candidates = _score_windows(read_seqs, target_set.sequences, max_start=max_start)
        chosen, status, warnings = _choose_candidate(candidates)
        _write_inferred_count_spec(out, status, assay_type, targets, reads, sample_id, chosen)
        report_data: dict[str, Any] = {
            "mode": mode,
            "assay_type": assay_type,
            "status": status,
            "chosen": chosen,
            "warnings": warnings,
            "candidates": candidates,
        }
    elif mode == "demux":
        if barcodes is None:
            raise AssaySpecError("--barcodes is required for demux inference")
        target_set = _read_target_sequences(barcodes)
        candidates = _score_windows(read_seqs, target_set.sequences, max_start=max_start)
        chosen, status, warnings = _choose_candidate(candidates)
        _write_inferred_demux_spec(out, status, assay_type, barcodes, reads, chosen)
        report_data = {
            "mode": mode,
            "assay_type": assay_type,
            "status": status,
            "chosen": chosen,
            "warnings": warnings,
            "candidates": candidates,
        }
    else:
        if left_targets is None or right_targets is None:
            raise AssaySpecError("--left-targets and --right-targets are required for pair-count inference")
        left_set = _read_target_sequences(left_targets)
        right_set = _read_target_sequences(right_targets)
        left_candidates = _score_windows(read_seqs, left_set.sequences, max_start=max_start)
        right_candidates = _score_windows(read_seqs, right_set.sequences, max_start=max_start)
        left_chosen, left_status, left_warnings = _choose_candidate(left_candidates)
        right_chosen, right_status, right_warnings = _choose_candidate(right_candidates)
        status = "ready" if left_status == "ready" and right_status == "ready" else "draft"
        _write_inferred_pair_spec(out, status, assay_type, left_targets, right_targets, reads, left_chosen, right_chosen)
        report_data = {
            "mode": mode,
            "assay_type": assay_type,
            "status": status,
            "left": {"chosen": left_chosen, "warnings": left_warnings, "candidates": left_candidates},
            "right": {"chosen": right_chosen, "warnings": right_warnings, "candidates": right_candidates},
            "warnings": left_warnings + right_warnings,
        }

    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps(report_data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    candidate_out = candidates_path or report.with_name("inference_candidates.tsv")
    _write_candidates_tsv(candidate_out, report_data)
    return {"spec": out, "report": report, "candidates": candidate_out}


def scaffold_assay_project(
    *,
    template: str,
    project_dir: Path,
    reads_dir: Path,
    targets: Path | None = None,
    barcodes: Path | None = None,
    link_reads: bool = False,
    threads: int = 1,
    max_reads: int = 50000,
    max_start: int = 32,
) -> dict[str, Path]:
    if template not in SCAFFOLD_TEMPLATES:
        raise AssaySpecError(
            f"template must be one of: {', '.join(sorted(SCAFFOLD_TEMPLATES))}; use dotmatch assay init for pair-count"
        )
    if project_dir.exists() and any(project_dir.iterdir()):
        raise AssaySpecError(f"refusing to scaffold into non-empty directory: {project_dir}")
    if not reads_dir.is_dir():
        raise AssaySpecError(f"reads-dir does not exist or is not a directory: {reads_dir}")

    config = SCAFFOLD_TEMPLATES[template]
    mode = config["mode"]
    assay_type = config["assay_type"]
    input_key = config["input_key"]
    table_path = barcodes if input_key == "barcodes" else targets
    if table_path is None:
        flag = "--barcodes" if input_key == "barcodes" else "--targets/--library"
        raise AssaySpecError(f"{flag} is required for template {template}")
    if not table_path.exists():
        raise AssaySpecError(f"{input_key} does not exist: {table_path}")

    fastqs = _discover_fastqs(reads_dir)
    if not fastqs:
        raise AssaySpecError(f"reads-dir contains no FASTQ files: {reads_dir}")
    if mode == "demux" and len(fastqs) != 1:
        raise AssaySpecError("demux scaffold requires exactly one FASTQ in reads-dir")

    project_dir.mkdir(parents=True, exist_ok=True)
    inputs_dir = project_dir / "inputs"
    reads_root = project_dir / "reads"
    inputs_dir.mkdir()
    reads_root.mkdir()

    if input_key == "barcodes":
        staged_table = inputs_dir / barcodes.name
        _stage_input_file(barcodes, staged_table, copy=True)
        staged_table_key = "barcodes"
    else:
        staged_table = inputs_dir / targets.name
        _stage_input_file(targets, staged_table, copy=True)
        staged_table_key = "targets"

    staged_samples: list[tuple[str, Path, Path]] = []
    used_ids: set[str] = set()
    for source_fastq in fastqs:
        sample_id = _unique_sample_id(_sample_id_from_fastq(source_fastq), used_ids)
        staged_fastq = reads_root / source_fastq.name
        _stage_input_file(source_fastq, staged_fastq, copy=not link_reads)
        staged_samples.append((sample_id, staged_fastq, source_fastq))

    infer_reads = staged_samples[0][1]
    read_seqs = _read_pooled_fastq_sequences([sample for _, sample, _ in staged_samples], max_reads=max_reads)
    if not read_seqs:
        raise AssaySpecError(f"reads contains no FASTQ records: {infer_reads}")

    spec_path = project_dir / "assay.toml"
    report_path = project_dir / "inference_report.json"
    candidates_path = project_dir / "inference_candidates.tsv"
    samples_tsv = project_dir / "samples.generated.tsv"

    if mode == "count":
        target_set = _read_target_sequences(staged_table)
        candidates = _score_windows(read_seqs, target_set.sequences, max_start=max_start)
        chosen, status, warnings = _choose_candidate(candidates)
        _write_scaffolded_count_spec(
            spec_path,
            project_dir=project_dir,
            status=status,
            assay_type=assay_type,
            targets=staged_table,
            samples=staged_samples,
            chosen=chosen,
            output_format=config["format"],
            threads=threads,
        )
        report_data: dict[str, Any] = {
            "schema_version": 1,
            "command": "assay new",
            "template": template,
            "mode": mode,
            "assay_type": assay_type,
            "status": status,
            "chosen": chosen,
            "warnings": warnings,
            "candidates": candidates,
            "samples": [
                {"sample_id": sample_id, "fastq": str(staged_fastq.relative_to(project_dir)), "source_fastq": str(source)}
                for sample_id, staged_fastq, source in staged_samples
            ],
            "inference_reads": str(infer_reads.relative_to(project_dir)),
            "inference_read_sources": [
                str(staged_fastq.relative_to(project_dir)) for _, staged_fastq, _ in staged_samples
            ],
        }
    else:
        target_set = _read_target_sequences(staged_table)
        candidates = _score_windows(read_seqs, target_set.sequences, max_start=max_start)
        chosen, status, warnings = _choose_candidate(candidates)
        sample_id, staged_fastq, source_fastq = staged_samples[0]
        _write_scaffolded_demux_spec(
            spec_path,
            project_dir=project_dir,
            status=status,
            assay_type=assay_type,
            barcodes=staged_table,
            reads=staged_fastq,
            chosen=chosen,
        )
        report_data = {
            "schema_version": 1,
            "command": "assay new",
            "template": template,
            "mode": mode,
            "assay_type": assay_type,
            "status": status,
            "chosen": chosen,
            "warnings": warnings,
            "candidates": candidates,
            "samples": [
                {
                    "sample_id": sample_id,
                    "fastq": str(staged_fastq.relative_to(project_dir)),
                    "source_fastq": str(source_fastq),
                }
            ],
            "inference_reads": str(infer_reads.relative_to(project_dir)),
        }

    report_path.write_text(json.dumps(report_data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_candidates_tsv(candidates_path, report_data)
    _write_generated_samples_tsv(samples_tsv, staged_samples, project_dir)
    _write_scaffold_readme(project_dir, template=template, status=status, warnings=warnings)
    _write_scaffold_run_script(project_dir, status=status, launcher=_detect_dotmatch_launcher(), native_cli=_detect_native_cli_for_scaffold())
    return {
        "project": project_dir,
        "spec": spec_path,
        "report": report_path,
        "candidates": candidates_path,
        "samples": samples_tsv,
        "readme": project_dir / "README.md",
        "run": project_dir / "run.sh",
        staged_table_key: staged_table,
    }


def run_autopsy(assay: AssaySpec, out_dir: Path) -> dict[str, Path]:
    native = find_native_cli()
    out_dir.mkdir(parents=True, exist_ok=True)
    artifacts: dict[str, Path] = {
        "autopsy": out_dir,
        "summary": out_dir / "autopsy_summary.json",
        "findings": out_dir / "findings.tsv",
    }
    findings: list[dict[str, str]] = []

    if assay.mode == "count":
        _autopsy_count(assay, native, out_dir, findings, artifacts)
    elif assay.mode == "demux":
        _autopsy_demux(assay, native, out_dir, findings, artifacts)
    else:
        _autopsy_pair(assay, native, out_dir, findings, artifacts)

    _write_findings(artifacts["findings"], findings)
    summary = {
        "mode": assay.mode,
        "assay_type": assay.assay_type,
        "findings_count": len(findings),
        "artifacts": {key: str(value) for key, value in artifacts.items()},
    }
    artifacts["summary"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return artifacts


@dataclass(frozen=True)
class TargetSet:
    sequences: list[str]
    lengths: list[int]


def _read_target_sequences(path: Path) -> TargetSet:
    if not path.exists():
        raise AssaySpecError(f"target table does not exist: {path}")
    sequences: list[str] = []
    with path.open("r", encoding="utf-8") as fh:
        first_data = True
        seq_col = 1
        for raw in fh:
            line = raw.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            delim = "," if "," in line and "\t" not in line else "\t"
            row = next(csv.reader([line], delimiter=delim))
            if first_data:
                header = {name.strip().lower(): i for i, name in enumerate(row)}
                for name in ["grna.sequence", "target_seq", "sequence", "seq", "barcode_seq", "guide_seq", "sgrna.sequence"]:
                    if name in header:
                        seq_col = header[name]
                        first_data = False
                        break
                if not first_data:
                    continue
            first_data = False
            if len(row) == 1:
                seq = row[0].strip().upper()
            elif seq_col < len(row):
                seq = row[seq_col].strip().upper()
            else:
                raise AssaySpecError(f"target table row does not contain sequence column: {path}")
            if seq:
                sequences.append(seq)
    if not sequences:
        raise AssaySpecError(f"target table contains no sequences: {path}")
    return TargetSet(sequences=sequences, lengths=sorted(set(len(seq) for seq in sequences)))


def _read_pooled_fastq_sequences(paths: Sequence[Path], *, max_reads: int) -> list[str]:
    sequences: list[str] = []
    for path in paths:
        if len(sequences) >= max_reads:
            break
        remaining = max_reads - len(sequences)
        sequences.extend(_read_fastq_sequences(path, max_reads=remaining))
    return sequences


def _read_fastq_sequences(path: Path, *, max_reads: int) -> list[str]:
    seqs: list[str] = []
    opener = gzip.open if str(path).endswith(".gz") else Path.open
    with opener(path, "rt", encoding="utf-8") as fh:
        while len(seqs) < max_reads:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not seq or not plus or not qual:
                raise AssaySpecError(f"truncated FASTQ record in {path}")
            if not header.startswith("@") or not plus.startswith("+"):
                raise AssaySpecError(f"invalid FASTQ record in {path}")
            seqs.append(seq.strip().upper())
    return seqs


def _score_windows(reads: Sequence[str], targets: Sequence[str], *, max_start: int) -> list[dict[str, Any]]:
    lengths = sorted(set(len(seq) for seq in targets))
    matcher_by_len = {length: Matcher([seq for seq in targets if len(seq) == length]) for length in lengths}
    candidates: list[dict[str, Any]] = []
    try:
        for length in lengths:
            matcher = matcher_by_len[length]
            upper = min(max_start, max((len(seq) - length for seq in reads), default=0))
            for start in range(upper + 1):
                observed = [seq[start : start + length] for seq in reads if start + length <= len(seq)]
                invalid = len(reads) - len(observed)
                if observed:
                    results = matcher.assign(observed, k=1)
                else:
                    results = []
                unique = sum(1 for result in results if result.status == MATCH_UNIQUE)
                exact = sum(1 for result in results if result.status == MATCH_UNIQUE and result.best_distance == 0)
                ambiguous = sum(1 for result in results if result.status == MATCH_AMBIGUOUS)
                no_match = sum(1 for result in results if result.status == MATCH_NONE)
                total = len(reads)
                valid = len(observed)
                assignment_rate = unique / total if total else 0.0
                exact_rate = exact / total if total else 0.0
                ambiguous_rate = ambiguous / total if total else 0.0
                no_match_rate = no_match / total if total else 0.0
                invalid_rate = invalid / total if total else 0.0
                score = assignment_rate - ambiguous_rate - invalid_rate
                candidates.append(
                    {
                        "start": start,
                        "length": length,
                        "sampled_reads": total,
                        "valid_reads": valid,
                        "unique": unique,
                        "exact": exact,
                        "ambiguous": ambiguous,
                        "no_match": no_match,
                        "invalid": invalid,
                        "assignment_rate": round(assignment_rate, 8),
                        "exact_rate": round(exact_rate, 8),
                        "ambiguous_rate": round(ambiguous_rate, 8),
                        "no_match_rate": round(no_match_rate, 8),
                        "invalid_rate": round(invalid_rate, 8),
                        "score": round(score, 8),
                    }
                )
    finally:
        for matcher in matcher_by_len.values():
            matcher.close()
    candidates.sort(key=lambda item: (-float(item["score"]), int(item["start"]), int(item["length"])))
    if len(candidates) >= 2:
        best = candidates[0]
        second = candidates[1]
        best["score_margin"] = round(float(best["score"]) - float(second["score"]), 8)
    elif candidates:
        candidates[0]["score_margin"] = float(candidates[0]["score"])
    return candidates


def _choose_candidate(candidates: list[dict[str, Any]]) -> tuple[dict[str, Any], str, list[str]]:
    if not candidates:
        raise AssaySpecError("inference found no candidate windows")
    chosen = dict(candidates[0])
    chosen.setdefault("score_margin", float(chosen["score"]))
    warnings: list[str] = []
    if float(chosen["assignment_rate"]) < 0.80:
        warnings.append("best candidate assignment_rate is below 0.80")
    if float(chosen.get("score_margin", 0.0)) < 0.10:
        warnings.append("best candidate is not well separated from the next candidate")
    status = "draft" if warnings else "ready"
    return chosen, status, warnings


def _write_inferred_count_spec(out: Path, status: str, assay_type: str, targets: Path, reads: Path, sample_id: str,
                               chosen: Mapping[str, Any]) -> None:
    _require_safe_identifier(sample_id, "sample_id")
    format_name = "mageck" if assay_type == "crispr" else "dotmatch"
    command = "crispr" if assay_type == "crispr" else assay_type
    _write_text_file(
        out,
        f"""schema_version = 1
status = {_toml_string(status)}
mode = "count"
assay_type = {_toml_string(assay_type)}
targets = {_toml_string(str(targets))}

[[samples]]
id = {_toml_string(sample_id)}
fastq = {_toml_string(str(reads))}

[run]
out_dir = {_toml_string(f"{out.with_suffix('').name}_out")}
threads = 1

[extract]
start = {chosen["start"]}
length = {chosen["length"]}

[assignment]
k = 1
metric = "hamming"
ambiguity_policy = "radius"
ambiguous = "discard"

[outputs]
format = "{format_name}"
assignments = true
ambiguous = true
unmatched = true
""",
    )
    _ = command


def _write_inferred_demux_spec(out: Path, status: str, assay_type: str, barcodes: Path, reads: Path,
                               chosen: Mapping[str, Any]) -> None:
    _write_text_file(
        out,
        f"""schema_version = 1
status = {_toml_string(status)}
mode = "demux"
assay_type = {_toml_string(assay_type)}
barcodes = {_toml_string(str(barcodes))}
reads = {_toml_string(str(reads))}

[run]
out_dir = {_toml_string(f"{out.with_suffix('').name}_out")}

[extract]
start = {chosen["start"]}
length = {chosen["length"]}

[assignment]
k = 1
metric = "hamming"
ambiguity_policy = "radius"

[outputs]
assignments = true
ambiguous = true
unmatched = true
""",
    )


def _write_inferred_pair_spec(out: Path, status: str, assay_type: str, left_targets: Path, right_targets: Path,
                              reads: Path, left: Mapping[str, Any], right: Mapping[str, Any]) -> None:
    _write_text_file(
        out,
        f"""schema_version = 1
status = {_toml_string(status)}
mode = "pair-count"
assay_type = {_toml_string(assay_type)}
left_targets = {_toml_string(str(left_targets))}
right_targets = {_toml_string(str(right_targets))}
reads = {_toml_string(str(reads))}

[run]
out_dir = {_toml_string(f"{out.with_suffix('').name}_out")}

[left]
start = {left["start"]}
length = {left["length"]}

[right]
start = {right["start"]}
length = {right["length"]}

[assignment]
k = 1
metric = "hamming"
ambiguity_policy = "radius"

[outputs]
assignments = true
""",
    )


def _discover_fastqs(reads_dir: Path) -> list[Path]:
    fastqs: list[Path] = []
    for path in sorted(reads_dir.iterdir()):
        if not path.is_file():
            continue
        lower = path.name.lower()
        if any(lower.endswith(suffix) for suffix in FASTQ_SUFFIXES):
            fastqs.append(path)
    return fastqs


def _sample_id_from_fastq(path: Path) -> str:
    stem = path.name
    lower = stem.lower()
    for suffix in FASTQ_SUFFIXES:
        if lower.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    sample_id = re.sub(r"[^A-Za-z0-9_.-]", "_", stem).strip("._-")
    if not sample_id:
        raise AssaySpecError(f"could not derive sample_id from FASTQ filename: {path.name}")
    _require_safe_identifier(sample_id, "sample_id")
    return sample_id


def _unique_sample_id(sample_id: str, used: set[str]) -> str:
    if sample_id not in used:
        used.add(sample_id)
        return sample_id
    index = 2
    while True:
        candidate = f"{sample_id}_{index}"
        if candidate not in used:
            used.add(candidate)
            return candidate
        index += 1


def _stage_input_file(source: Path, dest: Path, *, copy: bool) -> None:
    if dest.exists():
        raise AssaySpecError(f"refusing to overwrite staged input: {dest}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    if copy:
        shutil.copy2(source, dest)
        return
    dest.symlink_to(source.resolve())


def _write_scaffolded_count_spec(
    out: Path,
    *,
    project_dir: Path,
    status: str,
    assay_type: str,
    targets: Path,
    samples: Sequence[tuple[str, Path, Path]],
    chosen: Mapping[str, Any],
    output_format: str,
    threads: int,
) -> None:
    sample_blocks = []
    for sample_id, staged_fastq, _source in samples:
        _require_safe_identifier(sample_id, "sample_id")
        rel_fastq = staged_fastq.relative_to(project_dir).as_posix()
        sample_blocks.append(
            f"[[samples]]\nid = {_toml_string(sample_id)}\nfastq = {_toml_string(rel_fastq)}\n"
        )
    rel_targets = targets.relative_to(project_dir).as_posix()
    _write_text_file(
        out,
        f"""schema_version = 1
status = {_toml_string(status)}
mode = "count"
assay_type = {_toml_string(assay_type)}
targets = {_toml_string(rel_targets)}
inference_report = "inference_report.json"

{"".join(sample_blocks)}
[run]
out_dir = "assay_out"
threads = {threads}

[extract]
start = {chosen["start"]}
length = {chosen["length"]}

[assignment]
k = 1
metric = "hamming"
ambiguity_policy = "radius"
ambiguous = "discard"

[outputs]
format = "{output_format}"
assignments = true
ambiguous = true
unmatched = true
""",
    )


def _write_scaffolded_demux_spec(
    out: Path,
    *,
    project_dir: Path,
    status: str,
    assay_type: str,
    barcodes: Path,
    reads: Path,
    chosen: Mapping[str, Any],
) -> None:
    _write_text_file(
        out,
        f"""schema_version = 1
status = {_toml_string(status)}
mode = "demux"
assay_type = {_toml_string(assay_type)}
barcodes = {_toml_string(barcodes.relative_to(project_dir).as_posix())}
reads = {_toml_string(reads.relative_to(project_dir).as_posix())}
inference_report = "inference_report.json"

[run]
out_dir = "assay_out"

[extract]
start = {chosen["start"]}
length = {chosen["length"]}

[assignment]
k = 1
metric = "hamming"
ambiguity_policy = "radius"

[outputs]
assignments = true
ambiguous = true
unmatched = true
""",
    )


def _write_generated_samples_tsv(path: Path, samples: Sequence[tuple[str, Path, Path]], project_dir: Path) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "fastq", "source_fastq"])
        for sample_id, staged_fastq, source in samples:
            writer.writerow([sample_id, staged_fastq.relative_to(project_dir).as_posix(), str(source)])


def _write_scaffold_readme(project_dir: Path, *, template: str, status: str, warnings: Sequence[str]) -> None:
    warning_lines = [f"- {warning}" for warning in warnings] or ["- none"]
    promote = ""
    if status == "draft":
        promote = (
            "\n## Promote To Ready\n\n"
            "This scaffold wrote `status = \"draft\"` because inference confidence was low. "
            "Review `inference_report.json` and `inference_candidates.tsv`, adjust `[extract]` if needed, "
            "then set `status = \"ready\"` in `assay.toml` before production counting.\n"
        )
    text = f"""# DotMatch Assay Project

Template: `{template}`

Status: `{status}`

## What Was Generated

- `assay.toml` — AssaySpec for check/plan/run
- `inference_report.json` — chosen extract window and warnings
- `inference_candidates.tsv` — ranked window candidates
- `samples.generated.tsv` — sample_id to FASTQ mapping
- `inputs/` — staged target or barcode table
- `reads/` — staged FASTQ inputs (copied by default for a self-contained project)
- `run.sh` — runs `dotmatch assay start` (check, then run) and prints the reliability verdict

## Inference Warnings

{chr(10).join(warning_lines)}
{promote}
## Recommended Workflow

```bash
./run.sh
```

`run.sh` runs `dotmatch assay start assay.toml`: preflight `check`, production `run`, then prints
`reliability_report.html` and any suggested `assay_fixes.tsv` edits.

Optional dry-run steps:

```bash
dotmatch assay check assay.toml
dotmatch assay start --check-only assay.toml
dotmatch assay plan assay.toml
```

After the run, open `assay_out/reliability_report.html` first. If the verdict is not `passed`,
apply the suggested edits in `assay_out/assay_fixes.tsv`, then rerun `./run.sh`.
"""
    (project_dir / "README.md").write_text(text, encoding="utf-8")


def _detect_dotmatch_launcher() -> list[str]:
    argv0 = Path(sys.argv[0]).resolve()
    if argv0.name == "cli.py" and argv0.parent.name == "dotmatch":
        return [sys.executable, "-m", "dotmatch.cli"]
    if argv0.name in {"dotmatch", "quickdna"}:
        return [str(argv0)]
    which = shutil.which("dotmatch")
    if which:
        return [which]
    return [sys.executable, "-m", "dotmatch.cli"]


def _launcher_uses_python_module(launcher: Sequence[str]) -> bool:
    return len(launcher) >= 3 and launcher[1] == "-m" and launcher[2] == "dotmatch.cli"


def _detect_pythonpath_for_scaffold() -> str:
    # Embed only the absolute checkout/package root; run.sh appends runtime PYTHONPATH.
    return str(Path(__file__).resolve().parents[1])


def _detect_native_cli_for_scaffold() -> str | None:
    try:
        return str(find_native_cli())
    except FileNotFoundError:
        return None


def _write_scaffold_run_script(
    project_dir: Path,
    *,
    status: str,
    launcher: Sequence[str],
    native_cli: str | None,
) -> None:
    launcher_literal = " ".join(shlex.quote(part) for part in launcher)
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        'ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        'cd "$ROOT"',
        "",
        f"DOTMATCH_LAUNCHER=({launcher_literal})",
    ]
    if _launcher_uses_python_module(launcher):
        lines.extend(
            [
                f'export PYTHONPATH={shlex.quote(_detect_pythonpath_for_scaffold())}${{PYTHONPATH:+:$PYTHONPATH}}',
                "",
            ]
        )
    else:
        lines.append("")
    if native_cli:
        lines.append(f'export DOTMATCH_NATIVE_CLI={shlex.quote(native_cli)}')
    lines.extend(
        [
            'if [[ -z "${DOTMATCH_NATIVE_CLI:-}" || ! -x "${DOTMATCH_NATIVE_CLI}" ]]; then',
            "  unset DOTMATCH_NATIVE_CLI",
            '  for candidate in "${ROOT}/../dotmatch" "${ROOT}/../../dotmatch" "${ROOT}/../../../dotmatch"; do',
            '    if [[ -x "$candidate" ]]; then',
            '      export DOTMATCH_NATIVE_CLI="$candidate"',
            "      break",
            "    fi",
            "  done",
            "fi",
            "",
            '_dotmatch_ready() {',
            '  "${DOTMATCH_LAUNCHER[@]}" assay --help >/dev/null 2>&1',
            "}",
            "",
            "if ! _dotmatch_ready; then",
            '  if command -v dotmatch >/dev/null 2>&1; then',
            '    DOTMATCH_LAUNCHER=(dotmatch)',
            '  elif python3 -c "import dotmatch.cli" >/dev/null 2>&1; then',
            '    DOTMATCH_LAUNCHER=(python3 -m dotmatch.cli)',
            "  else",
            '    for pyroot in "${ROOT}/../../python" "${ROOT}/../python" "${ROOT}/../../../python"; do',
            '      if [[ -f "${pyroot}/dotmatch/cli.py" ]]; then',
            '        export PYTHONPATH="${pyroot}${PYTHONPATH:+:${PYTHONPATH}}"',
            '        if python3 -c "import dotmatch.cli" >/dev/null 2>&1; then',
            '          DOTMATCH_LAUNCHER=(python3 -m dotmatch.cli)',
            "          break",
            "        fi",
            "      fi",
            "    done",
            "  fi",
            '  if ! _dotmatch_ready; then',
            '    echo "dotmatch not found: install dotmatch (pip/conda) or rerun assay new from a source checkout" >&2',
            "    exit 127",
            "  fi",
            "fi",
            "",
        ]
    )
    if status == "draft":
        lines.extend(
            [
                'if grep -qE \'^status = "draft"\' assay.toml; then',
                '  echo "Draft assay.toml: ./run.sh will stop at preflight until you promote status to ready." >&2',
                '  echo "  1. Review inference_report.json and assay_out/reliability_report.html" >&2',
                '  echo "  2. Apply suggestions in assay_out/assay_fixes.tsv" >&2',
                '  echo "  3. Set status = \\"ready\\" in assay.toml" >&2',
                "fi",
            ]
        )
    lines.extend(
        [
            'set +e',
            '"${DOTMATCH_LAUNCHER[@]}" assay start assay.toml "$@"',
            "rc=$?",
            "set -e",
            "",
            'if [[ "$rc" -eq 0 ]]; then',
            '  echo "done: assay_out/reliability_report.html (passed)" >&2',
            'elif [[ "$rc" -eq 1 ]]; then',
            '  echo "done: assay_out/reliability_report.html (needs review)" >&2',
            "else",
            '  echo "done: assay_out/reliability_report.html (failed or blocked)" >&2',
            "fi",
            'exit "$rc"',
        ]
    )
    run_path = project_dir / "run.sh"
    run_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    run_path.chmod(0o755)


def _write_text_file(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _toml_string(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def _write_candidates_tsv(path: Path, report_data: Mapping[str, Any]) -> None:
    rows: list[tuple[str, Mapping[str, Any]]] = []
    if "candidates" in report_data:
        rows.extend(("candidate", item) for item in report_data["candidates"])
    else:
        rows.extend(("left", item) for item in report_data["left"]["candidates"])
        rows.extend(("right", item) for item in report_data["right"]["candidates"])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write("side\tstart\tlength\tscore\tassignment_rate\texact_rate\tambiguous_rate\tno_match_rate\tinvalid_rate\n")
        for side, row in rows:
            fh.write(
                f"{side}\t{row['start']}\t{row['length']}\t{row['score']}\t{row['assignment_rate']}\t"
                f"{row['exact_rate']}\t{row['ambiguous_rate']}\t{row['no_match_rate']}\t{row['invalid_rate']}\n"
            )


def _autopsy_count(assay: AssaySpec, native: Path, out_dir: Path, findings: list[dict[str, str]],
                   artifacts: dict[str, Path]) -> None:
    extract = _table(assay.data, "extract")
    assignment = _table(assay.data, "assignment")
    audit_dir = out_dir / "audit"
    artifacts["audit"] = audit_dir
    audit_dir.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        _resolve_native(_audit_cmd(_spec_path(assay, "targets"), audit_dir, assay.k), native),
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    _add_audit_findings(audit_dir / "audit_summary.json", findings, "targets")
    for sample in _samples(assay.data):
        sample_id = str(sample["id"])
        top = out_dir / f"top_unmatched.{sample_id}.tsv"
        artifacts[f"top_unmatched_{sample_id}"] = top
        cmd = [
            str(native),
            "inspect-unmatched",
            "--targets",
            str(_spec_path(assay, "targets")),
            "--reads",
            str(_path_from_spec(assay.path, str(sample["fastq"]), allow_absolute=True, name="samples.fastq")),
            "--target-start",
            str(extract["start"]),
            "--target-length",
            str(extract["length"]),
            "--k",
            str(min(int(assignment.get("k", 1)), 1)),
            "--offset-window",
            str(max(5, int(assignment.get("auto_offset", 0)))),
            "--low-quality-threshold",
            "20",
            "--top",
            "100",
            "--out",
            str(top),
        ]
        subprocess.run(cmd, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _add_top_unmatched_findings(top, findings, sample_id)


def _autopsy_demux(assay: AssaySpec, native: Path, out_dir: Path, findings: list[dict[str, str]],
                   artifacts: dict[str, Path]) -> None:
    extract = _table(assay.data, "extract")
    assignment = _table(assay.data, "assignment")
    audit_dir = out_dir / "audit"
    artifacts["audit"] = audit_dir
    subprocess.run(
        _resolve_native(_audit_cmd(_spec_path(assay, "barcodes"), audit_dir, assay.k), native),
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    _add_audit_findings(audit_dir / "audit_summary.json", findings, "barcodes")
    length = extract["length"]
    if length == "auto":
        target_set = _read_target_sequences(_spec_path(assay, "barcodes"))
        length = target_set.lengths[0]
    top = out_dir / "top_unmatched.reads.tsv"
    artifacts["top_unmatched_reads"] = top
    cmd = [
        str(native), "inspect-unmatched", "--targets", str(_spec_path(assay, "barcodes")), "--reads", str(_spec_path(assay, "reads")),
        "--target-start", str(extract["start"]), "--target-length", str(length), "--k", str(min(int(assignment.get("k", 1)), 1)),
        "--offset-window", "5", "--top", "100", "--out", str(top),
    ]
    subprocess.run(cmd, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _add_top_unmatched_findings(top, findings, "reads")


def _autopsy_pair(assay: AssaySpec, native: Path, out_dir: Path, findings: list[dict[str, str]],
                  artifacts: dict[str, Path]) -> None:
    assignment = _table(assay.data, "assignment")
    for side, target_key, extract_key in [
        ("left", "left_targets", "left"),
        ("right", "right_targets", "right"),
    ]:
        audit_dir = out_dir / "audit" / side
        artifacts[f"{side}_audit"] = audit_dir
        subprocess.run(
            _resolve_native(_audit_cmd(_spec_path(assay, target_key), audit_dir, assay.k), native),
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        _add_audit_findings(audit_dir / "audit_summary.json", findings, side)
        extract = _table(assay.data, extract_key)
        top = out_dir / f"top_unmatched.{side}.tsv"
        artifacts[f"top_unmatched_{side}"] = top
        cmd = [
            str(native), "inspect-unmatched", "--targets", str(_spec_path(assay, target_key)), "--reads", str(_spec_path(assay, "reads")),
            "--target-start", str(extract["start"]), "--target-length", str(extract["length"]), "--k", str(min(int(assignment.get("k", 1)), 1)),
            "--offset-window", "5", "--top", "100", "--out", str(top),
        ]
        subprocess.run(cmd, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _add_top_unmatched_findings(top, findings, side)


def _add_audit_findings(summary_path: Path, findings: list[dict[str, str]], sample: str) -> None:
    if not summary_path.exists():
        return
    try:
        data = json.loads(summary_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return
    k = data.get("k", 1)
    if data.get(f"safe_at_k{k}") is False:
        findings.append(
            {
                "sample": sample,
                "finding": "unsafe_target_collisions",
                "severity": "warning",
                "evidence": f"safe_at_k{k}=false",
                "artifact": str(summary_path),
            }
        )


def _add_top_unmatched_findings(top_path: Path, findings: list[dict[str, str]], sample: str) -> None:
    if not top_path.exists():
        return
    text = top_path.read_text(encoding="utf-8")
    reason_counts: Counter[str] = Counter()
    for line in text.splitlines()[1:]:
        cols = line.split("\t")
        if len(cols) >= 7:
            reason_counts[cols[6]] += int(cols[1] or 0)
    reason_map = {
        "offset_shift_candidate": "wrong_offset",
        "reverse_complement_candidate": "reverse_complement_issue",
        "adapter_or_primer_candidate": "adapter_or_primer_candidate",
        "low_quality_candidate": "low_quality_candidate",
        "contains_N": "contains_n",
        "wrong_length": "wrong_length",
    }
    for reason, count in reason_counts.most_common():
        finding = reason_map.get(reason)
        if finding is None:
            continue
        severity = "error" if finding in {"wrong_offset", "wrong_length"} else "warning"
        findings.append(
            {
                "sample": sample,
                "finding": finding,
                "severity": severity,
                "evidence": f"{count} top unmatched reads: {reason}",
                "artifact": str(top_path),
            }
        )


def _write_findings(path: Path, findings: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "finding", "severity", "evidence", "artifact"])
        for finding in findings:
            writer.writerow(
                [
                    finding["sample"],
                    finding["finding"],
                    finding["severity"],
                    finding["evidence"],
                    finding["artifact"],
                ]
            )


def _autopsy_trigger_reasons(plan: AssayPlan) -> list[str]:
    if plan.spec.mode != "count":
        return []
    sample_qc = plan.artifacts.get("sample_qc")
    if sample_qc is None or not sample_qc.exists():
        return []
    reasons: list[str] = []
    thresholds = _reliability_thresholds(plan.spec)
    include_representation = plan.spec.assay_type == "crispr"
    with sample_qc.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = row.get("sample_id", "")
            try:
                failed = _failed_sample_qc_checks(
                    row,
                    thresholds,
                    include_representation=include_representation,
                )
            except AssaySpecError:
                continue
            for _finding_id, metric, observed, op, threshold in failed:
                reasons.append(f"{sample}: {metric} {op} {threshold} ({observed:.6f})")
    return reasons


def _float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _required_float(row: Mapping[str, Any], key: str) -> float:
    value = row.get(key)
    if value is None or value == "":
        raise AssaySpecError(f"sample_qc missing {key}")
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise AssaySpecError(f"sample_qc has invalid {key}: {value!r}") from exc


def _optional_float(row: Mapping[str, Any], key: str) -> float | None:
    value = row.get(key)
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _sample_qc_representation_metrics(row: Mapping[str, Any]) -> dict[str, float] | None:
    observed = _optional_float(row, "targets_observed")
    zero_count = _optional_float(row, "zero_count_targets")
    if observed is None or zero_count is None:
        return None
    guide_count = observed + zero_count
    if guide_count <= 0:
        return None
    metrics: dict[str, float] = {
        "coverage_fraction": observed / guide_count,
        "zero_count_fraction": zero_count / guide_count,
    }
    gini = _optional_float(row, "gini_index")
    if gini is not None:
        metrics["gini_index"] = gini
    top_fraction = _optional_float(row, "top_1pct_read_fraction")
    if top_fraction is not None:
        metrics["top_1pct_fraction"] = top_fraction
    return metrics


def _failed_sample_qc_checks(
    row: Mapping[str, Any],
    thresholds: Mapping[str, float],
    *,
    include_representation: bool,
) -> list[tuple[str, str, float, str, float]]:
    total = _required_float(row, "total_reads")
    invalid = _required_float(row, "invalid_reads")
    assignment_rate = _required_float(row, "assignment_rate")
    ambiguous_rate = _required_float(row, "ambiguous_rate")
    no_match_rate = _required_float(row, "no_match_rate")
    invalid_rate = invalid / total if total else 0.0
    checks: list[tuple[str, str, float, str, float]] = [
        ("assignment_rate_below_min", "assignment_rate", assignment_rate, "<", thresholds["min_assignment_rate"]),
        ("ambiguous_rate_above_max", "ambiguous_rate", ambiguous_rate, ">", thresholds["max_ambiguous_rate"]),
        ("unmatched_rate_above_max", "no_match_rate", no_match_rate, ">", thresholds["max_unmatched_rate"]),
        ("invalid_rate_above_max", "invalid_rate", invalid_rate, ">", thresholds["max_invalid_rate"]),
    ]
    if include_representation:
        representation = _sample_qc_representation_metrics(row)
        if representation is not None:
            checks.extend(
                [
                    (
                        "coverage_fraction_below_min",
                        "coverage_fraction",
                        representation["coverage_fraction"],
                        "<",
                        thresholds["min_coverage_fraction"],
                    ),
                    (
                        "zero_count_fraction_above_max",
                        "zero_count_fraction",
                        representation["zero_count_fraction"],
                        ">",
                        thresholds["max_zero_count_fraction"],
                    ),
                ]
            )
            gini = representation.get("gini_index")
            if gini is not None:
                checks.append(
                    ("gini_index_above_max", "gini_index", gini, ">", thresholds["max_gini_index"])
                )
            top_fraction = representation.get("top_1pct_fraction")
            if top_fraction is not None:
                checks.append(
                    (
                        "top_1pct_fraction_above_max",
                        "top_1pct_fraction",
                        top_fraction,
                        ">",
                        thresholds["max_top_1pct_fraction"],
                    )
                )
    failed: list[tuple[str, str, float, str, float]] = []
    for finding_id, metric, observed, op, threshold in checks:
        if observed < threshold if op == "<" else observed > threshold:
            failed.append((finding_id, metric, observed, op, threshold))
    return failed


def _compile_count(assay: AssaySpec, steps: list[PlanStep], artifacts: dict[str, Path], samples_path: Path) -> None:
    data = assay.data
    assignment = _table(data, "assignment")
    extract = _table(data, "extract")
    outputs = _table(data, "outputs")
    format_name = str(outputs.get("format", "mageck" if assay.assay_type == "crispr" else "dotmatch"))
    counts_name = "counts.mageck.tsv" if format_name == "mageck" else "counts.tsv"
    out_dir = assay.out_dir

    artifacts.update(
        {
            "samples": samples_path,
            "counts": out_dir / counts_name,
            "target_counts_long": out_dir / "target_counts.long.tsv",
            "sample_qc": out_dir / "sample_qc.tsv",
            "summary": out_dir / "summary.json",
            "report": out_dir / "report.html",
        }
    )
    cmd = [
        "dotmatch-native",
        "crispr-count" if assay.assay_type == "crispr" else "count",
        "--library" if assay.assay_type == "crispr" else "--targets",
        str(_spec_path(assay, "targets")),
        "--samples",
        str(samples_path),
        "--guide-start" if assay.assay_type == "crispr" else "--target-start",
        str(extract["start"]),
        "--guide-length" if assay.assay_type == "crispr" else "--target-length",
        str(extract["length"]),
        "--k",
        str(assignment.get("k", 1)),
        "--metric",
        str(assignment.get("metric", "levenshtein")),
        "--ambiguous",
        str(assignment.get("ambiguous", "discard")),
        "--threads",
        str(_table(data, "run").get("threads", 1)),
        "--out",
        str(artifacts["counts"]),
        "--summary",
        str(artifacts["summary"]),
        "--sample-qc",
        str(artifacts["sample_qc"]),
        "--target-counts-long",
        str(artifacts["target_counts_long"]),
        "--report",
        str(artifacts["report"]),
        "--report-audit-dir",
        str(artifacts["audit"]),
    ]
    if assay.assay_type != "crispr":
        cmd.extend(["--format", format_name])
    _add_assignment_options(cmd, assignment)
    backend = assay.backend
    backend_mode = str(backend.get("mode", "auto"))
    if backend_mode != "auto":
        cmd.extend(["--backend", backend_mode])
        if backend_mode == "gpu-metal-experimental":
            cmd.append("--metal-validate")
    elif not bool(backend.get("allow_gpu", True)):
        cmd.extend(["--backend", "cpu"])
    if outputs.get("assignments"):
        artifacts["assignments"] = out_dir / "assignments.tsv"
        cmd.extend(["--assignments", str(artifacts["assignments"])])
    if outputs.get("ambiguous"):
        artifacts["ambiguous"] = out_dir / "ambiguous.tsv"
        cmd.extend(["--ambiguous-out", str(artifacts["ambiguous"])])
    if outputs.get("unmatched"):
        artifacts["unmatched"] = out_dir / "unmatched.tsv"
        cmd.extend(["--unmatched-out", str(artifacts["unmatched"])])
    steps.append(PlanStep("run", cmd))

    if assay.assay_type == "crispr":
        artifacts["crispr_qc"] = out_dir / "crispr_qc.json"
        artifacts["crispr_qc_summary"] = out_dir / "crispr_qc.summary.tsv"
        artifacts["crispr_qc_report"] = out_dir / "crispr_qc.html"
        steps.append(
            PlanStep(
                "crispr-qc",
                [
                    "dotmatch",
                    "crispr-qc",
                    "--counts",
                    str(artifacts["counts"]),
                    "--sample-qc",
                    str(artifacts["sample_qc"]),
                    "--library",
                    str(_spec_path(assay, "targets")),
                    "--k",
                    str(assignment.get("k", 1)),
                    "--out",
                    str(artifacts["crispr_qc"]),
                    "--summary-tsv",
                    str(artifacts["crispr_qc_summary"]),
                    "--report",
                    str(artifacts["crispr_qc_report"]),
                ],
            )
        )

    first_sample = _samples(data)[0]
    validate = [
        "dotmatch-native",
        "validate",
        "--targets",
        str(_spec_path(assay, "targets")),
        "--reads",
        str(_path_from_spec(assay.path, str(first_sample["fastq"]), allow_absolute=True, name="samples.fastq")),
        "--target-start",
        str(extract["start"]),
        "--target-length",
        str(extract["length"]),
        "--k",
        str(assignment.get("k", 1)),
        "--metric",
        str(assignment.get("metric", "levenshtein")),
        "--sample",
        "100000",
    ]
    if "indel_window" in assignment:
        validate.extend(["--indel-window", str(assignment["indel_window"])])
    steps.append(PlanStep("validate", validate))


def _compile_demux(assay: AssaySpec, steps: list[PlanStep], artifacts: dict[str, Path]) -> None:
    data = assay.data
    assignment = _table(data, "assignment")
    extract = _table(data, "extract")
    outputs = _table(data, "outputs")
    out_dir = assay.out_dir
    artifacts.update({"demuxed": out_dir / "demuxed", "summary": out_dir / "summary.json"})
    cmd = [
        "dotmatch-native",
        "demux",
        "--barcodes",
        str(_spec_path(assay, "barcodes")),
        "--reads",
        str(_spec_path(assay, "reads")),
        "--barcode-start",
        str(extract["start"]),
        "--barcode-length",
        str(extract["length"]),
        "--k",
        str(assignment.get("k", 1)),
        "--metric",
        str(assignment.get("metric", "levenshtein")),
        "--out-dir",
        str(artifacts["demuxed"]),
        "--summary",
        str(artifacts["summary"]),
    ]
    _add_assignment_options(cmd, assignment)
    if outputs.get("assignments"):
        artifacts["assignments"] = out_dir / "assignments.tsv"
        cmd.extend(["--assignments", str(artifacts["assignments"])])
    if outputs.get("ambiguous"):
        artifacts["ambiguous"] = out_dir / "ambiguous.fastq"
        cmd.extend(["--ambiguous-out", str(artifacts["ambiguous"])])
    if outputs.get("unmatched"):
        artifacts["unmatched"] = out_dir / "unmatched.fastq"
        cmd.extend(["--unmatched-out", str(artifacts["unmatched"])])
    steps.append(PlanStep("run", cmd))


def _compile_pair(assay: AssaySpec, steps: list[PlanStep], artifacts: dict[str, Path]) -> None:
    data = assay.data
    assignment = _table(data, "assignment")
    outputs = _table(data, "outputs")
    left = _table(data, "left")
    right = _table(data, "right")
    out_dir = assay.out_dir
    artifacts.update({"pair_counts": out_dir / "pair_counts.tsv", "pair_summary": out_dir / "pair_summary.json"})
    cmd = [
        "dotmatch-native",
        "pair-count",
        "--left-targets",
        str(_spec_path(assay, "left_targets")),
        "--right-targets",
        str(_spec_path(assay, "right_targets")),
        "--reads",
        str(_spec_path(assay, "reads")),
        "--left-start",
        str(left["start"]),
        "--left-length",
        str(left["length"]),
        "--right-start",
        str(right["start"]),
        "--right-length",
        str(right["length"]),
        "--k",
        str(assignment.get("k", 1)),
        "--metric",
        str(assignment.get("metric", "levenshtein")),
        "--out",
        str(artifacts["pair_counts"]),
        "--summary",
        str(artifacts["pair_summary"]),
    ]
    cmd.extend(["--ambiguity-policy", str(assignment.get("ambiguity_policy", "radius"))])
    if outputs.get("assignments"):
        artifacts["pair_assignments"] = out_dir / "pair_assignments.tsv"
        cmd.extend(["--assignments", str(artifacts["pair_assignments"])])
    steps.append(PlanStep("run", cmd))


def _audit_cmd(targets: Path, out_dir: Path, k: int) -> list[str]:
    return ["dotmatch-native", "audit", "--targets", str(targets), "--k", str(k), "--audit-mode", "auto", "--out-dir", str(out_dir)]


def _add_assignment_options(cmd: list[str], assignment: Mapping[str, Any], *, include_ambiguity_policy: bool = True) -> None:
    if include_ambiguity_policy:
        cmd.extend(["--ambiguity-policy", str(assignment.get("ambiguity_policy", "radius"))])
    if "indel_window" in assignment:
        cmd.extend(["--indel-window", str(assignment["indel_window"])])
    if "max_correction_qual" in assignment:
        cmd.extend(["--max-correction-qual", str(assignment["max_correction_qual"])])
    if "auto_offset" in assignment:
        cmd.extend(["--auto-offset", str(assignment["auto_offset"])])


def _write_generated_files(plan: AssayPlan) -> None:
    samples_path = plan.generated_files.get("samples")
    if samples_path is None:
        return
    samples_path.parent.mkdir(parents=True, exist_ok=True)
    with samples_path.open("w", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "fastq"])
        for sample in _samples(plan.spec.data):
            writer.writerow(
                [
                    sample["id"],
                    _path_from_spec(plan.spec.path, str(sample["fastq"]), allow_absolute=True, name="samples.fastq"),
                ]
            )


def _write_normalized_spec(plan: AssayPlan) -> None:
    path = plan.artifacts["normalized_spec"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(plan.spec.data, fh, indent=2, sort_keys=True)
        fh.write("\n")


def _write_citation_artifacts(plan: AssayPlan, manifest: Mapping[str, Any]) -> None:
    metadata = _citation_metadata()
    _write_methods_md(plan.artifacts["methods"], plan, manifest, metadata)
    _write_citation_bib(plan.artifacts["citation_bib"], metadata)
    _write_software_versions_yml(plan.artifacts["software_versions"], plan, manifest)


def _project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _citation_metadata() -> dict[str, str]:
    cff = _project_root() / "CITATION.cff"
    metadata: dict[str, str] = {
        "title": "DotMatch: deterministic known-target short-DNA assignment for sequencing workflows",
        "family_names": "O'Toole",
        "given_names": "Donncha",
        "version": PYTHON_PACKAGE_VERSION,
        "doi": "10.5281/zenodo.20541628",
        "url": "https://github.com/dnncha/dotmatch",
        "year": "2026",
    }
    if not cff.exists():
        return metadata
    current_author = False
    for raw_line in cff.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line.startswith("title:"):
            metadata["title"] = _unquote_cff_value(line.split(":", 1)[1].strip())
        elif line.startswith("version:"):
            metadata["version"] = _unquote_cff_value(line.split(":", 1)[1].strip())
        elif line.startswith("doi:"):
            metadata["doi"] = _unquote_cff_value(line.split(":", 1)[1].strip())
        elif line.startswith("repository-code:"):
            metadata["url"] = _unquote_cff_value(line.split(":", 1)[1].strip())
        elif line.startswith("- given-names:") and "given_names" in metadata:
            current_author = True
            metadata["given_names"] = _unquote_cff_value(line.split(":", 1)[1].strip())
        elif current_author and line.startswith("family-names:"):
            metadata["family_names"] = _unquote_cff_value(line.split(":", 1)[1].strip())
            current_author = False
    return metadata


def _unquote_cff_value(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        return value[1:-1]
    return value


def _write_methods_md(path: Path, plan: AssayPlan, manifest: Mapping[str, Any], metadata: Mapping[str, str]) -> None:
    assignment = _table(plan.spec.data, "assignment")
    reliability = plan.spec.reliability
    backend = plan.spec.backend
    extract_lines = _methods_extract_lines(plan.spec)
    target_lines = _methods_target_lines(plan.spec)
    sample_lines = _methods_sample_lines(plan.spec)
    command_lines = _methods_command_lines(manifest)
    warning_lines = [str(item) for item in manifest.get("production_warnings", []) or []] + [
        str(item) for item in manifest.get("warnings", []) or []
    ]
    if not warning_lines:
        warning_lines = ["No AssaySpec production warnings were recorded."]

    text = "\n".join(
        [
            "# DotMatch Methods and Citation",
            "",
            "## Methods",
            "",
            f"DotMatch {PYTHON_PACKAGE_VERSION} was used for deterministic known-target short-DNA assignment.",
            f"The assay mode was `{plan.spec.mode}` and the assay type was `{plan.spec.assay_type}`.",
            "Reads were assigned only when exactly one known target was compatible under the configured policy; ambiguous reads were not silently counted.",
            "",
            "## Assignment Configuration",
            "",
            f"- Edit radius (`k`): `{plan.spec.k}`",
            f"- Metric: `{assignment.get('metric', 'levenshtein')}`",
            f"- Ambiguity policy: `{assignment.get('ambiguity_policy', 'radius')}`",
            f"- Ambiguous output policy: `{assignment.get('ambiguous', 'discard')}`",
            f"- Reliability profile: `{reliability.get('profile', 'production')}`",
            f"- Backend mode: `{backend.get('mode', 'auto')}`",
            f"- Native version: `{manifest.get('native_version', '')}`",
            "",
            "## Extraction Windows",
            "",
            *extract_lines,
            "",
            "## Inputs",
            "",
            *target_lines,
            *sample_lines,
            "",
            "## Commands",
            "",
            *command_lines,
            "",
            "## Warnings",
            "",
            *(f"- {warning}" for warning in warning_lines),
            "",
            "## Citation",
            "",
            f"Cite DotMatch as: {metadata.get('family_names', '')} {metadata.get('given_names', '')}. {metadata.get('title', '')}. Software release v{metadata.get('version', PYTHON_PACKAGE_VERSION)}. {metadata.get('url', '')}.",
            f"DOI: https://doi.org/{metadata.get('doi', '')}",
            "",
            "BibTeX is written to `CITATION.bib`; software versions are written to `software_versions.yml`.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _methods_extract_lines(assay: AssaySpec) -> list[str]:
    if assay.mode in {"count", "demux"}:
        extract = _table(assay.data, "extract")
        return [f"- Primary window: start `{extract.get('start')}`, length `{extract.get('length')}`"]
    left = _table(assay.data, "left")
    right = _table(assay.data, "right")
    return [
        f"- Left window: start `{left.get('start')}`, length `{left.get('length')}`",
        f"- Right window: start `{right.get('start')}`, length `{right.get('length')}`",
    ]


def _methods_target_lines(assay: AssaySpec) -> list[str]:
    keys = ["targets"] if assay.mode == "count" else ["barcodes"] if assay.mode == "demux" else ["left_targets", "right_targets"]
    return [f"- `{key}`: `{assay.data.get(key, '')}`" for key in keys]


def _methods_sample_lines(assay: AssaySpec) -> list[str]:
    if assay.mode == "count":
        return [
            f"- Sample `{sample.get('id', '')}` FASTQ: `{sample.get('fastq', '')}`"
            for sample in _samples(assay.data)
        ]
    return [f"- Reads: `{assay.data.get('reads', '')}`"]


def _methods_command_lines(manifest: Mapping[str, Any]) -> list[str]:
    commands = manifest.get("commands", []) or []
    if not commands:
        return ["- No native commands were recorded for this preflight artifact."]
    lines = []
    for command in commands:
        if not isinstance(command, Mapping):
            continue
        argv = command.get("argv", []) or []
        lines.append(f"- `{shlex.join(str(part) for part in argv)}`")
    return lines or ["- No native commands were recorded for this preflight artifact."]


def _write_citation_bib(path: Path, metadata: Mapping[str, str]) -> None:
    family = metadata.get("family_names", "O'Toole")
    given = metadata.get("given_names", "Donncha")
    title = metadata.get("title", "DotMatch")
    version = metadata.get("version", PYTHON_PACKAGE_VERSION)
    doi = metadata.get("doi", "")
    url = metadata.get("url", "https://github.com/dnncha/dotmatch")
    year = metadata.get("year", "2026")
    bib = f"""@software{{dotmatch,
  author = {{{_bibtex_escape(family)}, {_bibtex_escape(given)}}},
  title = {{{_bibtex_escape(title)}}},
  version = {{{_bibtex_escape(version)}}},
  year = {{{_bibtex_escape(year)}}},
  doi = {{{_bibtex_escape(doi)}}},
  url = {{{_bibtex_escape(url)}}}
}}
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(bib, encoding="utf-8")


def _bibtex_escape(value: str) -> str:
    return value.replace("\\", "\\\\").replace("{", "\\{").replace("}", "\\}")


def _write_software_versions_yml(path: Path, plan: AssayPlan, manifest: Mapping[str, Any]) -> None:
    versions = {
        "dotmatch_python": PYTHON_PACKAGE_VERSION,
        "dotmatch_native": str(manifest.get("native_version", "")),
        "python": sys.version.split()[0],
        "assayspec_schema": "1",
        "assayspec_mode": plan.spec.mode,
        "assayspec_assay_type": plan.spec.assay_type,
        "workflow_wrapper": "dotmatch assay",
        "report_tool": "dotmatch assayspec report",
    }
    lines = ["software:"]
    for key, value in versions.items():
        lines.append(f"  {key}: {_yaml_scalar(value)}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _yaml_scalar(value: object) -> str:
    text = str(value)
    return "'" + text.replace("'", "''") + "'"


def _write_manifest(plan: AssayPlan, manifest: Mapping[str, Any]) -> None:
    _write_citation_artifacts(plan, manifest)
    path = plan.artifacts["manifest"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")
    _write_manifest_summary(plan, manifest)
    _write_assay_report(plan, manifest)


def _write_preflight_reliability(plan: AssayPlan) -> dict[str, Any]:
    summary = _build_reliability_summary(plan, stage="preflight")
    _write_reliability_artifacts(plan, summary)
    return summary


def _publish_reliability_summary(
    plan: AssayPlan,
    summary: dict[str, Any],
    *,
    extra_findings: Sequence[Mapping[str, str]] = (),
    manifest: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    if extra_findings:
        summary = dict(summary)
        summary["findings"] = list(summary.get("findings", []) or []) + list(extra_findings)
        summary["finding_counts"] = _finding_counts(summary["findings"])
        summary["overall_status"] = _overall_reliability_status(summary["finding_counts"])
        fixes = _build_assay_fixes(plan, summary["findings"], manifest)
        summary["assay_fixes"] = fixes
        summary["findings"] = _apply_recommended_actions(summary["findings"], fixes)
    _write_reliability_artifacts(plan, summary)
    return summary


def _preflight_for_assay_start(plan: AssayPlan) -> str:
    manifest: dict[str, Any] = {
        "commands": [],
        "warnings": [],
        "production_warnings": [],
    }
    try:
        native = find_native_cli()
    except FileNotFoundError as exc:
        summary = _build_reliability_summary(plan, stage="preflight", manifest=manifest)
        summary = _publish_reliability_summary(
            plan,
            summary,
            extra_findings=[
                _finding(
                    "native_cli_missing",
                    "blocked",
                    "preflight",
                    "",
                    "native_cli",
                    "missing",
                    "executable",
                    str(exc),
                    "Build with make dotmatch, install a wheel with the bundled native executable, or set DOTMATCH_NATIVE_CLI.",
                    "",
                )
            ],
            manifest=manifest,
        )
        return str(summary["overall_status"])

    for step in plan.steps:
        if not step.name.startswith("audit"):
            continue
        Path(step.argv[-1]).parent.mkdir(parents=True, exist_ok=True)
        argv = _resolve_native(step.argv, native)
        result = subprocess.run(argv, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        manifest["commands"].append(
            {
                "name": step.name,
                "argv": argv,
                "exit_code": result.returncode,
                "stdout": result.stdout,
                "stderr": result.stderr,
            }
        )
        if result.returncode != 0:
            summary = _build_reliability_summary(plan, stage="preflight", manifest=manifest)
            summary = _publish_reliability_summary(plan, summary, manifest=manifest)
            return str(summary["overall_status"])
        _append_audit_warnings(plan, step, manifest, warn=not _audit_step_unsafe(plan, step))
        if _audit_step_unsafe(plan, step) and _blocks_on_unsafe_targets(plan.spec):
            manifest["production_warnings"].append(
                f"{step.name}: unsafe target audit at k={plan.spec.k}; blocked by production reliability profile"
            )
            summary = _build_reliability_summary(plan, stage="preflight", manifest=manifest)
            summary = _publish_reliability_summary(plan, summary, manifest=manifest)
            return str(summary["overall_status"])

    summary = _build_reliability_summary(plan, stage="preflight", manifest=manifest)
    summary = _publish_reliability_summary(plan, summary, manifest=manifest)
    return str(summary["overall_status"])


def _build_reliability_summary(plan: AssayPlan, *, stage: str, manifest: Mapping[str, Any] | None = None) -> dict[str, Any]:
    findings: list[dict[str, str]] = []
    findings.extend(_base_reliability_findings(plan, stage))
    if manifest is not None or stage != "preflight":
        findings.extend(_audit_reliability_findings(plan, stage))
    if manifest is not None:
        findings.extend(_command_reliability_findings(manifest, stage))
        if stage != "preflight":
            findings.extend(_sample_qc_reliability_findings(plan, stage))
            findings.extend(_count_run_reliability_findings(plan, stage))
            findings.extend(_crispr_qc_reliability_findings(plan, stage))
            findings.extend(_autopsy_reliability_findings(manifest, stage))
    if stage == "preflight":
        findings.append(
            _finding(
                "read_qc_unavailable",
                "info",
                "preflight",
                "",
                "sample_qc",
                "unavailable",
                "",
                "Read-dependent QC is unavailable during assay check.",
                "Run dotmatch assay run to evaluate assignment, ambiguous, unmatched, and invalid read rates.",
                "",
            )
        )
    counts = _finding_counts(findings)
    fixes = _build_assay_fixes(plan, findings, manifest)
    findings = _apply_recommended_actions(findings, fixes)
    summary = {
        "schema_version": 1,
        "stage": stage,
        "overall_status": _overall_reliability_status(counts),
        "mode": plan.spec.mode,
        "assay_type": plan.spec.assay_type,
        "spec_status": plan.spec.status,
        "profile": str(plan.spec.reliability["profile"]),
        "thresholds": _reliability_thresholds(plan.spec),
        "backend": _backend_summary(plan.spec, _read_json_artifact(plan.artifacts.get("summary")) if stage != "preflight" else None),
        "backend_optimizer": optimize_assay_backend(plan.spec),
        "evidence_boundary": _evidence_boundary(plan.spec),
        "findings": findings,
        "finding_counts": counts,
        "assay_fixes": fixes,
        "artifacts": {key: str(value) for key, value in plan.artifacts.items()},
        "commands": list((manifest or {}).get("commands", []) or []),
    }
    return summary


def _command_reliability_findings(manifest: Mapping[str, Any], stage: str) -> list[dict[str, str]]:
    findings: list[dict[str, str]] = []
    for command in manifest.get("commands", []) or []:
        if not isinstance(command, dict):
            continue
        exit_code = command.get("exit_code")
        if exit_code not in (0, None):
            findings.append(
                _finding(
                    "command_failed",
                    "error",
                    stage,
                    "",
                    str(command.get("name", "command")),
                    str(exit_code),
                    "0",
                    "A native command failed during the assay run.",
                    "Inspect command stderr/stdout in assay_manifest.json before trusting outputs.",
                    str(manifest.get("artifacts", {}).get("manifest", "")),
                )
            )
    return findings


def _sample_qc_reliability_findings(plan: AssayPlan, stage: str) -> list[dict[str, str]]:
    sample_qc = plan.artifacts.get("sample_qc")
    if sample_qc is None or not sample_qc.exists():
        return [
            _finding(
                "read_qc_unavailable",
                "info",
                stage,
                "",
                "sample_qc",
                "unavailable",
                "",
                "Read-dependent QC was not available for this assay mode.",
                "Use a count workflow with sample_qc.tsv to evaluate read-level reliability thresholds.",
                "",
            )
        ]
    findings: list[dict[str, str]] = []
    thresholds = _reliability_thresholds(plan.spec)
    severity = _threshold_severity(plan.spec)
    with sample_qc.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"sample_id", "assignment_rate", "ambiguous_rate", "no_match_rate", "total_reads", "invalid_reads"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            return [
                _finding(
                    "sample_qc_malformed",
                    "error",
                    stage,
                    "",
                    "sample_qc",
                    f"missing columns: {', '.join(missing)}",
                    "required columns present",
                    "sample_qc.tsv is missing required reliability columns.",
                    "Regenerate assay outputs before trusting reliability metrics.",
                    str(sample_qc),
                )
            ]
        for row in reader:
            sample = row.get("sample_id", "")
            try:
                failed = _failed_sample_qc_checks(
                    row,
                    thresholds,
                    include_representation=plan.spec.assay_type == "crispr",
                )
            except AssaySpecError as exc:
                findings.append(
                    _finding(
                        "sample_qc_malformed",
                        "error",
                        stage,
                        sample,
                        "sample_qc",
                        str(exc),
                        "numeric reliability metrics",
                        "sample_qc.tsv contains malformed reliability data.",
                        "Regenerate assay outputs before trusting reliability metrics.",
                        str(sample_qc),
                    )
                )
                continue
            for finding_id, metric, observed, op, threshold in failed:
                findings.append(
                    _finding(
                        finding_id,
                        severity,
                        stage,
                        sample,
                        metric,
                        f"{observed:.8f}",
                        f"{threshold:.8f}",
                        f"{sample or 'sample'} has {metric} {op} configured reliability threshold.",
                        "Review the assay window, target set, correction radius, and autopsy outputs before using counts.",
                        str(sample_qc),
                    )
                )
    return findings


def _read_json_artifact(path: Path | None) -> dict[str, Any] | None:
    if path is None or not path.exists():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    return data if isinstance(data, dict) else None


def _count_run_reliability_findings(plan: AssayPlan, stage: str) -> list[dict[str, str]]:
    if plan.spec.mode != "count":
        return []
    run_summary = _read_json_artifact(plan.artifacts.get("summary"))
    if run_summary is None:
        return []
    findings: list[dict[str, str]] = []
    summary_path = str(plan.artifacts.get("summary", ""))
    metal_validation = run_summary.get("metal_validation")
    if metal_validation == "failed":
        findings.append(
            _finding(
                "metal_validation_failed",
                "error",
                stage,
                "",
                "metal_validation",
                "failed",
                "passed",
                "Experimental Metal counting did not match the CPU authority checksum.",
                "Re-run with --backend cpu or fix Metal eligibility before trusting counts.",
                summary_path,
            )
        )
    backend_effective = run_summary.get("backend_effective")
    if backend_effective == "gpu-metal-experimental" and metal_validation != "passed":
        findings.append(
            _finding(
                "experimental_gpu_backend_without_validation",
                "warning",
                stage,
                "",
                "backend_effective",
                str(backend_effective),
                "cpu",
                "Counting used the experimental Metal backend without a recorded CPU validation pass.",
                "Enable --metal-validate or compare against a CPU shadow run before downstream analysis.",
                summary_path,
            )
        )
    return findings


def _crispr_qc_reliability_findings(plan: AssayPlan, stage: str) -> list[dict[str, str]]:
    if plan.spec.assay_type != "crispr":
        return []
    crispr_qc = _read_json_artifact(plan.artifacts.get("crispr_qc"))
    if crispr_qc is None:
        return []
    findings: list[dict[str, str]] = []
    severity = _threshold_severity(plan.spec)
    source = str(plan.artifacts.get("crispr_qc", ""))
    thresholds = _reliability_thresholds(plan.spec)
    correlation_rows = crispr_qc.get("sample_correlations", crispr_qc.get("replicates", [])) or []
    for pair in correlation_rows:
        if not isinstance(pair, dict):
            continue
        pearson = pair.get("pearson_log2_count_plus_1")
        if not isinstance(pearson, (int, float)):
            continue
        threshold = thresholds["min_pairwise_sample_pearson"]
        if pearson >= threshold:
            continue
        sample_a = str(pair.get("sample_a", ""))
        sample_b = str(pair.get("sample_b", ""))
        scope = f"{sample_a}:{sample_b}" if sample_a and sample_b else ""
        findings.append(
            _finding(
                "pairwise_sample_correlation_below_min",
                severity,
                stage,
                scope,
                "pairwise_sample_pearson",
                f"{pearson:.8f}",
                f"{threshold:.8f}",
                f"Pairwise sample log2(count+1) Pearson correlation is {pearson:.3f}.",
                "Review replicate concordance and library representation before hit calling.",
                source,
            )
        )
    for warning in crispr_qc.get("warnings", []) or []:
        if not isinstance(warning, dict):
            continue
        code = str(warning.get("code", ""))
        if code == "low_pairwise_sample_correlation":
            continue
        if warning.get("scope") in {"inputs", "library"} or code in {
            "guide_collision",
            "non_acgt_guides",
            "library_collision_audit_radius",
            "sample_qc_not_provided",
            "library_not_provided",
        }:
            findings.append(
                _finding(
                    f"crispr_qc_{code}" if code else "crispr_qc_library_warning",
                    "warning" if code in {"guide_collision", "non_acgt_guides"} else severity,
                    stage,
                    str(warning.get("scope", "")),
                    code or "warning",
                    "review",
                    "pass",
                    str(warning.get("message", "CRISPR QC reported a library or input review warning.")),
                    "Inspect crispr_qc.json and the guide library before production counting.",
                    source,
                )
            )
    if crispr_qc.get("status") == "review" and not findings:
        findings.append(
            _finding(
                "crispr_qc_review",
                severity,
                stage,
                "",
                "status",
                "review",
                "pass",
                "CRISPR QC reported review status without per-sample detail in reliability findings.",
                "Inspect crispr_qc.json before downstream screen statistics.",
                source,
            )
        )
    return findings


def _audit_reliability_findings(plan: AssayPlan, stage: str) -> list[dict[str, str]]:
    findings: list[dict[str, str]] = []
    for key in ["audit", "left_audit", "right_audit"]:
        audit_dir = plan.artifacts.get(key)
        if audit_dir is None:
            continue
        summary = Path(audit_dir) / "audit_summary.json"
        if not summary.exists():
            continue
        try:
            data = json.loads(summary.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            findings.append(
                _finding(
                    "audit_summary_malformed",
                    "error",
                    stage,
                    "",
                    key,
                    "malformed",
                    "valid JSON",
                    "Target audit summary could not be parsed.",
                    "Re-run assay audit and inspect audit_summary.json.",
                    str(summary),
                )
            )
            continue
        safe_key = f"safe_at_k{plan.spec.k}"
        if data.get(safe_key) is False:
            severity = "blocked" if _blocks_on_unsafe_targets(plan.spec) else _threshold_severity(plan.spec)
            findings.append(
                _finding(
                    "unsafe_targets",
                    severity,
                    stage,
                    "",
                    safe_key,
                    "false",
                    "true",
                    f"Target audit reports {safe_key}=false for {key}.",
                    "Lower correction radius, use k=0, or redesign/fix colliding targets before production use.",
                    str(summary),
                )
            )
    return findings


def _blocks_on_unsafe_targets(assay: AssaySpec) -> bool:
    reliability = assay.reliability
    return str(reliability["profile"]) == "production" and bool(reliability["fail_on_unsafe_targets"])


def _blocks_on_draft_inference(assay: AssaySpec) -> bool:
    reliability = assay.reliability
    return str(reliability["profile"]) == "production" and bool(reliability["fail_on_draft_inference"])


def _audit_step_unsafe(plan: AssayPlan, step: PlanStep) -> bool:
    summary = Path(step.argv[-1]) / "audit_summary.json"
    if not summary.exists():
        return False
    try:
        data = json.loads(summary.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return False
    return data.get(f"safe_at_k{plan.spec.k}") is False


def _autopsy_reliability_findings(manifest: Mapping[str, Any], stage: str) -> list[dict[str, str]]:
    if not manifest.get("autopsy_triggered"):
        return []
    findings = [
        _finding(
            "autopsy_triggered",
            "warning",
            stage,
            "",
            "autopsy",
            "triggered",
            "not_triggered",
            "Automatic autopsy was triggered by conservative assay QC thresholds.",
            "Review autopsy findings before using the run for downstream interpretation.",
            str((manifest.get("autopsy_artifacts", {}) or {}).get("findings", "")),
        )
    ]
    autopsy_findings = (manifest.get("autopsy_artifacts", {}) or {}).get("findings")
    if autopsy_findings and Path(str(autopsy_findings)).exists():
        with Path(str(autopsy_findings)).open("r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                findings.append(
                    _finding(
                        f"autopsy_{row.get('finding', 'finding')}",
                        str(row.get("severity", "warning") or "warning"),
                        stage,
                        str(row.get("sample", "")),
                        "autopsy",
                        str(row.get("finding", "")),
                        "",
                        str(row.get("evidence", "")),
                        "Review the linked autopsy artifact before using the run.",
                        str(row.get("artifact", autopsy_findings)),
                    )
                )
    return findings


def _threshold_severity(assay: AssaySpec) -> str:
    return "error" if str(assay.reliability["profile"]) == "production" else "warning"


def _base_reliability_findings(plan: AssayPlan, stage: str) -> list[dict[str, str]]:
    findings: list[dict[str, str]] = []
    reliability = plan.spec.reliability
    if plan.spec.status == "draft":
        severity = "blocked" if _blocks_on_draft_inference(plan.spec) else _threshold_severity(plan.spec)
        findings.append(
            _finding(
                "draft_assayspec",
                severity,
                stage,
                "",
                "status",
                plan.spec.status,
                "ready",
                "Production reliability profile refuses draft inferred AssaySpec files.",
                "Review the inference report and promote status to ready before production runs.",
                str(plan.spec.path),
            )
        )
    evidence = _evidence_boundary(plan.spec)
    if reliability["require_public_evidence_boundary"] and evidence["status"] == "missing":
        findings.append(
            _finding(
                "evidence_boundary_missing",
                "error",
                stage,
                "",
                "assay_type",
                plan.spec.assay_type,
                "documented public evidence boundary",
                "No assay evidence boundary is recorded for this assay type.",
                "Add assay evidence metadata before making public claims for this workflow.",
                "docs/assay-evidence.json",
            )
        )
    elif reliability["require_public_evidence_boundary"] and evidence["status"] != "supported":
        findings.append(
            _finding(
                "evidence_boundary_not_supported",
                _threshold_severity(plan.spec),
                stage,
                "",
                "assay_type",
                plan.spec.assay_type,
                "supported public evidence boundary",
                f"Assay evidence boundary is {evidence['status']!r}, not supported.",
                "Treat this run as smoke or gated evidence until public comparator validation is recorded.",
                "docs/assay-evidence.json",
            )
        )
    backend = _backend_summary(plan.spec)
    if backend["mode"] == "gpu-metal-experimental":
        findings.append(
            _finding(
                "experimental_gpu_forced",
                "warning",
                stage,
                "",
                "backend.mode",
                "gpu-metal-experimental",
                "cpu or auto for production",
                "The Metal GPU backend is experimental and CPU remains the assignment authority.",
                "Use backend.mode = \"auto\" or \"cpu\" for production evidence until assay-specific GPU gates are promoted.",
                "",
            )
        )
    return findings


def _reliability_thresholds(assay: AssaySpec) -> dict[str, float]:
    reliability = assay.reliability
    return {
        "min_assignment_rate": float(reliability["min_assignment_rate"]),
        "max_ambiguous_rate": float(reliability["max_ambiguous_rate"]),
        "max_unmatched_rate": float(reliability["max_unmatched_rate"]),
        "max_invalid_rate": float(reliability["max_invalid_rate"]),
        "min_coverage_fraction": float(reliability["min_coverage_fraction"]),
        "max_zero_count_fraction": float(reliability["max_zero_count_fraction"]),
        "max_gini_index": float(reliability["max_gini_index"]),
        "max_top_1pct_fraction": float(reliability["max_top_1pct_fraction"]),
        "min_pairwise_sample_pearson": float(reliability["min_pairwise_sample_pearson"]),
    }


def _backend_summary(assay: AssaySpec, run_summary: Mapping[str, Any] | None = None) -> dict[str, str]:
    backend = assay.backend
    mode = str(backend["mode"])
    if mode == "cpu":
        gpu_status = "disabled_by_mode"
    elif not bool(backend["allow_gpu"]):
        gpu_status = "disabled_by_config"
    elif _gpu_eligible(assay) and _gpu_public_evidence_validated(assay):
        gpu_status = "eligible_but_not_used"
    elif _gpu_eligible(assay):
        gpu_status = "compute_compatible_no_public_gpu_gate"
    else:
        gpu_status = "not_eligible"
    summary: dict[str, str] = {
        "mode": mode,
        "authority": "cpu",
        "selected": "cpu",
        "gpu_status": gpu_status,
    }
    if run_summary is not None:
        backend_effective = run_summary.get("backend_effective")
        if isinstance(backend_effective, str) and backend_effective:
            summary["effective"] = backend_effective
            if backend_effective != "cpu":
                summary["selected"] = backend_effective
        metal_validation = run_summary.get("metal_validation")
        if isinstance(metal_validation, str) and metal_validation:
            summary["metal_validation"] = metal_validation
    return summary


def optimize_assay_backend(assay: AssaySpec) -> dict[str, Any]:
    features, reason_codes = _backend_optimizer_features(assay)
    speed_model = _gpu_speed_model(features)
    compute_compatible = not any(code.endswith("_not_gpu_supported") or code.endswith("_not_gpu_packable") for code in reason_codes)
    public_gate = _gpu_public_evidence_validated(assay)
    backend = assay.backend
    allow_gpu = bool(backend["allow_gpu"]) and str(backend["mode"]) != "cpu"
    cpu_strategy, route_reasons = _cpu_route_strategy(features, reason_codes)
    benchmark_confidence = _benchmark_confidence(features, reason_codes, compute_compatible, public_gate, allow_gpu)

    if not allow_gpu:
        candidate_backend = "cpu"
        recommendation = "cpu_required"
        expected_speedup_band = "1x"
        if str(backend["mode"]) == "cpu":
            reason_codes.append("gpu_disabled_by_mode")
            route_reasons.append("gpu_disabled_by_mode")
        else:
            reason_codes.append("gpu_disabled_by_config")
            route_reasons.append("gpu_disabled_by_config")
    elif compute_compatible and public_gate:
        candidate_backend = "gpu-metal-experimental"
        recommendation = "gpu_candidate_requires_cpu_validation"
        expected_speedup_band = _gpu_expected_speedup_band(speed_model["estimated_total_speedup"])
        reason_codes.append("public_gpu_gate_validated")
        route_reasons.append("gpu_candidate_public_gate")
    elif compute_compatible:
        candidate_backend = "gpu-metal-experimental"
        recommendation = "gpu_candidate_gated"
        expected_speedup_band = "unknown_until_public_gate"
        reason_codes.append("compute_compatible_no_public_gpu_gate")
        route_reasons.append("gpu_candidate_without_public_gate")
    else:
        candidate_backend = "cpu"
        recommendation = "cpu_required"
        expected_speedup_band = "1x"
        route_reasons.append("gpu_ineligible_cpu_only")

    accuracy_gates = [
        "cpu_assignment_authority",
        "cpu_count_checksum_required",
        "zero_mismatch_required_before_speed_claim",
    ]
    diagnostic_constraints = [
        "cpu_remains_assignment_authority",
        "cpu_count_checksum_required",
        "gpu_candidate_requires_zero_mismatch_diagnostic",
        "benchmark_priors_are_route_metadata_only",
    ]
    return {
        "schema_version": 1,
        "optimizer": "local_benchmark_informed_scorer_v1",
        "authority": "cpu",
        "selected_backend": "cpu",
        "candidate_backend": candidate_backend,
        "recommendation": recommendation,
        "expected_speedup_band": expected_speedup_band,
        "estimated_total_speedup": speed_model["estimated_total_speedup"],
        "speed_model": speed_model,
        "cpu_strategy": cpu_strategy,
        "thread_hint": _thread_hint(assay, features),
        "benchmark_prior_count": int(speed_model["training_rows"]),
        "benchmark_confidence": benchmark_confidence,
        "diagnostic_constraints": diagnostic_constraints,
        "route_reasons": sorted(dict.fromkeys(route_reasons)),
        "reason_codes": sorted(dict.fromkeys(reason_codes)),
        "accuracy_gates": accuracy_gates,
        "workload_features": features,
    }


def _cpu_route_strategy(features: Mapping[str, Any], reason_codes: Sequence[str]) -> tuple[str, list[str]]:
    metric = str(features.get("metric", "levenshtein"))
    reasons = [str(code) for code in reason_codes]
    fixed_length = bool(features.get("uniform_target_length")) and features.get("target_length") is not None
    acgt_packable = bool(features.get("acgt_packable"))

    if metric == "levenshtein":
        return "cpu_levenshtein_indexed", reasons + ["levenshtein_indexed_cpu"]

    if metric == "hamming":
        hamming_reasons = ["hamming_distance_cpu"]
        if fixed_length and acgt_packable:
            hamming_reasons.extend(["hamming_seed_index_available", "fixed_length_acgt_targets"])
            return "cpu_hamming_seed_index", reasons + hamming_reasons
        if not acgt_packable:
            hamming_reasons.append("non_acgt_targets_cpu_only")
        if not fixed_length:
            hamming_reasons.append("variable_length_targets_cpu_only")
        return "cpu_hamming_indexed", reasons + hamming_reasons

    return "cpu_generic_assignment", reasons + ["generic_cpu_assignment"]


def _benchmark_confidence(
    features: Mapping[str, Any],
    reason_codes: Sequence[str],
    compute_compatible: bool,
    public_gate: bool,
    allow_gpu: bool,
) -> str:
    if not allow_gpu:
        return "gpu_disabled"
    if not compute_compatible:
        return "unsupported_route"
    if public_gate and bool(features.get("public_gpu_evidence_validated")):
        return "public_prior"
    if reason_codes:
        return "nearest_prior_with_constraints"
    return "nearest_prior"


def _thread_hint(assay: AssaySpec, features: Mapping[str, Any]) -> dict[str, Any]:
    configured = _table(assay.data, "run").get("threads")
    configured_threads = configured if isinstance(configured, int) and configured > 0 else None
    max_threads = configured_threads or 4
    target_count = int(features.get("target_count") or 0)
    target_length = int(features.get("target_length") or 0)
    reasons: list[str] = []

    if configured_threads is not None:
        reasons.append("configured_threads_cap")
    else:
        reasons.append("default_threads_cap")

    if target_count <= 1024 or target_length <= 8:
        recommended = 1
        reasons.append("small_target_set")
    elif target_count <= 10000:
        recommended = min(4, max_threads)
        reasons.append("moderate_target_set")
    else:
        recommended = min(8, max_threads)
        reasons.append("large_target_set")

    recommended = max(1, min(recommended, max_threads))
    return {
        "recommended_threads": recommended,
        "max_threads": max_threads,
        "reason_codes": sorted(dict.fromkeys(reasons)),
    }


def _backend_optimizer_features(assay: AssaySpec) -> tuple[dict[str, Any], list[str]]:
    assignment = _table(assay.data, "assignment")
    metric = str(assignment.get("metric", "levenshtein"))
    k = int(assignment.get("k", 1))
    reason_codes: list[str] = []
    target_key = ""
    extract: Mapping[str, Any] = {}

    if assay.mode == "count":
        target_key = "targets"
        extract = _table(assay.data, "extract")
    elif assay.mode == "demux":
        target_key = "barcodes"
        extract = _table(assay.data, "extract")
    else:
        reason_codes.append("pair_count_not_gpu_supported")

    length = extract.get("length")
    target_count = 0
    target_lengths: list[int] = []
    acgt_packable = False
    uniform_target_length = False
    if target_key:
        try:
            target_set = _read_target_sequences(_spec_path(assay, target_key))
        except AssaySpecError:
            target_set = TargetSet(sequences=[], lengths=[])
            reason_codes.append("target_table_unreadable")
        target_count = len(target_set.sequences)
        target_lengths = target_set.lengths
        uniform_target_length = len(target_lengths) == 1
        acgt_packable = bool(target_set.sequences) and all(set(seq) <= {"A", "C", "G", "T"} for seq in target_set.sequences)
    if metric != "hamming":
        reason_codes.append("metric_not_gpu_supported")
    if k != 1:
        reason_codes.append("edit_radius_not_gpu_supported")
    if not isinstance(length, int) or not 1 <= length <= 32:
        reason_codes.append("target_length_not_gpu_supported")
    if target_lengths and target_lengths != [length]:
        reason_codes.append("variable_target_length_not_gpu_supported")
    if not uniform_target_length and target_key:
        reason_codes.append("variable_target_length_not_gpu_supported")
    if not acgt_packable and target_key:
        reason_codes.append("target_alphabet_not_gpu_packable")

    features = {
        "mode": assay.mode,
        "assay_type": assay.assay_type,
        "backend_mode": str(assay.backend["mode"]),
        "allow_gpu": bool(assay.backend["allow_gpu"]),
        "metric": metric,
        "k": k,
        "target_count": target_count,
        "target_length": int(length) if isinstance(length, int) else None,
        "target_lengths": target_lengths,
        "acgt_packable": acgt_packable,
        "uniform_target_length": uniform_target_length,
        "public_gpu_evidence_validated": _gpu_public_evidence_validated(assay),
    }
    return features, reason_codes


def _gpu_speed_model(features: Mapping[str, Any]) -> dict[str, Any]:
    target_count = max(int(features.get("target_count") or 0), 1)
    target_length = int(features.get("target_length") or 0)
    mode = str(features.get("mode", ""))
    assay_type = str(features.get("assay_type", ""))
    if mode == "count" and assay_type == "crispr" and bool(features.get("public_gpu_evidence_validated")):
        nearest = GPU_BENCHMARK_PRIORS[0]
        return {
            "model": "nearest_neighbor_benchmark_priors_v1",
            "training_rows": len(GPU_BENCHMARK_PRIORS),
            "estimated_total_speedup": round(float(nearest["total_speedup"]), 2),
            "nearest_workload": str(nearest["workload"]),
            "nearest_total_speedup": float(nearest["total_speedup"]),
        }
    ranked = []
    for row in GPU_BENCHMARK_PRIORS:
        target_penalty = abs(math.log10(target_count) - math.log10(int(row["target_count"])))
        length_penalty = abs(target_length - int(row["target_length"])) / 10.0
        mode_penalty = 0.0 if mode == row["mode"] else 1.0
        assay_penalty = 0.0 if assay_type == row["assay_type"] else 0.4
        distance = target_penalty + length_penalty + mode_penalty + assay_penalty
        ranked.append((distance, row))
    ranked.sort(key=lambda item: item[0])
    nearest = ranked[0][1]
    estimated = float(nearest["total_speedup"])
    return {
        "model": "nearest_neighbor_benchmark_priors_v1",
        "training_rows": len(GPU_BENCHMARK_PRIORS),
        "estimated_total_speedup": round(estimated, 2),
        "nearest_workload": str(nearest["workload"]),
        "nearest_total_speedup": float(nearest["total_speedup"]),
    }


def _gpu_expected_speedup_band(estimated_total_speedup: float) -> str:
    if estimated_total_speedup < 1.25:
        return "1x"
    if estimated_total_speedup < 3.0:
        return "1.5-3x"
    if estimated_total_speedup < 8.0:
        return "3-8x"
    return "8-15x"


def _write_backend_optimization(plan: AssayPlan, optimization: Mapping[str, Any]) -> None:
    path = plan.artifacts["backend_optimization"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(optimization, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _format_backend_optimization(optimization: Mapping[str, Any]) -> str:
    reason_codes = ", ".join(str(code) for code in optimization.get("reason_codes", []) or [])
    accuracy_gates = ", ".join(str(gate) for gate in optimization.get("accuracy_gates", []) or [])
    return "\n".join(
        [
            "DotMatch backend optimizer",
            f"authority: {optimization.get('authority', '')}",
            f"selected_backend: {optimization.get('selected_backend', '')}",
            f"candidate_backend: {optimization.get('candidate_backend', '')}",
            f"recommendation: {optimization.get('recommendation', '')}",
            f"expected_speedup_band: {optimization.get('expected_speedup_band', '')}",
            f"reason_codes: {reason_codes}",
            f"accuracy_gates: {accuracy_gates}",
        ]
    ) + "\n"


def _gpu_eligible(assay: AssaySpec) -> bool:
    assignment = _table(assay.data, "assignment")
    if int(assignment.get("k", 1)) != 1 or str(assignment.get("metric", "levenshtein")) != "hamming":
        return False
    extract: Mapping[str, Any]
    if assay.mode == "count":
        extract = _table(assay.data, "extract")
        target_key = "targets"
    elif assay.mode == "demux":
        extract = _table(assay.data, "extract")
        target_key = "barcodes"
    elif assay.mode == "pair-count":
        return False
    else:
        return False
    length = extract.get("length")
    if not isinstance(length, int) or not 1 <= length <= 32:
        return False
    try:
        target_set = _read_target_sequences(_spec_path(assay, target_key))
    except AssaySpecError:
        return False
    return target_set.lengths == [length] and all(set(seq) <= {"A", "C", "G", "T"} for seq in target_set.sequences)


def _gpu_public_evidence_validated(assay: AssaySpec) -> bool:
    return assay.mode == "count" and assay.assay_type == "crispr"


def _evidence_boundary(assay: AssaySpec) -> dict[str, str]:
    evidence_id = ASSAY_EVIDENCE_IDS.get(assay.assay_type, "")
    try:
        data = json.loads(importlib_resources.files("dotmatch").joinpath("data", "assay-evidence.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, ModuleNotFoundError, AttributeError):
        path = Path(__file__).resolve().parents[2] / "docs" / "assay-evidence.json"
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return _missing_evidence_boundary(evidence_id)
    for row in data.get("assays", []):
        if isinstance(row, dict) and row.get("id") == evidence_id:
            return {
                "id": str(row.get("id", "")),
                "status": str(row.get("status", "")),
                "label": str(row.get("label", "")),
                "claim_boundary": str(row.get("claim_boundary", "")),
                "biological_unit": str(row.get("biological_unit", "")),
                "unsupported_claims": "; ".join(str(item) for item in row.get("unsupported_claims", []) if item),
                "minimum_public_evidence": "; ".join(str(item) for item in row.get("minimum_public_evidence", []) if item),
            }
    return _missing_evidence_boundary(evidence_id)


def _missing_evidence_boundary(evidence_id: str) -> dict[str, str]:
    return {
        "id": evidence_id,
        "status": "missing",
        "label": "",
        "claim_boundary": "",
        "biological_unit": "",
        "unsupported_claims": "",
        "minimum_public_evidence": "",
    }


_FINDING_SEVERITY_ORDER = {"blocked": 0, "error": 1, "warning": 2, "info": 3}
_RELIABILITY_STATUS_HINTS = {
    "passed": "ready for downstream use",
    "needs_review": "review findings before interpreting counts",
    "failed": "QC thresholds failed; inspect the report and assay_fixes.tsv",
    "blocked": "resolve blockers before production use",
}
_PREFLIGHT_STATUS_HINTS = {
    "passed": "ready to run counting",
    "needs_review": "review findings before counting",
    "failed": "preflight checks failed; review assay_fixes.tsv",
    "blocked": "resolve blockers before production use",
}


def _reliability_status_hint(status: str, stage: str) -> str:
    hints = _PREFLIGHT_STATUS_HINTS if stage == "preflight" else _RELIABILITY_STATUS_HINTS
    return hints.get(status, _RELIABILITY_STATUS_HINTS.get(status, ""))


def _read_reliability_summary(plan: AssayPlan) -> dict[str, Any] | None:
    summary_path = plan.artifacts["reliability_summary"]
    if not summary_path.exists():
        return None
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    return summary if isinstance(summary, dict) else None


def _primary_reliability_reason(summary: Mapping[str, Any]) -> str:
    findings = list(summary.get("findings", []) or [])
    stage = str(summary.get("stage", ""))
    if stage:
        staged = [finding for finding in findings if str(finding.get("stage", "")) in {"", stage}]
        if staged:
            findings = staged
    for severity in ("blocked", "error", "warning"):
        for finding in findings:
            if str(finding.get("severity", "")) == severity and finding.get("finding_id"):
                return str(finding["finding_id"])
    return str(summary.get("overall_status", ""))


def _top_actionable_findings(summary: Mapping[str, Any], *, limit: int = 3) -> list[dict[str, str]]:
    findings = [finding for finding in (summary.get("findings", []) or []) if isinstance(finding, dict)]
    stage = str(summary.get("stage", ""))
    if stage:
        staged = [finding for finding in findings if str(finding.get("stage", "")) in {"", stage}]
        if staged:
            findings = staged
    ranked = sorted(
        findings,
        key=lambda finding: (
            _FINDING_SEVERITY_ORDER.get(str(finding.get("severity", "info")), 9),
            str(finding.get("finding_id", "")),
            str(finding.get("sample_id", "")),
        ),
    )
    actionable: list[dict[str, str]] = []
    for finding in ranked:
        if str(finding.get("severity", "info")) not in {"blocked", "error", "warning"}:
            continue
        actionable.append(finding)
        if len(actionable) >= limit:
            break
    return actionable


def _spec_user_label(spec_path: Path) -> str:
    try:
        return spec_path.resolve().relative_to(Path.cwd().resolve()).as_posix()
    except (OSError, ValueError):
        return spec_path.name


def _user_facing_path(plan: AssayPlan, path: str | Path) -> str:
    value = Path(path)
    anchor = plan.spec.path.parent.resolve()
    try:
        return value.resolve(strict=False).relative_to(anchor).as_posix()
    except (OSError, ValueError):
        return value.name or str(path)


def _reliability_next_step(plan: AssayPlan, summary: Mapping[str, Any]) -> str:
    status = str(summary.get("overall_status", ""))
    stage = str(summary.get("stage", ""))
    report = _user_facing_path(plan, plan.artifacts["reliability_report"])
    fixes = _user_facing_path(plan, plan.artifacts["assay_fixes"])
    if status == "passed":
        if stage == "preflight":
            return "Preflight passed; run dotmatch assay start or ./run.sh for check-and-run."
        return "Run passed reliability thresholds."
    if status == "needs_review":
        return f"Review {report} before downstream analysis."
    if status == "blocked" and str(summary.get("spec_status", "")) == "draft":
        return (
            'Promote status to "ready" in assay.toml after reviewing inference_report.json and '
            f"{fixes}."
        )
    if plan.artifacts["assay_fixes"].exists():
        return f"Apply edits from {fixes}, then rerun dotmatch assay start or ./run.sh."
    return f"Open {report} for findings and recommended actions."


def _command_assay_check(plan: AssayPlan) -> int:
    status = _preflight_for_assay_start(plan)
    return _finish_assay_preflight(plan, _spec_user_label(plan.spec.path), status, emit_ok_stdout=True)


def _finish_assay_preflight(
    plan: AssayPlan,
    spec_label: str | Path,
    status: str,
    *,
    emit_ok_stdout: bool = False,
    verb: str = "check",
) -> int:
    summary = _read_reliability_summary(plan) or {}
    reason = _primary_reliability_reason(summary)
    if status == "passed":
        print(f"{spec_label}: {verb} passed", file=sys.stderr)
        _print_reliability_verdict(plan)
        if emit_ok_stdout:
            print(f"{spec_label}: ok")
        return 0
    if status == "needs_review":
        print(f"{spec_label}: {verb} needs review ({reason})", file=sys.stderr)
        _print_reliability_verdict(plan)
        if emit_ok_stdout:
            print(f"{spec_label}: ok")
        return _reliability_exit_code(status)
    if status == "failed":
        print(f"{spec_label}: {verb} failed ({reason})", file=sys.stderr)
        _print_reliability_verdict(plan)
        return _reliability_exit_code(status)
    label = f"{verb} blocked ({reason})" if reason else f"{verb} blocked"
    print(f"{spec_label}: {label}", file=sys.stderr)
    _print_reliability_verdict(plan)
    return _reliability_exit_code(status)


def _read_reliability_status(plan: AssayPlan) -> str:
    summary_path = plan.artifacts["reliability_summary"]
    if not summary_path.exists():
        return ""
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return ""
    return str(summary.get("overall_status", ""))


def _reliability_exit_code(status: str) -> int:
    if status == "passed":
        return 0
    if status == "needs_review":
        return 1
    return 2


def _print_reliability_verdict(plan: AssayPlan, *, label: str = "", include_next: bool = True) -> None:
    summary = _read_reliability_summary(plan)
    if summary is None:
        return
    status = str(summary.get("overall_status", ""))
    report_path = Path(str(summary.get("artifacts", {}).get("reliability_report", plan.artifacts["reliability_report"])))
    fixes_path = Path(str(summary.get("artifacts", {}).get("assay_fixes", plan.artifacts["assay_fixes"])))
    report = _user_facing_path(plan, report_path)
    fixes = _user_facing_path(plan, fixes_path)
    stage = str(summary.get("stage", ""))
    display_stage = label or stage
    prefix = f"reliability ({display_stage})" if display_stage else "reliability"
    hint = _reliability_status_hint(status, display_stage if display_stage in {"preflight", "postrun"} else stage)
    print(f"{prefix}: {status}" + (f" — {hint}" if hint else ""), file=sys.stderr)
    for finding in _top_actionable_findings(summary):
        finding_id = str(finding.get("finding_id", "finding"))
        sample_id = str(finding.get("sample_id", ""))
        finding_label = f"{finding_id} ({sample_id})" if sample_id else finding_id
        detail = str(finding.get("recommended_action") or finding.get("message") or "").strip()
        if detail:
            print(f"finding: {finding_label}: {detail}", file=sys.stderr)
    print(f"report: {report}", file=sys.stderr)
    if summary.get("assay_fixes") or plan.artifacts["assay_fixes"].exists():
        print(f"fixes: {fixes}", file=sys.stderr)
    if include_next:
        print(f"next: {_reliability_next_step(plan, summary)}", file=sys.stderr)


def _assay_fix(
    fix_id: str,
    finding_id: str,
    section: str,
    key: str,
    current_value: str,
    suggested_value: str,
    rationale: str,
) -> dict[str, str]:
    return {
        "fix_id": fix_id,
        "finding_id": finding_id,
        "section": section,
        "key": key,
        "current_value": current_value,
        "suggested_value": suggested_value,
        "rationale": rationale,
    }


def _load_inference_report(plan: AssayPlan) -> dict[str, Any] | None:
    report_key = plan.spec.data.get("inference_report")
    if not report_key:
        return None
    report_path = _path_from_spec(plan.spec.path, str(report_key), allow_absolute=False, name="inference_report")
    if not report_path.exists():
        return None
    try:
        data = json.loads(report_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    return data if isinstance(data, dict) else None


def _load_inference_candidates(plan: AssayPlan) -> list[dict[str, Any]]:
    report = _load_inference_report(plan)
    if report is None:
        return []
    candidates = report.get("candidates", [])
    if not isinstance(candidates, list):
        return []
    return [item for item in candidates if isinstance(item, dict)]


def _best_inference_candidate(plan: AssayPlan) -> dict[str, Any] | None:
    report = _load_inference_report(plan)
    if report is None:
        return None
    chosen = report.get("chosen")
    if isinstance(chosen, dict):
        return chosen
    candidates = _load_inference_candidates(plan)
    return candidates[0] if candidates else None


def _extract_table(plan: AssayPlan) -> Mapping[str, Any]:
    if plan.spec.mode == "count":
        return _table(plan.spec.data, "extract")
    if plan.spec.mode == "demux":
        return _table(plan.spec.data, "extract")
    return {}


def _build_assay_fixes(
    plan: AssayPlan,
    findings: Sequence[Mapping[str, str]],
    manifest: Mapping[str, Any] | None,
) -> list[dict[str, str]]:
    fixes: list[dict[str, str]] = []
    seen: set[str] = set()
    finding_ids = {str(finding.get("finding_id", "")) for finding in findings}

    def add(
        fix_id: str,
        finding_id: str,
        section: str,
        key: str,
        current_value: Any,
        suggested_value: Any,
        rationale: str,
    ) -> None:
        if fix_id in seen:
            return
        seen.add(fix_id)
        fixes.append(
            _assay_fix(
                fix_id,
                finding_id,
                section,
                key,
                str(current_value),
                str(suggested_value),
                rationale,
            )
        )

    if "draft_assayspec" in finding_ids:
        add(
            "promote_status_ready",
            "draft_assayspec",
            "status",
            "status",
            plan.spec.status,
            "ready",
            "Promote the inferred AssaySpec after reviewing inference_report.json and inference_candidates.tsv.",
        )

    assignment = _table(plan.spec.data, "assignment")
    if "unsafe_targets" in finding_ids:
        add(
            "assignment_exact_matching",
            "unsafe_targets",
            "assignment",
            "k",
            assignment.get("k", 1),
            0,
            "The target library is unsafe at the configured correction radius; count exactly or redesign colliding targets.",
        )

    extract = _extract_table(plan)
    candidate = _best_inference_candidate(plan)
    extract_findings = finding_ids & {
        "autopsy_wrong_offset",
        "autopsy_wrong_length",
        "assignment_rate_below_min",
        "unmatched_rate_above_max",
        "invalid_rate_above_max",
    }
    if candidate and extract and extract_findings:
        current_start = extract.get("start")
        current_length = extract.get("length")
        candidate_start = candidate.get("start")
        candidate_length = candidate.get("length")
        if candidate_start is not None and str(candidate_start) != str(current_start):
            add(
                "extract_start_from_inference",
                next(iter(sorted(extract_findings))),
                "extract",
                "start",
                current_start,
                candidate_start,
                "Inference ranked a different extract start; update the fixed window before rerunning.",
            )
        if candidate_length is not None and str(candidate_length) != str(current_length):
            add(
                "extract_length_from_inference",
                next(iter(sorted(extract_findings))),
                "extract",
                "length",
                current_length,
                candidate_length,
                "Inference ranked a different extract length; update the fixed window before rerunning.",
            )

    if "ambiguous_rate_above_max" in finding_ids and int(assignment.get("k", 1)) > 0:
        add(
            "assignment_disable_correction",
            "ambiguous_rate_above_max",
            "assignment",
            "k",
            assignment.get("k", 1),
            0,
            "High ambiguous assignment usually means the correction radius is too permissive for this library.",
        )

    backend = plan.spec.data.get("backend")
    backend_table = backend if isinstance(backend, dict) else {}
    if "metal_validation_failed" in finding_ids:
        add(
            "backend_cpu_authority",
            "metal_validation_failed",
            "backend",
            "mode",
            backend_table.get("mode", "auto"),
            "cpu",
            "Experimental Metal counting did not match CPU authority; rerun on CPU before trusting counts.",
        )
    if "experimental_gpu_backend_without_validation" in finding_ids:
        add(
            "backend_cpu_until_validated",
            "experimental_gpu_backend_without_validation",
            "backend",
            "mode",
            backend_table.get("mode", "auto"),
            "cpu",
            "Use CPU assignment authority until Metal validation passes or is explicitly enabled.",
        )

    reliability = plan.spec.reliability
    if "pairwise_sample_correlation_below_min" in finding_ids:
        threshold = _reliability_thresholds(plan.spec)["min_pairwise_sample_pearson"]
        add(
            "reliability_lower_pairwise_correlation_threshold",
            "pairwise_sample_correlation_below_min",
            "reliability",
            "min_pairwise_sample_pearson",
            threshold,
            max(0.0, round(threshold - 0.05, 2)),
            "Lower the replicate concordance gate only after reviewing sample_qc.tsv and crispr_qc.json.",
        )

    representation_findings = sorted(
        finding_ids
        & {
            "coverage_fraction_below_min",
            "zero_count_fraction_above_max",
            "gini_index_above_max",
            "top_1pct_fraction_above_max",
        }
    )
    for representation_finding in representation_findings:
        add(
            f"reliability_exploratory_{representation_finding}",
            representation_finding,
            "reliability",
            "profile",
            reliability.get("profile", "production"),
            "exploratory",
            "Switch to exploratory review when representation metrics fail but you still need a diagnostic rerun.",
        )

    library_warning_findings = sorted(
        finding_ids
        & {
            "crispr_qc_guide_collision",
            "crispr_qc_non_acgt_guides",
            "crispr_qc_library_collision_audit_radius",
        }
    )
    if library_warning_findings:
        target_key = "targets" if plan.spec.mode == "count" else "barcodes"
        for library_finding in library_warning_findings:
            add(
                f"review_guide_library_{library_finding}",
                library_finding,
                "targets",
                "library",
                str(_spec_path(plan.spec, target_key)),
                "deduplicated library",
                "Resolve guide collisions or non-ACGT sequences in the staged target table before production counting.",
            )

    if "evidence_boundary_not_supported" in finding_ids and reliability.get("require_public_evidence_boundary"):
        add(
            "reliability_relax_evidence_boundary",
            "evidence_boundary_not_supported",
            "reliability",
            "require_public_evidence_boundary",
            True,
            False,
            "Disable the public-evidence gate only for internal exploratory runs with bounded claims.",
        )

    if manifest and manifest.get("autopsy_triggered"):
        autopsy_findings = (manifest.get("autopsy_artifacts", {}) or {}).get("findings")
        if autopsy_findings and Path(str(autopsy_findings)).exists():
            with Path(str(autopsy_findings)).open("r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    code = str(row.get("finding", ""))
                    if code == "wrong_offset" and candidate and extract:
                        candidate_start = candidate.get("start")
                        if candidate_start is not None and str(candidate_start) != str(extract.get("start")):
                            add(
                                "autopsy_extract_start_from_inference",
                                "autopsy_wrong_offset",
                                "extract",
                                "start",
                                extract.get("start", ""),
                                candidate_start,
                                "Autopsy flagged a likely offset shift; apply the top inference candidate start.",
                            )
    return fixes


def _format_fix_action(fix: Mapping[str, str]) -> str:
    section = fix.get("section", "")
    key = fix.get("key", "")
    current = fix.get("current_value", "")
    suggested = fix.get("suggested_value", "")
    rationale = fix.get("rationale", "")
    if section == "status":
        return f'In assay.toml set status = "{suggested}" ({rationale})'
    return f"In assay.toml [{section}] set {key} = {suggested} (was {current}). {rationale}"


def _apply_recommended_actions(
    findings: Sequence[Mapping[str, str]],
    fixes: Sequence[Mapping[str, str]],
) -> list[dict[str, str]]:
    fixes_by_finding: dict[str, list[Mapping[str, str]]] = {}
    for fix in fixes:
        finding_id = str(fix.get("finding_id", ""))
        fixes_by_finding.setdefault(finding_id, []).append(fix)
    enriched: list[dict[str, str]] = []
    for finding in findings:
        updated = dict(finding)
        finding_id = str(updated.get("finding_id", ""))
        related = fixes_by_finding.get(finding_id, [])
        if related:
            updated["recommended_action"] = " ".join(_format_fix_action(fix) for fix in related)
        enriched.append(updated)
    return enriched


def _write_assay_fixes(path: Path, fixes: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=ASSAY_FIX_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for fix in fixes:
            writer.writerow({field: str(fix.get(field, "")) for field in ASSAY_FIX_COLUMNS})


def _html_assay_fixes_table(fixes: Sequence[Mapping[str, str]]) -> str:
    if not fixes:
        return "<p class=\"ok\">No assay.toml edits were suggested.</p>"
    rows = ["<table><tr>"]
    for column in ASSAY_FIX_COLUMNS:
        rows.append(f"<th>{html.escape(column)}</th>")
    rows.append("</tr>")
    for fix in fixes:
        rows.append("<tr>")
        for column in ASSAY_FIX_COLUMNS:
            rows.append(f"<td>{html.escape(str(fix.get(column, '')))}</td>")
        rows.append("</tr>")
    rows.append("</table>")
    return "".join(rows)


def _finding(
    finding_id: str,
    severity: str,
    stage: str,
    sample_id: str,
    metric: str,
    observed: str,
    threshold: str,
    message: str,
    recommended_action: str,
    source_artifact: str,
) -> dict[str, str]:
    return {
        "finding_id": finding_id,
        "severity": severity,
        "stage": stage,
        "sample_id": sample_id,
        "metric": metric,
        "observed": observed,
        "threshold": threshold,
        "message": message,
        "recommended_action": recommended_action,
        "source_artifact": source_artifact,
    }


def _finding_counts(findings: Sequence[Mapping[str, str]]) -> dict[str, int]:
    counts = {"blocked": 0, "error": 0, "warning": 0, "info": 0}
    for finding in findings:
        severity = str(finding.get("severity", "info"))
        counts[severity] = counts.get(severity, 0) + 1
    return counts


def _overall_reliability_status(counts: Mapping[str, int]) -> str:
    if int(counts.get("blocked", 0)) > 0:
        return "blocked"
    if int(counts.get("error", 0)) > 0:
        return "failed"
    if int(counts.get("warning", 0)) > 0:
        return "needs_review"
    return "passed"


def _write_reliability_artifacts(plan: AssayPlan, summary: Mapping[str, Any]) -> None:
    path = plan.artifacts["reliability_summary"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_reliability_findings(plan.artifacts["reliability_findings"], summary.get("findings", []) or [])
    _write_assay_fixes(plan.artifacts["assay_fixes"], summary.get("assay_fixes", []) or [])
    _write_reliability_manifest_summary(plan.artifacts["reliability_manifest_summary"], summary)
    _write_reliability_report(plan.artifacts["reliability_report"], summary)


def _write_reliability_findings(path: Path, findings: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=RELIABILITY_FINDING_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for finding in findings:
            writer.writerow({field: str(finding.get(field, "")) for field in RELIABILITY_FINDING_COLUMNS})


def _write_reliability_manifest_summary(path: Path, summary: Mapping[str, Any]) -> None:
    counts = summary.get("finding_counts", {}) or {}
    header = ["overall_status", "profile", "finding_count", "blocked_count", "error_count", "warning_count", "report"]
    row = [
        str(summary.get("overall_status", "")),
        str(summary.get("profile", "")),
        str(len(summary.get("findings", []) or [])),
        str(counts.get("blocked", 0)),
        str(counts.get("error", 0)),
        str(counts.get("warning", 0)),
        str(summary.get("artifacts", {}).get("reliability_report", "")),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerow(row)


def _write_reliability_report(path: Path, summary: Mapping[str, Any]) -> None:
    evidence = summary.get("evidence_boundary", {}) or {}
    backend = summary.get("backend", {}) or {}
    sections = [
        "<!doctype html>",
        "<html><head><meta charset=\"utf-8\"><title>DotMatch Reliability Report</title>",
        "<style>body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;margin:32px;color:#18212f}"
        "table{border-collapse:collapse;width:100%;margin:12px 0}th,td{border:1px solid #d8dee4;padding:7px 9px;text-align:left;vertical-align:top}"
        "th{background:#eef2f7}.ok{color:#1a7f37}.warn{color:#9a6700}.bad{color:#cf222e}code{background:#eef2f7;padding:2px 4px;border-radius:4px}</style>",
        "</head><body>",
        "<h1>DotMatch Reliability Report</h1>",
        f"<p>Status: <strong>{html.escape(str(summary.get('overall_status', '')))}</strong></p>",
        "<h2>Backend</h2>",
        _mapping_table(backend),
        "<h2>Evidence Boundary</h2>",
        _mapping_table(evidence),
        "<h2>Recommended Assay Fixes</h2>",
        _html_assay_fixes_table(summary.get("assay_fixes", []) or []),
        "<h2>Findings</h2>",
        _html_findings_table(summary.get("findings", []) or []),
        "</body></html>\n",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(sections), encoding="utf-8")


def _html_findings_table(findings: Sequence[Mapping[str, str]]) -> str:
    if not findings:
        return "<p class=\"ok\">No reliability findings were recorded.</p>"
    rows = ["<table><tr>"]
    for column in RELIABILITY_FINDING_COLUMNS:
        rows.append(f"<th>{html.escape(column)}</th>")
    rows.append("</tr>")
    for finding in findings:
        rows.append("<tr>")
        for column in RELIABILITY_FINDING_COLUMNS:
            rows.append(f"<td>{html.escape(str(finding.get(column, '')))}</td>")
        rows.append("</tr>")
    rows.append("</table>")
    return "".join(rows)


def _write_manifest_summary(plan: AssayPlan, manifest: Mapping[str, Any]) -> None:
    path = plan.artifacts["manifest_summary"]
    header = [
        "schema_version",
        "mode",
        "assay_type",
        "status",
        "native_version",
        "autopsy_triggered",
        "warning_count",
        "production_warning_count",
        "sample_count",
        "primary_report",
        "manifest",
        "methods",
        "citation_bib",
        "software_versions",
    ]
    row = [
        str(manifest.get("schema_version", "")),
        str(manifest.get("mode", "")),
        str(manifest.get("assay_type", "")),
        str(manifest.get("status", plan.spec.status)),
        str(manifest.get("native_version", "")),
        "true" if manifest.get("autopsy_triggered") else "false",
        str(len(manifest.get("warnings", []) or [])),
        str(len(manifest.get("production_warnings", []) or [])),
        str(_sample_count(plan.spec)),
        str(plan.artifacts["assay_report"]),
        str(plan.artifacts["manifest"]),
        str(plan.artifacts["methods"]),
        str(plan.artifacts["citation_bib"]),
        str(plan.artifacts["software_versions"]),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerow(row)


def _write_assay_report(plan: AssayPlan, manifest: Mapping[str, Any]) -> None:
    path = plan.artifacts["assay_report"]
    artifacts = manifest.get("artifacts", {})
    warnings = list(manifest.get("warnings", []) or [])
    production_warnings = list(manifest.get("production_warnings", []) or [])
    autopsy_artifacts = manifest.get("autopsy_artifacts", {}) or {}
    failed_commands = [cmd for cmd in manifest.get("commands", []) if cmd.get("exit_code") not in (0, None)]
    status = "Needs Review" if failed_commands or warnings or production_warnings else "Ready"
    status_class = "warn" if status == "Needs Review" else "ok"

    sections: list[str] = [
        "<!doctype html>",
        "<html><head><meta charset=\"utf-8\"><title>DotMatch Assay Report</title>",
        "<style>",
        "body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;margin:0;color:#18212f;background:#f7f9fb;line-height:1.45}",
        "main{max-width:1180px;margin:0 auto;padding:32px}",
        "h1{font-size:32px;margin:0 0 8px}h2{margin-top:28px;border-bottom:1px solid #d8dee4;padding-bottom:6px}",
        ".lede{color:#4b5563;margin:0 0 20px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px}",
        ".card{background:#fff;border:1px solid #d8dee4;border-radius:8px;padding:14px}.label{font-size:12px;color:#57606a;text-transform:uppercase;letter-spacing:.04em}",
        ".value{font-size:20px;font-weight:650;margin-top:4px}.ok{color:#1a7f37}.warn{color:#9a6700}.bad{color:#cf222e}",
        "table{border-collapse:collapse;width:100%;background:#fff;margin:12px 0}th,td{border:1px solid #d8dee4;padding:7px 9px;text-align:left;vertical-align:top}th{background:#eef2f7}",
        "code{background:#eef2f7;padding:2px 4px;border-radius:4px}a{color:#0969da}.empty{color:#6e7781}",
        "</style></head><body><main>",
        "<h1>DotMatch Assay Report</h1>",
        "<p class=\"lede\">Workflow-facing summary for a fixed-window known-target assay. Ambiguous reads are not silently counted.</p>",
        "<h2>Run Status</h2>",
        "<div class=\"grid\">",
        _metric_card("Status", status, status_class),
        _metric_card("Mode", str(manifest.get("mode", ""))),
        _metric_card("Assay", str(manifest.get("assay_type", ""))),
        _metric_card("Spec Status", str(manifest.get("status", plan.spec.status))),
        _metric_card("Samples", str(_sample_count(plan.spec))),
        _metric_card("Autopsy", "Triggered" if manifest.get("autopsy_triggered") else "Not triggered", "warn" if manifest.get("autopsy_triggered") else "ok"),
        "</div>",
        "<h2>Inputs</h2>",
        _samples_table(plan.spec),
        "<h2>Reliability</h2>",
        _reliability_html(plan, path.parent),
        "<h2>Sample QC</h2>",
        _sample_qc_table(plan.artifacts.get("sample_qc")),
        "<h2>Warnings</h2>",
        _warnings_html(warnings + production_warnings),
        "<h2>Library Audit</h2>",
        _audit_html(plan),
        "<h2>Autopsy</h2>",
        _autopsy_html(autopsy_artifacts, path.parent),
        "<h2>Methods and Citation</h2>",
        _citation_artifacts_html(plan, path.parent),
        "<h2>Artifacts</h2>",
        _mapping_table(artifacts, path.parent),
        "<h2>Native Commands</h2>",
        _commands_table(manifest.get("commands", []) or [], path.parent),
        "</main></body></html>\n",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(sections), encoding="utf-8")


def _metric_card(label: str, value: str, css_class: str = "") -> str:
    cls = f" {css_class}" if css_class else ""
    return f"<div class=\"card\"><div class=\"label\">{html.escape(label)}</div><div class=\"value{cls}\">{html.escape(value)}</div></div>"


def _sample_count(assay: AssaySpec) -> int:
    if assay.mode == "count":
        return len(_samples(assay.data))
    return 1


def _samples_table(assay: AssaySpec) -> str:
    rows = ["<table><tr><th>Sample</th><th>FASTQ</th></tr>"]
    if assay.mode == "count":
        for sample in _samples(assay.data):
            rows.append(
                "<tr><td>{}</td><td>{}</td></tr>".format(
                    html.escape(str(sample.get("id", ""))),
                    html.escape(Path(str(sample.get("fastq", ""))).name),
                )
            )
    else:
        reads_key = "reads"
        rows.append(
            "<tr><td>{}</td><td>{}</td></tr>".format(
                html.escape(assay.mode),
                html.escape(Path(str(assay.data.get(reads_key, ""))).name),
            )
        )
    rows.append("</table>")
    return "".join(rows)


def _sample_qc_table(path: Path | None) -> str:
    if path is None or not path.exists():
        return "<p class=\"empty\">No sample QC table was produced for this mode.</p>"
    return _tsv_preview_table(path, 12)


def _reliability_html(plan: AssayPlan, report_dir: Path) -> str:
    artifacts = {
        "reliability_summary": plan.artifacts.get("reliability_summary", ""),
        "reliability_findings": plan.artifacts.get("reliability_findings", ""),
        "reliability_report": plan.artifacts.get("reliability_report", ""),
        "reliability_manifest_summary": plan.artifacts.get("reliability_manifest_summary", ""),
    }
    parts = [_mapping_table({key: str(value) for key, value in artifacts.items() if value}, report_dir)]
    summary = plan.artifacts.get("reliability_manifest_summary")
    if summary is not None and summary.exists():
        parts.append(_tsv_preview_table(summary, 4))
    return "".join(parts)


def _citation_artifacts_html(plan: AssayPlan, report_dir: Path) -> str:
    artifacts = {
        "methods": plan.artifacts.get("methods", ""),
        "citation_bib": plan.artifacts.get("citation_bib", ""),
        "software_versions": plan.artifacts.get("software_versions", ""),
    }
    parts = [_mapping_table({key: str(value) for key, value in artifacts.items() if value}, report_dir)]
    methods = plan.artifacts.get("methods")
    if methods is not None and methods.exists():
        preview = "\n".join(methods.read_text(encoding="utf-8").splitlines()[:12])
        parts.append(f"<pre><code>{html.escape(preview)}</code></pre>")
    return "".join(parts)


def _audit_html(plan: AssayPlan) -> str:
    paths = []
    for key in ["audit", "left_audit", "right_audit"]:
        artifact = plan.artifacts.get(key)
        if artifact is not None:
            paths.append(Path(artifact) / "audit_summary.tsv")
    blocks = []
    for path in paths:
        if path.exists():
            blocks.append(f"<h3>{html.escape(path.parent.name)}</h3>{_tsv_preview_table(path, 40)}")
    return "".join(blocks) if blocks else "<p class=\"empty\">No audit summary was available.</p>"


def _autopsy_html(artifacts: Mapping[str, str], report_dir: Path) -> str:
    if not artifacts:
        return "<p class=\"empty\">Autopsy was not triggered for this run.</p>"
    findings = artifacts.get("findings")
    parts = [_mapping_table(artifacts, report_dir)]
    if findings and Path(findings).exists():
        parts.append(_tsv_preview_table(Path(findings), 40))
    return "".join(parts)


def _warnings_html(warnings: Sequence[str]) -> str:
    if not warnings:
        return "<p class=\"ok\">No AssaySpec production warnings were recorded.</p>"
    items = "".join(f"<li>{html.escape(str(warning))}</li>" for warning in warnings)
    return f"<ul class=\"warn\">{items}</ul>"


def _mapping_table(mapping: Mapping[str, Any], report_dir: Path | None = None) -> str:
    if not mapping:
        return "<p class=\"empty\">No artifacts recorded.</p>"
    rows = ["<table><tr><th>Name</th><th>Path</th></tr>"]
    for key in sorted(mapping):
        value = str(mapping[key])
        rows.append(f"<tr><td>{html.escape(str(key))}</td><td>{_artifact_link(value, report_dir)}</td></tr>")
    rows.append("</table>")
    return "".join(rows)


def _commands_table(commands: Sequence[Mapping[str, Any]], report_dir: Path) -> str:
    if not commands:
        return "<p class=\"empty\">No native commands were recorded.</p>"
    rows = ["<table><tr><th>Step</th><th>Exit</th><th>Command</th></tr>"]
    for command in commands:
        argv = " ".join(shlex.quote(_report_path_label(str(part), report_dir)) for part in command.get("argv", []))
        rows.append(
            "<tr><td>{}</td><td>{}</td><td><code>{}</code></td></tr>".format(
                html.escape(str(command.get("name", ""))),
                html.escape(str(command.get("exit_code", ""))),
                html.escape(argv),
            )
        )
    rows.append("</table>")
    return "".join(rows)


def _tsv_preview_table(path: Path, max_rows: int) -> str:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return "<p class=\"empty\">Preview unavailable.</p>"
    if not lines:
        return "<p class=\"empty\">File is empty.</p>"
    rows = ["<table>"]
    for row_index, line in enumerate(lines[:max_rows]):
        tag = "th" if row_index == 0 else "td"
        cells = "".join(f"<{tag}>{html.escape(_html_table_cell(cell))}</{tag}>" for cell in line.split("\t"))
        rows.append(f"<tr>{cells}</tr>")
    rows.append("</table>")
    return "".join(rows)


def _html_table_cell(value: str) -> str:
    path = Path(value)
    if path.is_absolute():
        return path.name
    return value


def _artifact_link(value: str, report_dir: Path | None = None) -> str:
    path = Path(value)
    if report_dir is not None:
        try:
            relative = path.resolve(strict=False).relative_to(report_dir.resolve(strict=False))
        except (OSError, ValueError):
            return f"<code>{html.escape(path.name or '[external path]')}</code>"
        href = quote(relative.as_posix())
        label = html.escape(relative.as_posix())
        return f"<a href=\"{href}\">{label}</a>"
    escaped = html.escape(path.name or value)
    return f"<code>{escaped}</code>"


def _report_path_label(value: str, report_dir: Path) -> str:
    path = Path(value)
    if not path.is_absolute():
        return value
    try:
        return path.resolve(strict=False).relative_to(report_dir.resolve(strict=False)).as_posix()
    except (OSError, ValueError):
        return path.name or "[external path]"


def _append_audit_warnings(
    plan: AssayPlan,
    step: PlanStep,
    manifest: dict[str, Any],
    *,
    warn: bool = True,
) -> None:
    out_dir = Path(step.argv[-1])
    summary = out_dir / "audit_summary.json"
    if not summary.exists():
        return
    try:
        data = json.loads(summary.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return
    key = f"safe_at_k{plan.spec.k}"
    if data.get(key) is False:
        warning = (
            f"{step.name}: target audit reports {key}=false; "
            "review library safety in the reliability report before production use"
        )
        manifest["warnings"].append(warning)
        if warn:
            print(f"dotmatch assay: warning: {warning}", file=sys.stderr)


def _resolve_native(argv: Sequence[str], native: Path) -> list[str]:
    if not argv:
        return []
    if argv[0] == "dotmatch-native":
        return [str(native), *argv[1:]]
    if argv[0] == "dotmatch":
        return [sys.executable, "-m", "dotmatch.cli", *argv[1:]]
    return list(argv)


def _command_init(template: str, out: Path) -> int:
    text = _template_text(template)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    print(str(out))
    return 0


def _template_text(template: str) -> str:
    if template == "pair-count":
        return """schema_version = 1
mode = "pair-count"
assay_type = "generic"
left_targets = "left_targets.tsv"
right_targets = "right_targets.tsv"
reads = "reads.fastq.gz"

[run]
out_dir = "dotmatch_pair_out"

[left]
start = 0
length = 20

[right]
start = 24
length = 20

[assignment]
k = 1
metric = "hamming"

[outputs]
assignments = true
"""
    if template == "inline-barcode-demux":
        return """schema_version = 1
mode = "demux"
assay_type = "inline_barcode"
barcodes = "barcodes.tsv"
reads = "reads.fastq.gz"

[run]
out_dir = "dotmatch_demux_out"

[extract]
start = 0
length = "auto"

[assignment]
k = 1
metric = "hamming"

[outputs]
assignments = true
ambiguous = true
unmatched = true
"""
    assay_type = {
        "crispr": "crispr",
        "feature-barcode": "feature_barcode",
        "inline-barcode-count": "inline_barcode",
        "amplicon-panel": "amplicon_panel",
        "oligo-adapter": "oligo_adapter",
    }[template]
    format_name = "mageck" if template == "crispr" else "dotmatch"
    length = 20 if template != "inline-barcode-count" else 8
    return f"""schema_version = 1
mode = "count"
assay_type = "{assay_type}"
targets = "targets.tsv"

[[samples]]
id = "sample_1"
fastq = "sample_1.fastq.gz"

[run]
out_dir = "dotmatch_assay_out"
threads = 1

[extract]
start = 0
length = {length}

[assignment]
k = 1
metric = "levenshtein"
ambiguous = "discard"

[outputs]
format = "{format_name}"
assignments = true
ambiguous = true
unmatched = true
"""


def _table(data: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    value = data.get(name, {})
    if not isinstance(value, dict):
        raise AssaySpecError(f"{name} must be a table")
    return value


def _samples(data: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    samples = data.get("samples")
    if not isinstance(samples, list):
        raise AssaySpecError("samples must be an array of tables")
    return samples


def _require_equal(value: Any, expected: Any, name: str) -> None:
    if value != expected:
        raise AssaySpecError(f"{name} must be {expected!r}")


def _require_enum(value: Any, choices: set[str], name: str) -> None:
    if value is None:
        raise AssaySpecError(f"{name} is required")
    if str(value) not in choices:
        allowed = ", ".join(sorted(choices))
        raise AssaySpecError(f"{name} must be one of: {allowed}")


def _require_int_range(value: Any, low: int, high: int, name: str) -> None:
    if not isinstance(value, int) or not low <= value <= high:
        raise AssaySpecError(f"{name} must be an integer from {low} to {high}")


def _require_float_range(value: Any, low: float, high: float, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not low <= float(value) <= high:
        raise AssaySpecError(f"{name} must be a number from {low} to {high}")


def _require_bool(value: Any, name: str) -> None:
    if not isinstance(value, bool):
        raise AssaySpecError(f"{name} must be a boolean")


def _require_extract(data: Mapping[str, Any], name: str, *, allow_auto: bool = False) -> None:
    table = _table(data, name)
    if "start" not in table:
        raise AssaySpecError(f"{name}.start is required")
    _require_int_range(table["start"], 0, 10**9, f"{name}.start")
    length = table.get("length")
    if allow_auto and length == "auto":
        return
    _require_int_range(length, 1, 10**9, f"{name}.length")


def _require_path(assay: AssaySpec, key: str) -> None:
    _require_existing_path(assay, assay.data.get(key), key)


def _require_existing_path(assay: AssaySpec, value: Any, name: str) -> None:
    if not isinstance(value, str) or not value:
        raise AssaySpecError(f"{name} is required")
    path = _path_from_spec(assay.path, value, allow_absolute=True, name=name)
    if not path.exists():
        raise AssaySpecError(f"{name} does not exist: {path}")


def _spec_path(assay: AssaySpec, key: str) -> Path:
    return _path_from_spec(assay.path, str(assay.data[key]), allow_absolute=True, name=key)


def _path_from_spec(spec_path: Path, value: str, *, allow_absolute: bool = False, name: str = "path") -> Path:
    path = Path(value)
    if path.is_absolute():
        if not allow_absolute:
            raise AssaySpecError(f"{name} must be relative to the AssaySpec")
        return path
    base = spec_path.parent.resolve()
    resolved = (base / path).resolve(strict=False)
    try:
        resolved.relative_to(base)
    except ValueError as exc:
        raise AssaySpecError(f"{name} must stay inside the AssaySpec directory: {value}") from exc
    return resolved


def _require_safe_identifier(value: Any, name: str) -> None:
    text = str(value)
    if not SAMPLE_ID_RE.fullmatch(text):
        raise AssaySpecError(f"{name} must contain only letters, digits, '.', '_' or '-'")
