from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import re
import shlex
import subprocess
import sys
from importlib import resources as importlib_resources
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence
from urllib.parse import quote

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
    "require_public_evidence_boundary": True,
}
BACKEND_DEFAULTS: dict[str, Any] = {
    "mode": "auto",
    "allow_gpu": True,
}
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
    for key in ["min_assignment_rate", "max_ambiguous_rate", "max_unmatched_rate", "max_invalid_rate"]:
        if key in reliability:
            _require_float_range(reliability[key], 0.0, 1.0, f"reliability.{key}")

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
        ]
    )
    return "\n".join(lines) + "\n"


def run_assay_plan(plan: AssayPlan) -> int:
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
        "autopsy_thresholds": AUTOPSY_THRESHOLDS,
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
            _append_audit_warnings(plan, step, manifest)
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

    autopsy = sub.add_parser("autopsy", help="diagnose suspicious fixed-window assay runs")
    autopsy.add_argument("spec")
    autopsy.add_argument("--out-dir", required=True)

    for name in ["check", "plan", "run"]:
        child = sub.add_parser(name)
        child.add_argument("spec")

    args = parser.parse_args(list(argv))
    try:
        if args.command == "init":
            return _command_init(args.template, Path(args.out))
        if args.command == "infer":
            return _command_infer(args)
        assay = load_assay_spec(args.spec)
        if args.command == "autopsy":
            run_autopsy(assay, Path(args.out_dir))
            return 0
        if args.command == "check":
            plan = compile_assay_plan(assay)
            _write_preflight_reliability(plan)
            print(f"{assay.path}: ok")
            return 0
        if args.command == "run" and assay.status == "draft" and _blocks_on_draft_inference(assay):
            plan = compile_assay_plan(assay)
            reliability = _build_reliability_summary(plan, stage="preflight")
            _write_reliability_artifacts(plan, reliability)
            raise AssaySpecError("refusing to run draft AssaySpec; review inference report and promote status to 'ready'")
        plan = compile_assay_plan(assay)
        if args.command == "plan":
            print(format_plan(plan), end="")
            return 0
        return run_assay_plan(plan)
    except AssaySpecError as exc:
        print(f"dotmatch assay: {exc}", file=sys.stderr)
        return 2
    except FileNotFoundError as exc:
        print(f"dotmatch assay: {exc}", file=sys.stderr)
        return 2


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
    subprocess.run(_resolve_native(_audit_cmd(_spec_path(assay, "targets"), audit_dir, assay.k), native), check=False)
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
    subprocess.run(_resolve_native(_audit_cmd(_spec_path(assay, "barcodes"), audit_dir, assay.k), native), check=False)
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
        subprocess.run(_resolve_native(_audit_cmd(_spec_path(assay, target_key), audit_dir, assay.k), native), check=False)
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
    with sample_qc.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = row.get("sample_id", "")
            try:
                total = _required_float(row, "total_reads")
                invalid = _required_float(row, "invalid_reads")
                assignment_rate = _required_float(row, "assignment_rate")
                ambiguous_rate = _required_float(row, "ambiguous_rate")
                no_match_rate = _required_float(row, "no_match_rate")
            except AssaySpecError:
                continue
            invalid_rate = invalid / total if total else 0.0
            checks = [
                ("assignment_rate", assignment_rate, "<", thresholds["min_assignment_rate"]),
                ("ambiguous_rate", ambiguous_rate, ">", thresholds["max_ambiguous_rate"]),
                ("no_match_rate", no_match_rate, ">", thresholds["max_unmatched_rate"]),
                ("invalid_rate", invalid_rate, ">", thresholds["max_invalid_rate"]),
            ]
            for metric, value, op, threshold in checks:
                if (op == "<" and value < threshold) or (op == ">" and value > threshold):
                    reasons.append(f"{sample}: {metric} {op} {threshold} ({value:.6f})")
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


def _write_manifest(plan: AssayPlan, manifest: Mapping[str, Any]) -> None:
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


def _build_reliability_summary(plan: AssayPlan, *, stage: str, manifest: Mapping[str, Any] | None = None) -> dict[str, Any]:
    findings: list[dict[str, str]] = []
    findings.extend(_base_reliability_findings(plan, stage))
    findings.extend(_audit_reliability_findings(plan, stage))
    if manifest is not None:
        findings.extend(_command_reliability_findings(manifest, stage))
        findings.extend(_sample_qc_reliability_findings(plan, stage))
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
    summary = {
        "schema_version": 1,
        "stage": stage,
        "overall_status": _overall_reliability_status(counts),
        "mode": plan.spec.mode,
        "assay_type": plan.spec.assay_type,
        "spec_status": plan.spec.status,
        "profile": str(plan.spec.reliability["profile"]),
        "thresholds": _reliability_thresholds(plan.spec),
        "backend": _backend_summary(plan.spec),
        "evidence_boundary": _evidence_boundary(plan.spec),
        "findings": findings,
        "finding_counts": counts,
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
                total = _required_float(row, "total_reads")
                invalid = _required_float(row, "invalid_reads")
                assignment_rate = _required_float(row, "assignment_rate")
                ambiguous_rate = _required_float(row, "ambiguous_rate")
                no_match_rate = _required_float(row, "no_match_rate")
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
            invalid_rate = invalid / total if total else 0.0
            checks = [
                ("assignment_rate_below_min", "assignment_rate", assignment_rate, "<", thresholds["min_assignment_rate"]),
                ("ambiguous_rate_above_max", "ambiguous_rate", ambiguous_rate, ">", thresholds["max_ambiguous_rate"]),
                ("unmatched_rate_above_max", "no_match_rate", no_match_rate, ">", thresholds["max_unmatched_rate"]),
                ("invalid_rate_above_max", "invalid_rate", invalid_rate, ">", thresholds["max_invalid_rate"]),
            ]
            for finding_id, metric, observed, op, threshold in checks:
                failed = observed < threshold if op == "<" else observed > threshold
                if not failed:
                    continue
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
    if plan.spec.status == "draft" and reliability["fail_on_draft_inference"]:
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
    }


def _backend_summary(assay: AssaySpec) -> dict[str, str]:
    backend = assay.backend
    mode = str(backend["mode"])
    if mode == "cpu":
        gpu_status = "disabled_by_mode"
    elif not bool(backend["allow_gpu"]):
        gpu_status = "disabled_by_config"
    elif _gpu_eligible(assay):
        gpu_status = "eligible_but_not_used"
    else:
        gpu_status = "not_eligible"
    return {
        "mode": mode,
        "authority": "cpu",
        "selected": "cpu",
        "gpu_status": gpu_status,
    }


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


def _evidence_boundary(assay: AssaySpec) -> dict[str, str]:
    evidence_id = ASSAY_EVIDENCE_IDS.get(assay.assay_type, "")
    try:
        data = json.loads(importlib_resources.files("dotmatch").joinpath("data", "assay-evidence.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, ModuleNotFoundError, AttributeError):
        path = Path(__file__).resolve().parents[2] / "docs" / "assay-evidence.json"
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return {"id": evidence_id, "status": "missing", "label": "", "claim_boundary": ""}
    for row in data.get("assays", []):
        if isinstance(row, dict) and row.get("id") == evidence_id:
            return {
                "id": str(row.get("id", "")),
                "status": str(row.get("status", "")),
                "label": str(row.get("label", "")),
                "claim_boundary": str(row.get("claim_boundary", "")),
            }
    return {"id": evidence_id, "status": "missing", "label": "", "claim_boundary": ""}


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


def _append_audit_warnings(plan: AssayPlan, step: PlanStep, manifest: dict[str, Any]) -> None:
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
        warning = f"{step.name}: target audit reports {key}=false; continuing with explicit ambiguity handling"
        manifest["warnings"].append(warning)
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
