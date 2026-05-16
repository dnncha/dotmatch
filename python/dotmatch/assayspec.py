from __future__ import annotations

import argparse
import csv
import json
import shlex
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

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
TEMPLATES = {
    "crispr",
    "feature-barcode",
    "inline-barcode-count",
    "inline-barcode-demux",
    "amplicon-panel",
    "oligo-adapter",
    "pair-count",
}


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
        return _path_from_spec(self.path, str(_table(self.data, "run").get("out_dir", "dotmatch_assay_out")))

    @property
    def k(self) -> int:
        return int(_table(self.data, "assignment").get("k", 1))

    @property
    def status(self) -> str:
        return str(self.data.get("status", "ready"))


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
    _require_enum(assignment.get("ambiguity_policy", "best"), AMBIGUITY_POLICIES, "assignment.ambiguity_policy")
    _require_enum(assignment.get("ambiguous", "discard"), AMBIGUOUS_OUTPUT, "assignment.ambiguous")
    if int(assignment.get("k", 1)) == 2 and assignment.get("metric", "levenshtein") == "hamming":
        raise AssaySpecError("assignment.k=2 is only valid with assignment.metric='levenshtein'")

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
        "normalized_spec": out_dir / "assay.normalized.json",
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
    return "\n".join(shlex.join(step.argv) for step in plan.steps) + "\n"


def run_assay_plan(plan: AssayPlan) -> int:
    native = find_native_cli()
    out_dir = plan.spec.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_generated_files(plan)
    _write_normalized_spec(plan)

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "mode": plan.spec.mode,
        "assay_type": plan.spec.assay_type,
        "spec_path": str(plan.spec.path),
        "native_cli": str(native),
        "commands": [],
        "artifacts": {key: str(value) for key, value in plan.artifacts.items()},
        "inference_report": str(plan.spec.data.get("inference_report", "")),
        "autopsy_triggered": False,
        "autopsy_thresholds": AUTOPSY_THRESHOLDS,
        "autopsy_artifacts": {},
        "production_warnings": [],
        "warnings": [],
    }
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
        if result.returncode != 0:
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
            print(f"{assay.path}: ok")
            return 0
        if args.command == "run" and assay.status == "draft":
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
    with path.open("r", encoding="utf-8") as fh:
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
    format_name = "mageck" if assay_type == "crispr" else "dotmatch"
    command = "crispr" if assay_type == "crispr" else assay_type
    _write_text_file(
        out,
        f"""schema_version = 1
status = "{status}"
mode = "count"
assay_type = "{assay_type}"
targets = "{targets}"

[[samples]]
id = "{sample_id}"
fastq = "{reads}"

[run]
out_dir = "{out.with_suffix('').name}_out"
threads = 1

[extract]
start = {chosen["start"]}
length = {chosen["length"]}

[assignment]
k = 1
metric = "hamming"
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
status = "{status}"
mode = "demux"
assay_type = "{assay_type}"
barcodes = "{barcodes}"
reads = "{reads}"

[run]
out_dir = "{out.with_suffix('').name}_out"

[extract]
start = {chosen["start"]}
length = {chosen["length"]}

[assignment]
k = 1
metric = "hamming"

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
status = "{status}"
mode = "pair-count"
assay_type = "{assay_type}"
left_targets = "{left_targets}"
right_targets = "{right_targets}"
reads = "{reads}"

[run]
out_dir = "{out.with_suffix('').name}_out"

[left]
start = {left["start"]}
length = {left["length"]}

[right]
start = {right["start"]}
length = {right["length"]}

[assignment]
k = 1
metric = "hamming"

[outputs]
assignments = true
""",
    )


def _write_text_file(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


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
            str(_path_from_spec(assay.path, str(sample["fastq"]))),
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
        fh.write("sample\tfinding\tseverity\tevidence\tartifact\n")
        for finding in findings:
            fh.write(
                f"{finding['sample']}\t{finding['finding']}\t{finding['severity']}\t"
                f"{finding['evidence']}\t{finding['artifact']}\n"
            )


def _autopsy_trigger_reasons(plan: AssayPlan) -> list[str]:
    if plan.spec.mode != "count":
        return []
    sample_qc = plan.artifacts.get("sample_qc")
    if sample_qc is None or not sample_qc.exists():
        return []
    reasons: list[str] = []
    with sample_qc.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = row.get("sample_id", "")
            total = _float(row.get("total_reads"))
            invalid = _float(row.get("invalid_reads"))
            invalid_rate = invalid / total if total else 0.0
            checks = [
                ("assignment_rate", _float(row.get("assignment_rate")), "<", AUTOPSY_THRESHOLDS["assignment_rate_min"]),
                ("ambiguous_rate", _float(row.get("ambiguous_rate")), ">", AUTOPSY_THRESHOLDS["ambiguous_rate_max"]),
                ("no_match_rate", _float(row.get("no_match_rate")), ">", AUTOPSY_THRESHOLDS["no_match_rate_max"]),
                ("invalid_rate", invalid_rate, ">", AUTOPSY_THRESHOLDS["invalid_rate_max"]),
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

    first_sample = _samples(data)[0]
    validate = [
        "dotmatch-native",
        "validate",
        "--targets",
        str(_spec_path(assay, "targets")),
        "--reads",
        str(_path_from_spec(assay.path, str(first_sample["fastq"]))),
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
    _add_assignment_options(cmd, assignment, include_ambiguity_policy=False)
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
    if outputs.get("assignments"):
        artifacts["pair_assignments"] = out_dir / "pair_assignments.tsv"
        cmd.extend(["--assignments", str(artifacts["pair_assignments"])])
    steps.append(PlanStep("run", cmd))


def _audit_cmd(targets: Path, out_dir: Path, k: int) -> list[str]:
    return ["dotmatch-native", "audit", "--targets", str(targets), "--k", str(k), "--audit-mode", "auto", "--out-dir", str(out_dir)]


def _add_assignment_options(cmd: list[str], assignment: Mapping[str, Any], *, include_ambiguity_policy: bool = True) -> None:
    if include_ambiguity_policy and "ambiguity_policy" in assignment:
        cmd.extend(["--ambiguity-policy", str(assignment["ambiguity_policy"])])
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
        fh.write("sample_id\tfastq\n")
        for sample in _samples(plan.spec.data):
            fh.write(f"{sample['id']}\t{_path_from_spec(plan.spec.path, str(sample['fastq']))}\n")


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
    path = _path_from_spec(assay.path, value)
    if not path.exists():
        raise AssaySpecError(f"{name} does not exist: {path}")


def _spec_path(assay: AssaySpec, key: str) -> Path:
    return _path_from_spec(assay.path, str(assay.data[key]))


def _path_from_spec(spec_path: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return spec_path.parent / path
