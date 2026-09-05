"""Agent-friendly JSON commands, bounded I/O and stable exit codes."""
from __future__ import annotations

import argparse
import json
import platform
import sys
from importlib.resources import files
from pathlib import Path
from typing import Any, NoReturn

from pydantic import BaseModel, ValidationError

from ._version import MODEL_VERSION, __version__
from .engine import analyze
from .fasta import init_from_fasta
from .io import (
    InputError, atomic_write, canonical_json, check_destinations, digest,
    load_manifest, verify_result,
)
from .models import Analysis, Manifest, ScanResult
from .report import render_report
from .scan import scan_deletions


class Parser(argparse.ArgumentParser):
    def error(self, message: str) -> NoReturn:
        raise InputError(message)


def _output_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-o", "--output", default="-", help="JSON destination; '-' is stdout")
    parser.add_argument("--pretty", action="store_true", help="Indent JSON for human reading")
    parser.add_argument("--force", action="store_true", help="Explicitly replace existing output files")


def build_parser() -> Parser:
    parser = Parser(prog="editwitness", description="Show what a CRISPR validation assay cannot establish.")
    parser.add_argument("--version", action="version", version=f"editwitness {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("analyze", "validate", "scan", "witness"):
        p = sub.add_parser(command)
        p.add_argument("manifest", help="Local manifest JSON; '-' reads stdin")
        _output_options(p)
        if command == "analyze":
            p.add_argument("--html", help="Also write a self-contained HTML report")
            p.add_argument("--compact", action="store_true", help="Return a token-efficient summary without sequences")
            p.add_argument("--fail-on-ambiguity", action="store_true", help="Exit 4 after emitting a result with counterexamples")
        if command == "witness":
            p.add_argument("--hypothesis", required=True)
            p.add_argument("--include-sequences", action="store_true")
    p = sub.add_parser("verify", help="Verify result integrity; optionally replay an analysis")
    p.add_argument("result")
    p.add_argument("--manifest", help="Also reproduce and compare the result from this manifest")
    _output_options(p)
    p = sub.add_parser("schema")
    p.add_argument("kind", choices=("manifest", "analysis", "scan"))
    _output_options(p)
    for command in ("doctor", "capabilities"):
        _output_options(sub.add_parser(command))
    p = sub.add_parser("demo", help="Write an explicitly synthetic example manifest")
    p.add_argument("--paired-end", action="store_true", help="Demonstrate a read-gap blind spot too")
    _output_options(p)
    p = sub.add_parser("init", help="Create a starting manifest from local FASTA and exact primer sites")
    p.add_argument("--fasta", required=True)
    p.add_argument("--left-primer", required=True)
    p.add_argument("--right-primer", required=True)
    p.add_argument("--edit-position", required=True, type=int, help="Local 0-based substitution position")
    p.add_argument("--alternate", required=True)
    _output_options(p)
    return parser


def compact_summary(result: Analysis) -> dict[str, Any]:
    return {
        "kind": "editwitness.summary", "schema_version": result.schema_version,
        "package_version": result.package_version, "model_version": result.model_version,
        "analysis_sha256": result.result_sha256, "manifest_sha256": result.manifest_sha256,
        "conclusion": result.conclusion, "validation_status": result.validation_status,
        "expected_hypothesis": result.expected_hypothesis,
        "declared_hypotheses": len(result.hypotheses),
        "equivalent_alternatives": [w.hypothesis_id for w in result.witnesses],
        "plan": result.plan.model_dump(mode="json"),
        "notice_codes": [n.code for n in result.notices],
        "caveat": "Finite declared hypotheses; idealized sequence-presence model. "
                  "Not a defect diagnosis or safety certificate. Use full output for evidence and assumptions.",
    }


def schema_for(kind: str) -> dict[str, Any]:
    models: dict[str, type[BaseModel]] = {
        "manifest": Manifest, "analysis": Analysis, "scan": ScanResult,
    }
    schema = models[kind].model_json_schema()
    schema["$schema"] = "https://json-schema.org/draft/2020-12/schema"
    schema["$id"] = f"urn:editwitness:schema:{kind}:1.0"
    return schema


def capabilities() -> dict[str, Any]:
    return {
        "kind": "editwitness.capabilities", "schema_version": "1.0", "version": __version__,
        "model_version": MODEL_VERSION, "network_required": False, "telemetry": False,
        "commands": ["demo", "init", "validate", "analyze", "witness", "scan", "schema", "verify", "doctor", "capabilities"],
        "input": "Strict JSON manifest; local 0-based half-open reference coordinates.",
        "output": "JSON to stdout by default; structured JSON errors to stderr; no progress text on stdout.",
        "exit_codes": {"0": "completed (may demonstrate ambiguity)", "2": "invalid input or usage",
                       "3": "I/O failure", "4": "ambiguity with --fail-on-ambiguity", "5": "integrity or replay mismatch"},
        "limits": {"reference_bases": 20000, "alleles": 128, "hypotheses": 1000,
                   "existing_assays": 16, "candidate_assays": 24, "scan_grid_pairs": 500000,
                   "exact_planner_useful_candidates": 18},
        "supports": ["explicit allele replacements", "diploid clonal hypotheses",
                     "full-insert or post-trim paired-end sequence presence", "counterexample witnesses",
                     "candidate-panel selection", "streaming local-deletion geometry scan"],
        "does_not_support": ["raw-read analysis", "probabilistic PCR", "allele-fraction inference",
                             "empirical sensitivity", "genome-wide primer specificity", "clinical interpretation",
                             "mosaic samples", "inversions or translocations", "copy-number assay simulation"],
        "agent_rule": "Never translate a successful process exit or absence of a declared counterexample into 'safe', 'validated' or 'biallelic confirmed'.",
    }


def run(args: argparse.Namespace) -> tuple[dict[str, Any], int]:
    command = args.command
    if command == "schema":
        return schema_for(args.kind), 0
    if command == "capabilities":
        return capabilities(), 0
    if command == "doctor":
        import pydantic
        return {"kind": "editwitness.doctor", "version": __version__,
                "python": platform.python_version(), "pydantic": pydantic.__version__,
                "model_version": MODEL_VERSION, "network_used": False,
                "note": "Environment information, not an empirical validation check."}, 0
    if command == "demo":
        filename = "paired_end.json" if args.paired_end else "demo.json"
        data = json.loads(files("editwitness").joinpath("data", filename).read_text(encoding="utf-8"))
        return Manifest.model_validate(data).model_dump(mode="json"), 0
    if command == "init":
        initialized = init_from_fasta(args.fasta, args.left_primer, args.right_primer, args.edit_position, args.alternate)
        return initialized.model_dump(mode="json"), 0
    if command == "verify":
        verified = verify_result(args.result)
        if args.manifest is not None:
            manifest = load_manifest(args.manifest)
            replay = analyze(manifest) if isinstance(verified, Analysis) else scan_deletions(manifest)
            if replay.result_sha256 != verified.result_sha256:
                raise InputError("result replay mismatch (inputs or model/package version may differ)")
        return {"kind": "editwitness.integrity", "checksum_matches": True,
                "replayed": args.manifest is not None, "result_sha256": verified.result_sha256,
                "caveat": "Content integrity only; this does not establish scientific validity or authentic authorship."}, 0
    manifest = load_manifest(args.manifest)
    if command == "validate":
        return {"kind": "editwitness.manifest_validation", "valid": True,
                "manifest_sha256": digest(manifest.model_dump(mode="json")),
                "caveat": "Valid input syntax and model invariants; not experimental validation."}, 0
    if command == "scan":
        return scan_deletions(manifest).model_dump(mode="json"), 0
    result = analyze(manifest)
    if command == "witness":
        witness = next((w for w in result.witnesses if w.hypothesis_id == args.hypothesis), None)
        if witness is None:
            raise InputError("hypothesis is not an equivalent alternative in this analysis")
        response: dict[str, Any] = {"kind": "editwitness.witness", "schema_version": "1.0",
                                    "analysis_sha256": result.result_sha256,
                                    "witness": witness.model_dump(mode="json"),
                                    "assumptions": list(result.assumptions)}
        if args.include_sequences:
            allele_ids = set(witness.expected_alleles + witness.alternative_alleles)
            response["allele_observations"] = [o.model_dump(mode="json") for o in result.allele_observations
                                               if o.allele_id in allele_ids]
        return response, 0
    if args.html:
        atomic_write(args.html, render_report(result), force=args.force)
    data = compact_summary(result) if args.compact else result.model_dump(mode="json")
    return data, 4 if args.fail_on_ambiguity and result.witnesses else 0


def main(argv: list[str] | None = None) -> int:
    args: argparse.Namespace | None = None
    try:
        args = build_parser().parse_args(argv)
        input_path = getattr(args, "manifest", None) or getattr(args, "result", None) or getattr(args, "fasta", None)
        if args.command == "verify" and args.output != "-" and Path(args.output).resolve() == Path(args.result).resolve():
            raise InputError("refusing to replace the result being verified")
        check_destinations([args.output, getattr(args, "html", None)], force=args.force, input_path=input_path)
        if getattr(args, "html", None) == "-":
            raise InputError("--html requires a file path, not stdout")
        data, status = run(args)
        text = (json.dumps(data, indent=2, ensure_ascii=True, allow_nan=False) if args.pretty else canonical_json(data)) + "\n"
        if args.output == "-":
            sys.stdout.write(text)
        else:
            atomic_write(args.output, text, force=args.force)
        return status
    except ValidationError as error:
        payload = {"kind": "editwitness.error", "code": "INVALID_MANIFEST_OR_RESULT",
                   "message": "Input violates the strict data contract.",
                   "details": error.errors(include_input=False, include_context=False, include_url=False)}
        status = 2
    except InputError as error:
        integrity = args is not None and args.command == "verify" and (
            "checksum mismatch" in str(error) or "replay mismatch" in str(error)
        )
        payload = {"kind": "editwitness.error", "code": "INTEGRITY_MISMATCH" if integrity else "INVALID_INPUT",
                   "message": str(error)}
        status = 5 if integrity else 2
    except (OSError, UnicodeError) as error:
        payload = {"kind": "editwitness.error", "code": "IO_ERROR", "message": str(error)}
        status = 3
    sys.stderr.write(canonical_json(payload) + "\n")
    return status
