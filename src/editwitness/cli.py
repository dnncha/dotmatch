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

from ._version import EXACT_MODEL_VERSION, MODEL_VERSION, SCHEMA_VERSION, __version__
from .compare import compare_models
from .generate import expand_deletions
from .engine import analyze
from .fasta import init_from_fasta
from .io import (
    InputError, atomic_write, canonical_json, check_destinations, digest,
    load_manifest, verify_result,
)
from .models import Analysis, DeletionScan, Manifest, ScanResult
from .sequence import apply_edits
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
    for command in ("analyze", "validate", "scan", "witness", "expand-deletions", "compare-models"):
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
    p.add_argument("--legacy-model", action="store_true", help="Explicitly select historical original-site eligibility")
    _output_options(p)
    p = sub.add_parser("init", help="Create a starting manifest from local FASTA and exact primer sites")
    p.add_argument("--fasta", required=True)
    p.add_argument("--left-primer", required=True)
    p.add_argument("--right-primer", required=True)
    p.add_argument("--edit-position", required=True, type=int, help="Local 0-based substitution position")
    p.add_argument("--alternate", required=True)
    p.add_argument("--deletion-radius", type=int, help="Add a deletion grid extending this many reference bases around the edit")
    p.add_argument("--deletion-step", type=int, default=50, help="Grid step in bases (default 50; not a biological prior)")
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
        "distinct_alternatives": result.distinct_alternatives,
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
    schema["$id"] = f"urn:editwitness:schema:{kind}:{SCHEMA_VERSION}"
    return schema


def capabilities() -> dict[str, Any]:
    return {
        "kind": "editwitness.capabilities", "schema_version": SCHEMA_VERSION, "version": __version__,
        "default_for_omitted_model": MODEL_VERSION, "models": [MODEL_VERSION, EXACT_MODEL_VERSION], "network_required": False, "telemetry": False,
        "commands": ["demo", "init", "validate", "analyze", "witness", "scan", "schema", "verify", "doctor", "capabilities", "expand-deletions", "compare-models"],
        "input": "Strict JSON manifest; local 0-based half-open reference coordinates.",
        "output": "JSON to stdout by default; structured JSON errors to stderr; no progress text on stdout.",
        "exit_codes": {"0": "completed (may demonstrate ambiguity)", "2": "invalid input or usage",
                       "3": "I/O failure", "4": "ambiguity with --fail-on-ambiguity", "5": "integrity or replay mismatch"},
        "limits": {"reference_bases": 20000, "alleles": 128, "hypotheses": 1000,
                   "existing_assays": 16, "candidate_assays": 24, "scan_grid_pairs": 500000,
                   "exact_planner_useful_candidates": 18, "exact_site_matches": 512,
                   "exact_products_per_allele_assay": 128, "exact_product_bases": 20000000,
                   "generation_endpoint_pairs": 5000},
        "supports": ["explicit allele replacements", "diploid clonal hypotheses",
                     "full-insert or post-trim paired-end sequence presence", "counterexample witnesses",
                     "candidate-panel selection", "streaming local-deletion geometry scan",
                     "both-orientation exact local primer rematching", "bounded deletion-hypothesis generation",
                     "observation-model sensitivity comparison"],
        "does_not_support": ["raw-read analysis", "probabilistic PCR", "allele-fraction inference",
                             "empirical sensitivity", "genome-wide primer specificity", "clinical interpretation",
                             "mosaic samples", "inversions or translocations", "copy-number assay simulation"],
        "command_contracts": {
            "analyze": {"input": "manifest path or -", "output": "editwitness.analysis", "compact": "editwitness.summary"},
            "expand-deletions": {"input": "manifest with deletion_scan and homozygous expectation", "output": "manifest JSON"},
            "compare-models": {"input": "manifest", "output": "editwitness.model_comparison"},
            "verify": {"input": "full result and optional --manifest", "output": "editwitness.integrity"},
        },
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
                "default_for_omitted_model": MODEL_VERSION, "models": [MODEL_VERSION, EXACT_MODEL_VERSION], "network_used": False,
                "note": "Environment information, not an empirical validation check."}, 0
    if command == "demo":
        filename = "paired_end.json" if args.paired_end else "demo.json"
        data = json.loads(files("editwitness").joinpath("data", filename).read_text(encoding="utf-8"))
        data["schema_version"] = SCHEMA_VERSION
        data["observation_model"] = MODEL_VERSION if args.legacy_model else EXACT_MODEL_VERSION
        return Manifest.model_validate(data).model_dump(mode="json"), 0
    if command == "init":
        initialized = init_from_fasta(args.fasta, args.left_primer, args.right_primer, args.edit_position, args.alternate)
        data = initialized.model_dump(mode="json")
        data["observation_model"] = EXACT_MODEL_VERSION
        if args.deletion_radius is not None:
            if args.deletion_radius < 1 or args.deletion_step < 1:
                raise InputError("deletion radius and step must be positive integers")
            pos, n = args.edit_position, len(initialized.reference.sequence)
            data["deletion_scan"] = DeletionScan(
                start_min=max(0, pos-args.deletion_radius), start_max=pos,
                end_min=pos+1, end_max=min(n, pos+1+args.deletion_radius), step=args.deletion_step,
            ).model_dump(mode="json")
        return Manifest.model_validate(data).model_dump(mode="json"), 0
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
    if command == "expand-deletions":
        return expand_deletions(manifest).model_dump(mode="json"), 0
    if command == "compare-models":
        return compare_models(manifest), 0
    if command == "validate":
        return {"kind": "editwitness.manifest_validation", "valid": True,
                "manifest_sha256": digest(manifest.model_dump(mode="json")),
                "caveat": "Valid input syntax and model invariants; not experimental validation."}, 0
    if command == "scan":
        return scan_deletions(manifest).model_dump(mode="json"), 0
    result = analyze(manifest)
    if command == "witness":
        assessment = next((h for h in result.hypotheses if h.hypothesis_id == args.hypothesis), None)
        representative = assessment.representative_hypothesis if assessment else args.hypothesis
        witness = next((w for w in result.witnesses if w.hypothesis_id == representative), None)
        if witness is None:
            raise InputError("hypothesis is not an equivalent alternative in this analysis")
        response: dict[str, Any] = {"kind": "editwitness.witness", "schema_version": SCHEMA_VERSION,
                                    "analysis_sha256": result.result_sha256,
                                    "requested_hypothesis": args.hypothesis,
                                    "witness": witness.model_dump(mode="json"),
                                    "assumptions": list(result.assumptions)}
        response["allele_edits"] = [a.model_dump(mode="json") for a in result.alleles
                                    if a.id in set(witness.expected_alleles + witness.alternative_alleles)]
        if args.include_sequences:
            allele_ids = set(witness.expected_alleles + witness.alternative_alleles)
            response["local_allele_sequences"] = {
                a.id: apply_edits(manifest.reference.sequence, a.edits)
                for a in manifest.alleles if a.id in allele_ids
            }
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
