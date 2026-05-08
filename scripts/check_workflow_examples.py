#!/usr/bin/env python3
"""Audit local workflow-manager examples for DotMatch.

This gate checks that the in-repository Nextflow, nf-core-style, Snakemake,
Galaxy, and MultiQC examples are complete enough to serve as release evidence.
It does not claim upstream or external workflow adoption.
"""

from __future__ import annotations

import argparse
import json
import xml.etree.ElementTree as ET
from pathlib import Path


WORKFLOW_FIXTURES = [
    "README.md",
    "crispr_library.csv",
    "sample_a.fastq",
    "sample_b.fastq",
    "expected_counts.mageck.tsv",
    "expected_sample_qc.tsv",
]
GALAXY_TEST_DATA = [
    "crispr_library.csv",
    "sample_a.fastq",
    "sample_b.fastq",
    "expected_counts.mageck.tsv",
]


class WorkflowAudit:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def _read(path: Path, result: WorkflowAudit) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except Exception as exc:
        result.failures.append(f"{path.as_posix()} could not be read: {exc}")
        return ""


def _require(text: str, needle: str, message: str, result: WorkflowAudit) -> None:
    if needle not in text:
        result.failures.append(message)


def check_snakemake(root: Path, result: WorkflowAudit) -> None:
    config_path = root / "examples" / "workflows" / "snakemake" / "config.json"
    snakefile_path = root / "examples" / "workflows" / "snakemake" / "Snakefile"
    try:
        config = json.loads(_read(config_path, result))
    except json.JSONDecodeError as exc:
        result.failures.append(f"Snakemake config.json is invalid JSON: {exc}")
        return

    required_keys = ["library", "samples", "guide_start", "guide_length", "metric", "outdir"]
    for key in required_keys:
        if key not in config:
            result.failures.append(f"Snakemake config.json missing {key}")
    if config.get("metric") not in {"hamming", "levenshtein"}:
        result.failures.append("Snakemake config.json metric must be hamming or levenshtein")
    if not isinstance(config.get("samples"), dict) or not config.get("samples"):
        result.failures.append("Snakemake config.json must define at least one sample")

    snakefile = _read(snakefile_path, result)
    _require(snakefile, "rule dotmatch_crispr_count", "Snakemake Snakefile missing rule dotmatch_crispr_count", result)
    _require(snakefile, "dotmatch crispr-count", "Snakemake Snakefile must run dotmatch crispr-count", result)
    _require(snakefile, "--ambiguous discard", "Snakemake Snakefile must keep ambiguity policy explicit", result)
    _require(snakefile, "--sample-qc", "Snakemake Snakefile must emit sample_qc.tsv for MultiQC", result)
    _require(snakefile, "sample_qc", "Snakemake Snakefile must declare sample_qc output", result)

    if not any("Snakemake" in failure for failure in result.failures):
        result.passed.append("Snakemake CRISPR workflow example present")


def check_nextflow(root: Path, result: WorkflowAudit) -> None:
    config = _read(root / "examples" / "workflows" / "nextflow" / "nextflow.config", result)
    workflow = _read(root / "examples" / "workflows" / "nextflow" / "main.nf", result)

    for needle in [
        "library = 'examples/crispr_guides/data/yusa_library.csv'",
        "samples = 'examples/workflows/nextflow/samples.tsv'",
        "guide_start = 23",
        "guide_length = 19",
        "metric = 'levenshtein'",
        "outdir = 'examples/workflows/nextflow/output'",
    ]:
        _require(config, needle, f"Nextflow config missing {needle}", result)
    for needle, message in [
        ("nextflow.enable.dsl=2", "Nextflow workflow must enable DSL2"),
        ("process DOTMATCH_CRISPR_COUNT", "Nextflow workflow missing DOTMATCH_CRISPR_COUNT process"),
        ("dotmatch crispr-count", "Nextflow workflow must run dotmatch crispr-count"),
        ("--ambiguous discard", "Nextflow workflow must keep ambiguity policy explicit"),
        ("--sample-qc", "Nextflow workflow must emit sample_qc.tsv for MultiQC"),
        ("path \"sample_qc.tsv\", emit: sample_qc", "Nextflow workflow must declare sample_qc output"),
        ("publishDir params.outdir", "Nextflow workflow must publish outputs to params.outdir"),
    ]:
        _require(workflow, needle, message, result)

    if not any("Nextflow" in failure for failure in result.failures):
        result.passed.append("Nextflow CRISPR workflow example present")


def check_nfcore(root: Path, result: WorkflowAudit) -> None:
    base = root / "examples" / "workflows" / "nf-core"
    readme = _read(base / "README.md", result)
    module = _read(base / "modules" / "local" / "dotmatch" / "crispr_count" / "main.nf", result)
    meta = _read(base / "modules" / "local" / "dotmatch" / "crispr_count" / "meta.yml", result)
    test_path = base / "modules" / "local" / "dotmatch" / "crispr_count" / "tests" / "main.nf.test"
    nf_test = _read(test_path, result) if test_path.is_file() else ""

    _require(readme, "nf-core-style module candidate", "nf-core README must describe a module candidate", result)
    _require(readme, "not an upstream nf-core module", "nf-core README must avoid upstream adoption claims", result)
    if "external adoption" in readme and "not external adoption" not in readme:
        result.failures.append("nf-core README must not claim external adoption")

    for needle, message in [
        ("process DOTMATCH_CRISPR_COUNT", "nf-core module missing DOTMATCH_CRISPR_COUNT process"),
        ("tuple val(meta), path(reads), path(library)", "nf-core module missing expected input tuple"),
        ("dotmatch crispr-count", "nf-core module must run dotmatch crispr-count"),
        ("--ambiguous discard", "nf-core module must keep ambiguity policy explicit"),
        ("--sample-qc", "nf-core module must emit sample_qc.tsv for MultiQC"),
        ("emit: sample_qc", "nf-core module must declare sample_qc output"),
        ("versions.yml", "nf-core module must emit versions.yml"),
        ("dotmatch --version", "nf-core module must record dotmatch --version"),
        ("task.ext.args", "nf-core module must expose task.ext.args"),
    ]:
        _require(module, needle, message, result)
    for needle in [
        "name: dotmatch_crispr_count",
        "Count CRISPR guides with DotMatch",
        "- dotmatch",
        "- crispr",
        "counts:",
        "summary:",
        "sample_qc:",
        "versions:",
    ]:
        _require(meta, needle, f"nf-core module metadata missing {needle}", result)
    if not nf_test:
        result.failures.append("nf-core module must include an nf-test candidate at tests/main.nf.test")
    else:
        for needle, message in [
            ("nextflow_process", "nf-core nf-test candidate must define a nextflow_process"),
            ('script "../main.nf"', "nf-core nf-test candidate must reference ../main.nf"),
            ('process "DOTMATCH_CRISPR_COUNT"', "nf-core nf-test candidate must test DOTMATCH_CRISPR_COUNT"),
            ("Channel.of", "nf-core nf-test candidate must build input channels"),
            ("examples/workflows/fixtures/crispr_library.csv", "nf-core nf-test candidate must use shared workflow fixtures"),
            ("sample_qc", "nf-core nf-test candidate must assert sample_qc output"),
            ("versions.yml", "nf-core nf-test candidate must assert versions.yml output"),
        ]:
            _require(nf_test, needle, message, result)

    if not any("nf-core" in failure for failure in result.failures):
        result.passed.append("nf-core-style module candidate present without adoption claim")


def check_multiqc(root: Path, result: WorkflowAudit) -> None:
    config = _read(root / "examples" / "workflows" / "multiqc" / "multiqc_config.yaml", result)
    sample_qc_path = root / "examples" / "workflows" / "multiqc" / "data" / "sample_qc.tsv"
    sample_qc = _read(sample_qc_path, result)

    for needle in [
        "custom_data:",
        "dotmatch_sample_qc:",
        'plot_type: "table"',
        'fn: "*sample_qc.tsv"',
        "assignment_rate:",
        "ambiguous_rate:",
        "no_match_rate:",
    ]:
        _require(config, needle, f"MultiQC config missing {needle}", result)

    header = sample_qc.splitlines()[0].split("\t") if sample_qc.splitlines() else []
    required_columns = [
        "sample_id",
        "fastq",
        "total_reads",
        "valid_extracted_reads",
        "assigned_reads",
        "exact_reads",
        "assignment_rate",
        "ambiguous_rate",
        "no_match_rate",
        "candidates_verified",
    ]
    for column in required_columns:
        if column not in header:
            result.failures.append(f"MultiQC sample_qc.tsv missing {column}")

    if not any("MultiQC" in failure for failure in result.failures):
        result.passed.append("MultiQC custom-content example present")


def check_galaxy(root: Path, result: WorkflowAudit) -> None:
    readme = _read(root / "examples" / "workflows" / "galaxy" / "README.md", result)
    test_data = root / "examples" / "workflows" / "galaxy" / "test-data"
    wrapper_path = root / "examples" / "workflows" / "galaxy" / "dotmatch_crispr_count.xml"
    try:
        wrapper = ET.parse(wrapper_path).getroot()
    except Exception as exc:
        result.failures.append(f"Galaxy wrapper XML could not be parsed: {exc}")
        return

    if wrapper.tag != "tool" or wrapper.attrib.get("id") != "dotmatch_crispr_count":
        result.failures.append("Galaxy wrapper must be tool id dotmatch_crispr_count")
    command = wrapper.findtext("command") or ""
    _require(command, "dotmatch crispr-count", "Galaxy wrapper command must run dotmatch crispr-count", result)
    _require(command, "--ambiguous", "Galaxy wrapper command must expose --ambiguous", result)
    _require(command, "--summary", "Galaxy wrapper command must include --summary", result)
    _require(command, "--sample-qc", "Galaxy wrapper command must include --sample-qc", result)
    requirements = [node.text for node in wrapper.findall("./requirements/requirement")]
    if "dotmatch" not in requirements:
        result.failures.append("Galaxy wrapper requirements must include dotmatch")
    output_names = {node.attrib.get("name", "") for node in wrapper.findall("./outputs/data")}
    if not {"counts", "summary", "sample_qc"} <= output_names:
        result.failures.append("Galaxy wrapper outputs must include counts, summary, and sample_qc")
    test = wrapper.find("./tests/test")
    if test is None:
        result.failures.append("Galaxy wrapper must include a Planemo test with tiny CRISPR fixtures")
    else:
        params = {node.attrib.get("name", ""): node.attrib.get("value", "") for node in test.findall("param")}
        for name, value in [
            ("library", "crispr_library.csv"),
            ("sample1_fastq", "sample_a.fastq"),
            ("sample1_label", "sample_a"),
            ("sample2_fastq", "sample_b.fastq"),
            ("sample2_label", "sample_b"),
        ]:
            if params.get(name) != value:
                result.failures.append(f"Galaxy Planemo test must set {name}={value}")
        test_outputs = {node.attrib.get("name", ""): node for node in test.findall("output")}
        counts = test_outputs.get("counts")
        if counts is None or counts.attrib.get("file") != "expected_counts.mageck.tsv":
            result.failures.append("Galaxy Planemo test must compare counts to expected_counts.mageck.tsv")
        sample_qc = test_outputs.get("sample_qc")
        if sample_qc is None:
            result.failures.append("Galaxy Planemo test must assert sample_qc output")
        elif sample_qc.find("./assert_contents/has_text[@text='assignment_rate']") is None:
            result.failures.append("Galaxy Planemo test must assert sample_qc assignment_rate content")
    for filename in GALAXY_TEST_DATA:
        if not (test_data / filename).is_file():
            result.failures.append(f"Galaxy Planemo test-data file is missing: {filename}")
    _require(readme, "example wrapper", "Galaxy README must describe an example wrapper", result)
    _require(readme, "not a ToolShed release", "Galaxy README must avoid ToolShed release claims", result)
    _require(readme, "planemo", "Galaxy README must mention planemo linting", result)

    if not any("Galaxy" in failure for failure in result.failures):
        result.passed.append("Galaxy wrapper example present without ToolShed claim")


def check_submission_fixtures(root: Path, result: WorkflowAudit) -> None:
    fixtures = root / "examples" / "workflows" / "fixtures"
    for filename in WORKFLOW_FIXTURES:
        if not (fixtures / filename).is_file():
            result.failures.append(f"workflow test fixture is missing: examples/workflows/fixtures/{filename}")
    readme = _read(fixtures / "README.md", result) if (fixtures / "README.md").is_file() else ""
    for outcome in ["unique", "ambiguous", "unmatched", "invalid"]:
        _require(readme, outcome, f"workflow fixture README must describe {outcome} outcome", result)
    sample_qc_path = fixtures / "expected_sample_qc.tsv"
    sample_qc = _read(sample_qc_path, result) if sample_qc_path.is_file() else ""
    header = sample_qc.splitlines()[0].split("\t") if sample_qc.splitlines() else []
    for column in ["sample_id", "assignment_rate", "ambiguous_rate", "no_match_rate", "invalid_reads"]:
        if column not in header:
            result.failures.append(f"workflow expected_sample_qc.tsv missing {column}")
    counts_path = fixtures / "expected_counts.mageck.tsv"
    counts = _read(counts_path, result) if counts_path.is_file() else ""
    _require(counts, "sgRNA\tGene\tsample_a\tsample_b", "workflow expected_counts.mageck.tsv must be MAGeCK-compatible", result)

    if not any("workflow test fixture" in failure or "workflow fixture" in failure or "expected_sample_qc" in failure for failure in result.failures):
        result.passed.append("shared workflow submission fixtures present")


def check_adoption_submission_dossier(root: Path, result: WorkflowAudit) -> None:
    dossier = _read(root / "docs" / "workflow-adoption-submission.md", result)
    required = [
        ("Workflow Adoption Submission Dossier", "workflow adoption submission dossier missing title"),
        ("not external adoption", "workflow adoption submission dossier must avoid external adoption claims"),
        ("make workflow-examples-ready", "workflow adoption submission dossier must require workflow-examples-ready"),
        ("make workflow-adoption-status", "workflow adoption submission dossier must mention post-adoption status gate"),
        ("docs/workflow-adoption.json", "workflow adoption submission dossier must name workflow-adoption.json"),
        ("nf-core/modules", "workflow adoption submission dossier must include nf-core/modules target"),
        ("Galaxy ToolShed", "workflow adoption submission dossier must include Galaxy ToolShed target"),
        ("MultiQC", "workflow adoption submission dossier must include MultiQC target"),
        ("Snakemake", "workflow adoption submission dossier must include Snakemake target"),
        ("Nextflow", "workflow adoption submission dossier must include Nextflow target"),
        ("adoption_url", "workflow adoption submission dossier must include adoption_url manifest field"),
        ("evidence_url", "workflow adoption submission dossier must include evidence_url manifest field"),
        ("validation_notes", "workflow adoption submission dossier must include validation_notes manifest field"),
        ("recorded_date", "workflow adoption submission dossier must include recorded_date manifest field"),
        ("stable external HTTPS adoption and evidence URLs", "workflow adoption submission dossier must require stable external HTTPS URLs before ready status"),
    ]
    for needle, message in required:
        _require(dossier, needle, message, result)

    if not any("workflow adoption submission dossier" in failure for failure in result.failures):
        result.passed.append("workflow adoption submission dossier present")


def audit(root: Path) -> WorkflowAudit:
    root = root.resolve()
    result = WorkflowAudit()
    check_submission_fixtures(root, result)
    check_snakemake(root, result)
    check_nextflow(root, result)
    check_nfcore(root, result)
    check_multiqc(root, result)
    check_galaxy(root, result)
    check_adoption_submission_dossier(root, result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="repository root")
    args = parser.parse_args()

    result = audit(Path(args.root))
    for item in result.passed:
        print(f"PASS: {item}")
    for item in result.failures:
        print(f"FAIL: {item}")
    if result.ok:
        print("WORKFLOW EXAMPLES: PASS")
        return 0
    print("WORKFLOW EXAMPLES: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
