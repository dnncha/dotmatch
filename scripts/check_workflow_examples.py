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
    "crispr_assay.toml",
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
    _require(snakefile, "rule dotmatch_assay_run", "Snakemake Snakefile missing rule dotmatch_assay_run", result)
    _require(snakefile, "dotmatch assay run", "Snakemake AssaySpec rule must run dotmatch assay run", result)
    _require(snakefile, "assay_report.html", "Snakemake AssaySpec rule must expose assay_report.html", result)
    _require(snakefile, "assay_manifest.json", "Snakemake AssaySpec rule must expose assay_manifest.json", result)
    _require(snakefile, "assay_manifest.summary.tsv", "Snakemake AssaySpec rule must expose assay_manifest.summary.tsv", result)
    _require(snakefile, "crispr_qc.html", "Snakemake AssaySpec rule must expose crispr_qc.html", result)
    _require(snakefile, "crispr_qc.json", "Snakemake AssaySpec rule must expose crispr_qc.json", result)
    _require(snakefile, "crispr_qc.summary.tsv", "Snakemake AssaySpec rule must expose crispr_qc.summary.tsv", result)
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
        ("process DOTMATCH_ASSAY_RUN", "Nextflow workflow missing DOTMATCH_ASSAY_RUN process"),
        ("dotmatch assay run", "Nextflow AssaySpec workflow must run dotmatch assay run"),
        ("path \"assay_report.html\", emit: assay_report", "Nextflow AssaySpec workflow must emit assay_report.html"),
        ("path \"assay_manifest.json\", emit: assay_manifest", "Nextflow AssaySpec workflow must emit assay_manifest.json"),
        ("path \"assay_manifest.summary.tsv\", emit: assay_manifest_summary", "Nextflow AssaySpec workflow must emit assay_manifest.summary.tsv"),
        ("path \"crispr_qc.html\", emit: assay_crispr_qc_report", "Nextflow AssaySpec workflow must emit crispr_qc.html"),
        ("path \"crispr_qc.json\", emit: assay_crispr_qc_json", "Nextflow AssaySpec workflow must emit crispr_qc.json"),
        ("path \"crispr_qc.summary.tsv\", emit: assay_crispr_qc_summary", "Nextflow AssaySpec workflow must emit crispr_qc.summary.tsv"),
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
    _require(readme, "not been submitted to or", "nf-core README must avoid upstream adoption claims", result)
    _require(readme, "dotmatch_assay_run", "nf-core README must describe the AssaySpec module candidate", result)
    if "external adoption" in readme:
        result.failures.append("nf-core README must not use external-adoption maintainer language")

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

    assay_base = base / "modules" / "local" / "dotmatch" / "assay_run"
    assay_module = _read(assay_base / "main.nf", result)
    assay_meta = _read(assay_base / "meta.yml", result)
    assay_test_path = assay_base / "tests" / "main.nf.test"
    assay_nf_test = _read(assay_test_path, result) if assay_test_path.is_file() else ""
    for needle, message in [
        ("process DOTMATCH_ASSAY_RUN", "nf-core AssaySpec module missing DOTMATCH_ASSAY_RUN process"),
        ("tuple val(meta), path(assay_spec), path(assay_inputs)", "nf-core AssaySpec module missing expected input tuple"),
        ("basename", "nf-core AssaySpec module must stage assay input files by basename"),
        ("dotmatch assay run", "nf-core AssaySpec module must run dotmatch assay run"),
        ("emit: assay_report", "nf-core AssaySpec module must emit assay_report"),
        ("emit: assay_manifest", "nf-core AssaySpec module must emit assay_manifest"),
        ("emit: assay_manifest_summary", "nf-core AssaySpec module must emit assay_manifest_summary"),
        ("emit: sample_qc", "nf-core AssaySpec module must emit sample_qc"),
        ("emit: crispr_qc_report", "nf-core AssaySpec module must emit crispr_qc_report"),
        ("emit: crispr_qc_json", "nf-core AssaySpec module must emit crispr_qc_json"),
        ("emit: crispr_qc_summary", "nf-core AssaySpec module must emit crispr_qc_summary"),
        ("versions.yml", "nf-core AssaySpec module must emit versions.yml"),
    ]:
        _require(assay_module, needle, message, result)
    for needle in [
        "name: dotmatch_assay_run",
        "Run a DotMatch AssaySpec",
        "- assayspec",
        "assay_report:",
        "assay_manifest:",
        "assay_manifest_summary:",
        "sample_qc:",
        "crispr_qc_report:",
        "crispr_qc_json:",
        "crispr_qc_summary:",
        "versions:",
    ]:
        _require(assay_meta, needle, f"nf-core AssaySpec module metadata missing {needle}", result)
    if not assay_nf_test:
        result.failures.append("nf-core AssaySpec module must include an nf-test candidate at tests/main.nf.test")
    else:
        for needle, message in [
            ("nextflow_process", "nf-core AssaySpec nf-test candidate must define a nextflow_process"),
            ('script "../main.nf"', "nf-core AssaySpec nf-test candidate must reference ../main.nf"),
            ('process "DOTMATCH_ASSAY_RUN"', "nf-core AssaySpec nf-test candidate must test DOTMATCH_ASSAY_RUN"),
            ("examples/workflows/fixtures/crispr_assay.toml", "nf-core AssaySpec nf-test candidate must use shared AssaySpec fixture"),
            ("examples/workflows/fixtures/crispr_library.csv", "nf-core AssaySpec nf-test candidate must stage target table"),
            ("examples/workflows/fixtures/sample_a.fastq", "nf-core AssaySpec nf-test candidate must stage FASTQ inputs"),
            ("assay_report", "nf-core AssaySpec nf-test candidate must assert assay_report output"),
            ("assay_manifest_summary", "nf-core AssaySpec nf-test candidate must assert assay_manifest_summary output"),
            ("sample_qc", "nf-core AssaySpec nf-test candidate must assert sample_qc output"),
            ("crispr_qc_summary", "nf-core AssaySpec nf-test candidate must assert crispr_qc_summary output"),
            ("crispr_qc_report", "nf-core AssaySpec nf-test candidate must assert crispr_qc_report output"),
        ]:
            _require(assay_nf_test, needle, message, result)

    if not any("nf-core" in failure for failure in result.failures):
        result.passed.append("nf-core-style module candidate present without adoption claim")


def check_multiqc(root: Path, result: WorkflowAudit) -> None:
    config = _read(root / "examples" / "workflows" / "multiqc" / "multiqc_config.yaml", result)
    sample_qc_path = root / "examples" / "workflows" / "multiqc" / "data" / "sample_qc.tsv"
    sample_qc = _read(sample_qc_path, result)

    for needle in [
        "custom_data:",
        "dotmatch_sample_qc:",
        "dotmatch_assay_manifest:",
        "dotmatch_crispr_qc:",
        'plot_type: "table"',
        'fn: "*sample_qc.tsv"',
        'fn: "*assay_manifest.summary.tsv"',
        'fn: "*crispr_qc.summary.tsv"',
        "assignment_rate:",
        "ambiguous_rate:",
        "no_match_rate:",
        "primary_report:",
        "autopsy_triggered:",
        "warning_count:",
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

    manifest_summary_path = root / "examples" / "workflows" / "multiqc" / "data" / "assay_manifest.summary.tsv"
    manifest_summary = _read(manifest_summary_path, result)
    manifest_header = manifest_summary.splitlines()[0].split("\t") if manifest_summary.splitlines() else []
    for column in [
        "schema_version",
        "mode",
        "assay_type",
        "status",
        "autopsy_triggered",
        "warning_count",
        "production_warning_count",
        "sample_count",
        "primary_report",
        "manifest",
    ]:
        if column not in manifest_header:
            result.failures.append(f"MultiQC assay_manifest.summary.tsv missing {column}")

    crispr_qc_path = root / "examples" / "workflows" / "multiqc" / "data" / "crispr_qc.summary.tsv"
    crispr_qc = _read(crispr_qc_path, result)
    crispr_qc_header = crispr_qc.splitlines()[0].split("\t") if crispr_qc.splitlines() else []
    for column in [
        "sample_id",
        "qc_status",
        "total_count",
        "coverage_fraction",
        "zero_count_fraction",
        "gini_index",
        "top_1pct_fraction",
        "assignment_rate",
        "ambiguous_rate",
        "no_match_rate",
        "invalid_rate",
    ]:
        if column not in crispr_qc_header:
            result.failures.append(f"MultiQC crispr_qc.summary.tsv missing {column}")

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
    _require(readme, "not been published to a Galaxy ToolShed", "Galaxy README must avoid ToolShed release claims", result)
    _require(readme, "AssaySpec", "Galaxy README must describe the AssaySpec example wrapper", result)
    _require(readme, "planemo", "Galaxy README must mention planemo linting", result)

    assay_wrapper_path = root / "examples" / "workflows" / "galaxy" / "dotmatch_assay_run.xml"
    try:
        assay_wrapper = ET.parse(assay_wrapper_path).getroot()
    except Exception as exc:
        result.failures.append(f"Galaxy AssaySpec wrapper XML could not be parsed: {exc}")
        assay_wrapper = None
    if assay_wrapper is not None:
        if assay_wrapper.tag != "tool" or assay_wrapper.attrib.get("id") != "dotmatch_assay_run":
            result.failures.append("Galaxy AssaySpec wrapper must be tool id dotmatch_assay_run")
        assay_command = assay_wrapper.findtext("command") or ""
        _require(assay_command, "cat > assay.toml", "Galaxy AssaySpec wrapper command must generate an AssaySpec from staged inputs", result)
        _require(assay_command, "dotmatch assay run assay.toml", "Galaxy AssaySpec wrapper command must run dotmatch assay run", result)
        assay_input_names = {node.attrib.get("name", "") for node in assay_wrapper.findall("./inputs/param")}
        required_inputs = {"library", "sample1_fastq", "sample1_label", "sample2_fastq", "sample2_label", "guide_start", "guide_length", "k", "metric", "ambiguous"}
        if not required_inputs <= assay_input_names:
            result.failures.append("Galaxy AssaySpec wrapper inputs must stage library, FASTQs, labels, window, metric, k, and ambiguity policy")
        assay_output_names = {node.attrib.get("name", "") for node in assay_wrapper.findall("./outputs/data")}
        required_outputs = {
            "assay_report",
            "assay_manifest",
            "assay_manifest_summary",
            "sample_qc",
            "crispr_qc_report",
            "crispr_qc_json",
            "crispr_qc_summary",
            "counts",
            "summary",
        }
        if not required_outputs <= assay_output_names:
            result.failures.append("Galaxy AssaySpec wrapper outputs must include report, manifest, manifest summary, sample QC, CRISPR QC, counts, and summary")
        assay_test = assay_wrapper.find("./tests/test")
        if assay_test is None:
            result.failures.append("Galaxy AssaySpec wrapper must include a Planemo test with tiny AssaySpec fixtures")
        else:
            params = {node.attrib.get("name", ""): node.attrib.get("value", "") for node in assay_test.findall("param")}
            for name, value in [
                ("library", "crispr_library.csv"),
                ("sample1_fastq", "sample_a.fastq"),
                ("sample1_label", "sample_a"),
                ("sample2_fastq", "sample_b.fastq"),
                ("sample2_label", "sample_b"),
            ]:
                if params.get(name) != value:
                    result.failures.append(f"Galaxy AssaySpec Planemo test must set {name}={value}")
            test_outputs = {node.attrib.get("name", ""): node for node in assay_test.findall("output")}
            report = test_outputs.get("assay_report")
            if report is None or report.find("./assert_contents/has_text[@text='DotMatch Assay Report']") is None:
                result.failures.append("Galaxy AssaySpec Planemo test must assert DotMatch Assay Report content")
            manifest_summary = test_outputs.get("assay_manifest_summary")
            if manifest_summary is None or manifest_summary.find("./assert_contents/has_text[@text='primary_report']") is None:
                result.failures.append("Galaxy AssaySpec Planemo test must assert manifest summary content")
            crispr_qc_summary = test_outputs.get("crispr_qc_summary")
            if crispr_qc_summary is None or crispr_qc_summary.find("./assert_contents/has_text[@text='qc_status']") is None:
                result.failures.append("Galaxy AssaySpec Planemo test must assert CRISPR QC summary content")
            crispr_qc_report = test_outputs.get("crispr_qc_report")
            if crispr_qc_report is None or crispr_qc_report.find("./assert_contents/has_text[@text='DotMatch CRISPR QC']") is None:
                result.failures.append("Galaxy AssaySpec Planemo test must assert CRISPR QC report content")

    if not any("Galaxy" in failure for failure in result.failures):
        result.passed.append("Galaxy wrapper example present without ToolShed claim")


def check_workflow_fixtures(root: Path, result: WorkflowAudit) -> None:
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
        result.passed.append("shared workflow test fixtures present")


def audit(root: Path) -> WorkflowAudit:
    root = root.resolve()
    result = WorkflowAudit()
    check_workflow_fixtures(root, result)
    check_snakemake(root, result)
    check_nextflow(root, result)
    check_nfcore(root, result)
    check_multiqc(root, result)
    check_galaxy(root, result)
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
