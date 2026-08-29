#!/usr/bin/env python3
"""Audit local workflow-manager examples for DotMatch."""

from __future__ import annotations

import argparse
import json
import xml.etree.ElementTree as ET
from pathlib import Path


WORKFLOW_FIXTURES = [
    "README.md",
    "crispr_assay.toml",
    "crispr_library.csv",
    "barcodes.tsv",
    "barcode_reads.fastq",
    "panel_barcodes.tsv",
    "sample_a.fastq",
    "sample_b.fastq",
    "expected_counts.mageck.tsv",
    "expected_sample_qc.tsv",
]
GALAXY_TEST_DATA = [
    "crispr_library.csv",
    "barcodes.tsv",
    "barcode_reads.fastq",
    "panel_barcodes.tsv",
    "sample_a.fastq",
    "sample_b.fastq",
    "expected_counts.mageck.tsv",
]
NFCORE_MODULES = ["count", "demux", "audit", "panel_check", "crispr_count", "assay_run"]
NFCORE_CONTAINER_TAG = "0.2.2--py311h13f8228_1"


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

    required_keys = [
        "library",
        "samples",
        "guide_start",
        "guide_length",
        "metric",
        "ambiguity_policy",
        "ambiguous",
        "outdir",
        "barcode_table",
        "barcode_reads",
        "panel_table",
    ]
    for key in required_keys:
        if key not in config:
            result.failures.append(f"Snakemake config.json missing {key}")
    if config.get("metric") not in {"hamming", "levenshtein"}:
        result.failures.append("Snakemake config.json metric must be hamming or levenshtein")
    if config.get("ambiguity_policy") not in {"radius", "best"}:
        result.failures.append("Snakemake config.json ambiguity_policy must be radius or best")
    if config.get("ambiguous") not in {"discard", "include", "separate"}:
        result.failures.append("Snakemake config.json ambiguous must be discard, include, or separate")
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
    _require(snakefile, "ambiguity_policy = config.get", "Snakemake Snakefile must read ambiguity policy from config", result)
    _require(snakefile, "ambiguous = config.get", "Snakemake Snakefile must read ambiguous-output handling from config", result)
    _require(snakefile, "--ambiguity-policy {params.ambiguity_policy}", "Snakemake Snakefile must keep assignment ambiguity policy explicit", result)
    _require(snakefile, "--ambiguous {params.ambiguous}", "Snakemake Snakefile must keep ambiguous-output handling explicit", result)
    _require(snakefile, "--sample-qc", "Snakemake Snakefile must emit sample_qc.tsv for MultiQC", result)
    _require(snakefile, "sample_qc", "Snakemake Snakefile must declare sample_qc output", result)
    _require(snakefile, "rule dotmatch_demux", "Snakemake Snakefile missing rule dotmatch_demux", result)
    _require(snakefile, "dotmatch demux", "Snakemake demux rule must run dotmatch demux", result)
    _require(snakefile, "--assignments {output.assignments}", "Snakemake demux rule must emit assignments", result)
    _require(snakefile, "rule dotmatch_panel_check", "Snakemake Snakefile missing rule dotmatch_panel_check", result)
    _require(snakefile, "dotmatch panel check", "Snakemake panel-check rule must run dotmatch panel check", result)
    _require(snakefile, "panel_summary.json", "Snakemake panel-check rule must emit panel_summary.json", result)

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
        "k = 1",
        "metric = 'levenshtein'",
        "ambiguity_policy = 'radius'",
        "ambiguous = 'discard'",
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
        ("--ambiguity-policy ${params.ambiguity_policy}", "Nextflow workflow must keep assignment ambiguity policy explicit"),
        ("--ambiguous ${params.ambiguous}", "Nextflow workflow must keep ambiguous-output handling explicit"),
        ("--sample-qc", "Nextflow workflow must emit sample_qc.tsv for MultiQC"),
        ("path \"sample_qc.tsv\", emit: sample_qc", "Nextflow workflow must declare sample_qc output"),
        ("publishDir params.outdir", "Nextflow workflow must publish outputs to params.outdir"),
    ]:
        _require(workflow, needle, message, result)

    if not any("Nextflow" in failure for failure in result.failures):
        result.passed.append("Nextflow CRISPR workflow example present")


def check_nfcore(root: Path, result: WorkflowAudit) -> None:
    base = root / "examples" / "workflows" / "nf-core"
    if not (base / "README.md").is_file():
        result.failures.append("nf-core README.md missing")
    module = _read(base / "modules" / "local" / "dotmatch" / "crispr_count" / "main.nf", result)
    meta = _read(base / "modules" / "local" / "dotmatch" / "crispr_count" / "meta.yml", result)
    test_path = base / "modules" / "local" / "dotmatch" / "crispr_count" / "tests" / "main.nf.test"
    nf_test = _read(test_path, result) if test_path.is_file() else ""

    for needle, message in [
        ("process DOTMATCH_CRISPR_COUNT", "nf-core module missing DOTMATCH_CRISPR_COUNT process"),
        ("tuple val(meta), path(reads), path(library)", "nf-core module missing expected input tuple"),
        ("dotmatch crispr-count", "nf-core module must run dotmatch crispr-count"),
        ("--ambiguity-policy radius", "nf-core module must keep assignment ambiguity policy explicit"),
        ("--ambiguous discard", "nf-core module must keep ambiguous-output handling explicit"),
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
            ("${System.getenv('PWD')}/../../../../../fixtures/crispr_library.csv", "nf-core nf-test candidate must use shared workflow fixtures"),
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
            ("${System.getenv('PWD')}/../../../../../fixtures/crispr_assay.toml", "nf-core AssaySpec nf-test candidate must use shared AssaySpec fixture"),
            ("${System.getenv('PWD')}/../../../../../fixtures/crispr_library.csv", "nf-core AssaySpec nf-test candidate must stage target table"),
            ("${System.getenv('PWD')}/../../../../../fixtures/sample_a.fastq", "nf-core AssaySpec nf-test candidate must stage FASTQ inputs"),
            ("assay_report", "nf-core AssaySpec nf-test candidate must assert assay_report output"),
            ("assay_manifest_summary", "nf-core AssaySpec nf-test candidate must assert assay_manifest_summary output"),
            ("sample_qc", "nf-core AssaySpec nf-test candidate must assert sample_qc output"),
            ("crispr_qc_summary", "nf-core AssaySpec nf-test candidate must assert crispr_qc_summary output"),
            ("crispr_qc_report", "nf-core AssaySpec nf-test candidate must assert crispr_qc_report output"),
        ]:
            _require(assay_nf_test, needle, message, result)

    if not any("nf-core" in failure for failure in result.failures):
        result.passed.append("nf-core-style module candidate present")


def check_nfcore_container_pins(root: Path, result: WorkflowAudit) -> None:
    base = root / "examples" / "workflows" / "nf-core"
    module_roots = [
        base / "modules" / "local" / "dotmatch",
        base / "upstream" / "modules" / "nf-core" / "dotmatch",
    ]
    expected_singularity = f"https://depot.galaxyproject.org/singularity/dotmatch:{NFCORE_CONTAINER_TAG}"
    expected_docker = f"biocontainers/dotmatch:{NFCORE_CONTAINER_TAG}"

    for module_root in module_roots:
        for module_name in NFCORE_MODULES:
            module_path = module_root / module_name / "main.nf"
            module = _read(module_path, result)
            if expected_singularity not in module or expected_docker not in module:
                result.failures.append(
                    f"{module_path.as_posix()} must use the verified immutable nf-core container tag "
                    f"{NFCORE_CONTAINER_TAG} for both Singularity and Docker"
                )

    if not any("immutable nf-core container tag" in failure for failure in result.failures):
        result.passed.append(f"nf-core modules use immutable container tag {NFCORE_CONTAINER_TAG}")


def _check_nfcore_tool_module(
    root: Path,
    result: WorkflowAudit,
    module_name: str,
    process_name: str,
    command: str,
    required_module_needles: list[str],
    required_meta_needles: list[str],
    required_test_needles: list[str],
) -> None:
    base = root / "examples" / "workflows" / "nf-core" / "modules" / "local" / "dotmatch" / module_name
    module = _read(base / "main.nf", result)
    meta = _read(base / "meta.yml", result)
    test_path = base / "tests" / "main.nf.test"
    nf_test = _read(test_path, result) if test_path.is_file() else ""
    for needle, message in [
        (f"process {process_name}", f"nf-core {module_name} module missing {process_name} process"),
        (command, f"nf-core {module_name} module must run {command}"),
        ("task.ext.args", f"nf-core {module_name} module must expose task.ext.args"),
        ("versions.yml", f"nf-core {module_name} module must emit versions.yml"),
        ("dotmatch --version", f"nf-core {module_name} module must record dotmatch --version"),
    ]:
        _require(module, needle, message, result)
    for needle in required_module_needles:
        _require(module, needle, f"nf-core {module_name} module missing {needle}", result)
    for needle in required_meta_needles:
        _require(meta, needle, f"nf-core {module_name} metadata missing {needle}", result)
    if not nf_test:
        result.failures.append(f"nf-core {module_name} module must include an nf-test candidate at tests/main.nf.test")
    else:
        for needle, message in [
            ("nextflow_process", f"nf-core {module_name} nf-test candidate must define a nextflow_process"),
            ('script "../main.nf"', f"nf-core {module_name} nf-test candidate must reference ../main.nf"),
            (f'process "{process_name}"', f"nf-core {module_name} nf-test candidate must test {process_name}"),
            ("Channel.of", f"nf-core {module_name} nf-test candidate must build input channels"),
        ]:
            _require(nf_test, needle, message, result)
        for needle in required_test_needles:
            _require(nf_test, needle, f"nf-core {module_name} nf-test candidate missing {needle}", result)


def check_nfcore_dotmatch_modules(root: Path, result: WorkflowAudit) -> None:
    _check_nfcore_tool_module(
        root,
        result,
        "count",
        "DOTMATCH_COUNT",
        "dotmatch count",
        [
            "tuple val(meta), path(reads), path(targets)",
            "--target-start ${target_start}",
            "--target-length ${target_length}",
            "--sample-qc",
            "--assignments",
            "emit: sample_qc",
            "emit: assignments",
        ],
        [
            "name: dotmatch_count",
            "Count fixed-window known targets with DotMatch",
            "counts:",
            "sample_qc:",
            "assignments:",
            "versions:",
        ],
        [
            "${System.getenv('PWD')}/../../../../../fixtures/barcode_reads.fastq",
            "${System.getenv('PWD')}/../../../../../fixtures/barcodes.tsv",
            "sample_qc",
            "assignments",
        ],
    )
    _check_nfcore_tool_module(
        root,
        result,
        "demux",
        "DOTMATCH_DEMUX",
        "dotmatch demux",
        [
            "tuple val(meta), path(reads), path(barcodes)",
            "--barcode-start ${barcode_start}",
            "--barcode-length ${barcode_length}",
            "--out-dir demuxed",
            "--assignments",
            "emit: demuxed",
            "emit: assignments",
        ],
        [
            "name: dotmatch_demux",
            "Demultiplex fixed-window inline barcodes with DotMatch",
            "demuxed:",
            "summary:",
            "assignments:",
            "versions:",
        ],
        [
            "${System.getenv('PWD')}/../../../../../fixtures/barcode_reads.fastq",
            "${System.getenv('PWD')}/../../../../../fixtures/barcodes.tsv",
            "demuxed",
            "assigned_unique",
        ],
    )
    _check_nfcore_tool_module(
        root,
        result,
        "audit",
        "DOTMATCH_AUDIT",
        "dotmatch audit",
        [
            "tuple val(meta), path(targets)",
            "--audit-mode ${audit_mode}",
            "--out-dir audit",
            "emit: target_safety",
            "emit: collision_pairs",
        ],
        [
            "name: dotmatch_audit",
            "Audit a DotMatch target library",
            "audit_dir:",
            "target_safety:",
            "collision_pairs:",
            "versions:",
        ],
        [
            "${System.getenv('PWD')}/../../../../../fixtures/barcodes.tsv",
            "target_safety",
            "min_hamming_distance",
        ],
    )
    _check_nfcore_tool_module(
        root,
        result,
        "panel_check",
        "DOTMATCH_PANEL_CHECK",
        "dotmatch panel check",
        [
            "tuple val(meta), path(panel)",
            "--metric ${metric}",
            "--out-dir panel_check",
            "emit: panel_summary",
            "emit: target_safety",
            "emit: collision_pairs",
        ],
        [
            "name: dotmatch_panel_check",
            "Check a barcode panel for DotMatch assignment safety",
            "panel_check_dir:",
            "panel_summary:",
            "target_safety:",
            "versions:",
        ],
        [
            "${System.getenv('PWD')}/../../../../../fixtures/panel_barcodes.tsv",
            "panel_summary",
            "target_safety",
        ],
    )
    if not any("nf-core count" in failure or "nf-core demux" in failure or "nf-core audit" in failure or "nf-core panel_check" in failure for failure in result.failures):
        result.passed.append("nf-core DotMatch count/demux/audit/panel_check modules present")


def check_multiqc(root: Path, result: WorkflowAudit) -> None:
    config = _read(root / "examples" / "workflows" / "multiqc" / "multiqc_config.yaml", result)
    data_dir = root / "examples" / "workflows" / "multiqc" / "data"
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

    summary_path = data_dir / "assignment_summary.json"
    try:
        summary = json.loads(_read(summary_path, result))
    except json.JSONDecodeError as exc:
        result.failures.append(f"MultiQC assignment_summary.json is invalid JSON: {exc}")
        summary = {}
    for field in [
        "sample_label",
        "total_reads",
        "assigned_unique",
        "ambiguous",
        "unmatched",
        "assignment_rate",
    ]:
        if field not in summary:
            result.failures.append(f"MultiQC assignment_summary.json missing {field}")

    top_unmatched = _read(data_dir / "sample_top_unmatched.tsv", result)
    top_unmatched_header = top_unmatched.splitlines()[0].split("\t") if top_unmatched.splitlines() else []
    for column in ["sequence", "count"]:
        if column not in top_unmatched_header:
            result.failures.append(f"MultiQC sample_top_unmatched.tsv missing {column}")

    if not any("MultiQC" in failure for failure in result.failures):
        result.passed.append("MultiQC custom-content example present")


def check_galaxy(root: Path, result: WorkflowAudit) -> None:
    if not (root / "examples" / "workflows" / "galaxy" / "README.md").is_file():
        result.failures.append("Galaxy README.md missing")
    test_data = root / "examples" / "workflows" / "galaxy" / "test-data"
    wrapper_path = root / "examples" / "workflows" / "galaxy" / "dotmatch_crispr_count.xml"
    try:
        wrapper = ET.parse(wrapper_path).getroot()
    except Exception as exc:
        result.failures.append(f"Galaxy wrapper XML could not be parsed: {exc}")
        return

    if wrapper.tag != "tool" or wrapper.attrib.get("id") != "dotmatch_crispr_count":
        result.failures.append("Galaxy wrapper must be tool id dotmatch_crispr_count")
    if wrapper.attrib.get("version") != "0.2.1+galaxy0":
        result.failures.append("Galaxy CRISPR wrapper must track the public Bioconda 0.2.1 package")
    command = wrapper.findtext("command") or ""
    _require(command, "dotmatch crispr-count", "Galaxy wrapper command must run dotmatch crispr-count", result)
    _require(command, "--ambiguity-policy radius", "Galaxy wrapper command must keep assignment ambiguity policy explicit", result)
    _require(command, "--ambiguous", "Galaxy wrapper command must expose --ambiguous", result)
    _require(command, "--summary", "Galaxy wrapper command must include --summary", result)
    _require(command, "--sample-qc", "Galaxy wrapper command must include --sample-qc", result)
    _require(command, "element_identifier", "Galaxy wrapper command must derive sample IDs from Galaxy datasets", result)
    _require(command, "ln -s", "Galaxy wrapper command must stage input FASTQs", result)
    requirements = {node.text: node.attrib.get("version", "") for node in wrapper.findall("./requirements/requirement")}
    if requirements.get("dotmatch") != "0.2.1":
        result.failures.append("Galaxy wrapper must require public Bioconda dotmatch=0.2.1")
    reads = wrapper.find("./inputs/param[@name='reads']")
    if reads is None or reads.attrib.get("multiple") != "true":
        result.failures.append("Galaxy wrapper must accept one or more FASTQ datasets through reads")
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
            ("reads", "sample_a.fastq,sample_b.fastq"),
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
    expected_counts = test_data / "expected_counts.mageck.tsv"
    if expected_counts.is_file() and "guide_a\tGENEA\t0\t0" not in expected_counts.read_text(encoding="utf-8"):
        result.failures.append("Galaxy expected counts must match the pinned dotmatch=0.2.1 guide_a assignment")
    for filename in GALAXY_TEST_DATA:
        if not (test_data / filename).is_file():
            result.failures.append(f"Galaxy Planemo test-data file is missing: {filename}")
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
        _require(assay_command, "printf '%s\\n'", "Galaxy AssaySpec wrapper command must generate an AssaySpec from staged inputs", result)
        _require(assay_command, 'ambiguity_policy = "radius"', "Galaxy AssaySpec wrapper command must keep assignment ambiguity policy explicit", result)
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

    demux_wrapper_path = root / "examples" / "workflows" / "galaxy" / "dotmatch_demux.xml"
    try:
        demux_wrapper = ET.parse(demux_wrapper_path).getroot()
    except Exception as exc:
        result.failures.append(f"Galaxy demux wrapper XML could not be parsed: {exc}")
        demux_wrapper = None
    if demux_wrapper is not None:
        if demux_wrapper.tag != "tool" or demux_wrapper.attrib.get("id") != "dotmatch_demux":
            result.failures.append("Galaxy demux wrapper must be tool id dotmatch_demux")
        demux_command = demux_wrapper.findtext("command") or ""
        _require(demux_command, "dotmatch demux", "Galaxy demux wrapper command must run dotmatch demux", result)
        _require(demux_command, "--summary", "Galaxy demux wrapper command must include --summary", result)
        _require(demux_command, "--assignments", "Galaxy demux wrapper command must include --assignments", result)
        demux_outputs = {node.attrib.get("name", "") for node in demux_wrapper.findall("./outputs/data")}
        if not {"demuxed_fastqs", "summary", "assignments"} <= demux_outputs:
            result.failures.append("Galaxy demux wrapper outputs must include demuxed_fastqs, summary, and assignments")
        demux_test = demux_wrapper.find("./tests/test")
        if demux_test is None:
            result.failures.append("Galaxy demux wrapper must include a Planemo test")
        else:
            params = {node.attrib.get("name", ""): node.attrib.get("value", "") for node in demux_test.findall("param")}
            if params.get("barcodes") != "barcodes.tsv" or params.get("reads") != "barcode_reads.fastq":
                result.failures.append("Galaxy demux Planemo test must use barcode fixtures")

    panel_wrapper_path = root / "examples" / "workflows" / "galaxy" / "dotmatch_panel_check.xml"
    try:
        panel_wrapper = ET.parse(panel_wrapper_path).getroot()
    except Exception as exc:
        result.failures.append(f"Galaxy panel-check wrapper XML could not be parsed: {exc}")
        panel_wrapper = None
    if panel_wrapper is not None:
        if panel_wrapper.tag != "tool" or panel_wrapper.attrib.get("id") != "dotmatch_panel_check":
            result.failures.append("Galaxy panel-check wrapper must be tool id dotmatch_panel_check")
        panel_command = panel_wrapper.findtext("command") or ""
        _require(panel_command, "dotmatch panel check", "Galaxy panel-check wrapper command must run dotmatch panel check", result)
        _require(panel_command, "panel_summary.json", "Galaxy panel-check wrapper command must expose panel_summary.json", result)
        panel_outputs = {node.attrib.get("name", "") for node in panel_wrapper.findall("./outputs/data")}
        if not {"panel_summary", "target_safety", "collision_pairs", "panel_report"} <= panel_outputs:
            result.failures.append("Galaxy panel-check wrapper outputs must include summary, safety, collision, and report files")
        panel_test = panel_wrapper.find("./tests/test")
        if panel_test is None:
            result.failures.append("Galaxy panel-check wrapper must include a Planemo test")
        else:
            params = {node.attrib.get("name", ""): node.attrib.get("value", "") for node in panel_test.findall("param")}
            if params.get("panel") != "panel_barcodes.tsv":
                result.failures.append("Galaxy panel-check Planemo test must use panel_barcodes.tsv")

    if not any("Galaxy" in failure for failure in result.failures):
        result.passed.append("Galaxy wrapper example present")


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
    check_nfcore_container_pins(root, result)
    check_nfcore_dotmatch_modules(root, result)
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
