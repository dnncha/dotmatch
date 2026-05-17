from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_snakemake_crispr_config_is_complete() -> None:
    config_path = ROOT / "examples" / "workflows" / "snakemake" / "config.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))

    assert config["library"] == "examples/crispr_guides/data/yusa_library.csv"
    assert config["samples"] == {
        "plasmid": "examples/crispr_guides/data/ERR376998.fastq.gz",
        "ESC1": "examples/crispr_guides/data/ERR376999.fastq.gz",
    }
    assert config["guide_start"] == 23
    assert config["guide_length"] == 19
    assert config["metric"] in {"hamming", "levenshtein"}
    assert config["outdir"] == "examples/workflows/snakemake/output"


def test_snakemake_example_runs_dotmatch_crispr_count() -> None:
    snakefile = ROOT / "examples" / "workflows" / "snakemake" / "Snakefile"
    text = snakefile.read_text(encoding="utf-8")

    assert "rule dotmatch_crispr_count" in text
    assert "dotmatch crispr-count" in text
    assert "--format" not in text
    assert "--ambiguous discard" in text
    assert "--sample-qc" in text
    assert "sample_qc" in text


def test_nextflow_crispr_config_is_complete() -> None:
    config_path = ROOT / "examples" / "workflows" / "nextflow" / "nextflow.config"
    config = config_path.read_text(encoding="utf-8")

    assert "library = 'examples/crispr_guides/data/yusa_library.csv'" in config
    assert "samples = 'examples/workflows/nextflow/samples.tsv'" in config
    assert "guide_start = 23" in config
    assert "guide_length = 19" in config
    assert "metric = 'levenshtein'" in config
    assert "outdir = 'examples/workflows/nextflow/output'" in config


def test_nextflow_example_runs_dotmatch_crispr_count() -> None:
    workflow = ROOT / "examples" / "workflows" / "nextflow" / "main.nf"
    text = workflow.read_text(encoding="utf-8")

    assert "nextflow.enable.dsl=2" in text
    assert "process DOTMATCH_CRISPR_COUNT" in text
    assert "dotmatch crispr-count" in text
    assert "--ambiguous discard" in text
    assert "--sample-qc" in text
    assert 'path "sample_qc.tsv", emit: sample_qc' in text
    assert "publishDir params.outdir" in text


def test_nfcore_module_candidate_has_expected_metadata() -> None:
    meta_path = ROOT / "examples" / "workflows" / "nf-core" / "modules" / "local" / "dotmatch" / "crispr_count" / "meta.yml"
    text = meta_path.read_text(encoding="utf-8")

    assert "name: dotmatch_crispr_count" in text
    assert "description: Count CRISPR guides with DotMatch" in text
    assert "- dotmatch" in text
    assert "- crispr" in text
    assert "counts:" in text
    assert "summary:" in text
    assert "sample_qc:" in text
    assert "versions:" in text


def test_nfcore_module_candidate_runs_dotmatch_crispr_count() -> None:
    module_path = ROOT / "examples" / "workflows" / "nf-core" / "modules" / "local" / "dotmatch" / "crispr_count" / "main.nf"
    text = module_path.read_text(encoding="utf-8")

    assert "process DOTMATCH_CRISPR_COUNT" in text
    assert "tuple val(meta), path(reads), path(library)" in text
    assert "dotmatch crispr-count" in text
    assert "--ambiguous discard" in text
    assert "--sample-qc" in text
    assert "emit: sample_qc" in text
    assert "versions.yml" in text
    assert "dotmatch --version" in text
    assert "importlib.metadata" not in text
    assert "task.ext.args" in text


def test_nfcore_module_candidate_has_nf_test_fixture() -> None:
    test_path = (
        ROOT
        / "examples"
        / "workflows"
        / "nf-core"
        / "modules"
        / "local"
        / "dotmatch"
        / "crispr_count"
        / "tests"
        / "main.nf.test"
    )
    text = test_path.read_text(encoding="utf-8")

    assert "nextflow_process" in text
    assert "process \"DOTMATCH_CRISPR_COUNT\"" in text
    assert "examples/workflows/fixtures/crispr_library.csv" in text
    assert "sample_qc" in text
    assert "versions.yml" in text


def test_multiqc_config_targets_dotmatch_sample_qc() -> None:
    config_path = ROOT / "examples" / "workflows" / "multiqc" / "multiqc_config.yaml"
    config = config_path.read_text(encoding="utf-8")

    assert "custom_data:" in config
    assert "dotmatch_sample_qc:" in config
    assert "dotmatch_crispr_qc:" in config
    assert 'plot_type: "table"' in config
    assert 'fn: "*sample_qc.tsv"' in config
    assert 'fn: "*crispr_qc.summary.tsv"' in config
    assert "assignment_rate:" in config
    assert "ambiguous_rate:" in config
    assert "no_match_rate:" in config


def test_multiqc_example_fixture_matches_dotmatch_schema() -> None:
    fixture_path = ROOT / "examples" / "workflows" / "multiqc" / "data" / "sample_qc.tsv"
    lines = fixture_path.read_text(encoding="utf-8").splitlines()

    header = lines[0].split("\t")
    assert header[:6] == [
        "sample_id",
        "fastq",
        "total_reads",
        "valid_extracted_reads",
        "assigned_reads",
        "exact_reads",
    ]
    assert "assignment_rate" in header
    assert "ambiguous_rate" in header
    assert "candidates_verified" in header
    assert len(lines) == 3

    crispr_qc = ROOT / "examples" / "workflows" / "multiqc" / "data" / "crispr_qc.summary.tsv"
    crispr_header = crispr_qc.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert crispr_header == [
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
    ]


def test_galaxy_wrapper_has_dotmatch_crispr_count_surface() -> None:
    wrapper_path = ROOT / "examples" / "workflows" / "galaxy" / "dotmatch_crispr_count.xml"
    root = ET.parse(wrapper_path).getroot()

    assert root.tag == "tool"
    assert root.attrib["id"] == "dotmatch_crispr_count"

    requirement_names = [node.text for node in root.findall("./requirements/requirement")]
    assert "dotmatch" in requirement_names

    command = root.findtext("command") or ""
    assert "dotmatch crispr-count" in command
    assert "--ambiguous" in command
    assert "--summary" in command
    assert "--sample-qc" in command

    output_names = {node.attrib["name"] for node in root.findall("./outputs/data")}
    assert {"counts", "summary", "sample_qc"} <= output_names


def test_galaxy_wrapper_has_planemo_fixture_test() -> None:
    wrapper_path = ROOT / "examples" / "workflows" / "galaxy" / "dotmatch_crispr_count.xml"
    root = ET.parse(wrapper_path).getroot()
    test = root.find("./tests/test")

    assert test is not None
    params = {node.attrib["name"]: node.attrib.get("value", "") for node in test.findall("param")}
    assert params["library"] == "crispr_library.csv"
    assert params["sample1_fastq"] == "sample_a.fastq"
    assert params["sample2_fastq"] == "sample_b.fastq"
    outputs = {node.attrib["name"]: node.attrib.get("file", "") for node in test.findall("output")}
    assert outputs["counts"] == "expected_counts.mageck.tsv"
    sample_qc = next(node for node in test.findall("output") if node.attrib["name"] == "sample_qc")
    assert sample_qc.find("./assert_contents/has_text[@text='assignment_rate']") is not None
    assert sample_qc.find("./assert_contents/has_text[@text='sample_a']") is not None


def test_workflow_fixtures_cover_core_outcomes() -> None:
    fixtures = ROOT / "examples" / "workflows" / "fixtures"
    readme = (fixtures / "README.md").read_text(encoding="utf-8")
    library = (fixtures / "crispr_library.csv").read_text(encoding="utf-8")
    sample_a = (fixtures / "sample_a.fastq").read_text(encoding="utf-8")
    expected_qc = (fixtures / "expected_sample_qc.tsv").read_text(encoding="utf-8")

    for outcome in ["unique", "ambiguous", "unmatched", "invalid"]:
        assert outcome in readme
        assert outcome in sample_a or outcome in expected_qc
    assert "id,gRNA.sequence,Gene" in library
    assert "assignment_rate" in expected_qc
