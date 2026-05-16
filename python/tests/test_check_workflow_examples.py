import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_workflow_examples.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_workflow_examples", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_workflow_repo(root: Path) -> None:
    files = {
        "examples/workflows/snakemake/config.json": (
            '{"library": "examples/crispr_guides/data/yusa_library.csv", '
            '"samples": {"plasmid": "examples/crispr_guides/data/ERR376998.fastq.gz"}, '
            '"guide_start": 23, "guide_length": 19, "metric": "levenshtein", '
            '"outdir": "examples/workflows/snakemake/output"}\n'
        ),
        "examples/workflows/snakemake/Snakefile": (
            "rule all:\n"
            "    input: 'counts.tsv', 'assay_out/assay_report.html', 'assay_out/assay_manifest.summary.tsv'\n"
            "rule dotmatch_assay_spec:\n"
            "    output: 'assay.toml'\n"
            "    shell: 'dotmatch assay init --template crispr --out {output}'\n"
            "rule dotmatch_assay_run:\n"
            "    input: spec='assay.toml'\n"
            "    output: report='assay_out/assay_report.html', manifest='assay_out/assay_manifest.json', manifest_summary='assay_out/assay_manifest.summary.tsv', sample_qc='assay_out/sample_qc.tsv'\n"
            "    shell: 'dotmatch assay run {input.spec}'\n"
            "rule dotmatch_crispr_count:\n"
            "    output: counts='counts.mageck.tsv', summary='summary.json', sample_qc='sample_qc.tsv'\n"
            "    shell: 'dotmatch crispr-count --ambiguous discard --summary {output.summary} --sample-qc {output.sample_qc}'\n"
        ),
        "examples/workflows/nextflow/nextflow.config": (
            "params {\n"
            "  library = 'examples/crispr_guides/data/yusa_library.csv'\n"
            "  samples = 'examples/workflows/nextflow/samples.tsv'\n"
            "  guide_start = 23\n"
            "  guide_length = 19\n"
            "  metric = 'levenshtein'\n"
            "  outdir = 'examples/workflows/nextflow/output'\n"
            "}\n"
        ),
        "examples/workflows/nextflow/main.nf": (
            "nextflow.enable.dsl=2\n"
            "process DOTMATCH_CRISPR_COUNT {\n"
            "  publishDir params.outdir\n"
            "  output:\n"
            "  path \"sample_qc.tsv\", emit: sample_qc\n"
            "  script:\n"
            "  \"\"\"\n"
            "  dotmatch crispr-count --ambiguous discard --summary summary.json --sample-qc sample_qc.tsv\n"
            "  \"\"\"\n"
            "}\n"
            "process DOTMATCH_ASSAY_RUN {\n"
            "  publishDir params.outdir\n"
            "  output:\n"
            "  path \"assay_report.html\", emit: assay_report\n"
            "  path \"assay_manifest.json\", emit: assay_manifest\n"
            "  path \"assay_manifest.summary.tsv\", emit: assay_manifest_summary\n"
            "  path \"sample_qc.tsv\", emit: sample_qc\n"
            "  script:\n"
            "  \"\"\"\n"
            "  dotmatch assay run assay.toml\n"
            "  cp assay_out/assay_report.html assay_report.html\n"
            "  cp assay_out/assay_manifest.json assay_manifest.json\n"
            "  cp assay_out/assay_manifest.summary.tsv assay_manifest.summary.tsv\n"
            "  cp assay_out/sample_qc.tsv sample_qc.tsv\n"
            "  \"\"\"\n"
            "}\n"
        ),
        "examples/workflows/nf-core/README.md": (
            "# nf-core-style Module Candidate\n\n"
            "This is a local nf-core-style module candidate. It has not been submitted to or "
            "accepted by nf-core. It includes a dotmatch_assay_run AssaySpec module candidate.\n"
        ),
        "examples/workflows/nf-core/modules/local/dotmatch/crispr_count/main.nf": (
            "process DOTMATCH_CRISPR_COUNT {\n"
            "  input:\n"
            "  tuple val(meta), path(reads), path(library)\n"
            "  output:\n"
            "  path 'versions.yml'\n"
            "  tuple val(meta), path('sample_qc.tsv'), emit: sample_qc\n"
            "  script:\n"
            "  \"\"\"\n"
            "  dotmatch crispr-count --ambiguous discard --summary summary.json --sample-qc sample_qc.tsv\n"
            "  dotmatch --version > versions.yml\n"
            "  $task.ext.args\n"
            "  \"\"\"\n"
            "}\n"
        ),
        "examples/workflows/nf-core/modules/local/dotmatch/crispr_count/meta.yml": (
            "name: dotmatch_crispr_count\n"
            "description: Count CRISPR guides with DotMatch\n"
            "keywords:\n"
            "  - dotmatch\n"
            "  - crispr\n"
            "output:\n"
            "  - counts:\n"
            "  - summary:\n"
            "  - sample_qc:\n"
            "  - versions:\n"
        ),
        "examples/workflows/nf-core/modules/local/dotmatch/crispr_count/tests/main.nf.test": (
            "nextflow_process {\n"
            "  name \"Test Process DOTMATCH_CRISPR_COUNT\"\n"
            "  script \"../main.nf\"\n"
            "  process \"DOTMATCH_CRISPR_COUNT\"\n"
            "  test(\"runs tiny CRISPR fixture\") {\n"
            "    when {\n"
            "      process {\n"
            "        \"\"\"\n"
            "        input[0] = Channel.of([ [id:'fixture'], file('examples/workflows/fixtures/sample_a.fastq'), file('examples/workflows/fixtures/crispr_library.csv') ])\n"
            "        input[1] = 0\n"
            "        input[2] = 4\n"
            "        input[3] = 1\n"
            "        input[4] = 'hamming'\n"
            "        \"\"\"\n"
            "      }\n"
            "    }\n"
            "    then {\n"
            "      assert process.success\n"
            "      assert path(process.out.counts[0][1]).name.endsWith('.counts.mageck.tsv')\n"
            "      assert path(process.out.sample_qc[0][1]).text.contains('assignment_rate')\n"
            "      assert path(process.out.versions[0]).name == 'versions.yml'\n"
            "      assert path(process.out.versions[0]).text.contains('dotmatch')\n"
            "    }\n"
            "  }\n"
            "}\n"
        ),
        "examples/workflows/nf-core/modules/local/dotmatch/assay_run/main.nf": (
            "process DOTMATCH_ASSAY_RUN {\n"
            "  input:\n"
            "  tuple val(meta), path(assay_spec), path(assay_inputs)\n"
            "  output:\n"
            "  tuple val(meta), path('assay_report.html'), emit: assay_report\n"
            "  tuple val(meta), path('assay_manifest.json'), emit: assay_manifest\n"
            "  tuple val(meta), path('assay_manifest.summary.tsv'), emit: assay_manifest_summary\n"
            "  tuple val(meta), path('sample_qc.tsv'), emit: sample_qc\n"
            "  path 'versions.yml', emit: versions\n"
            "  script:\n"
            "  \"\"\"\n"
            "  basename\n"
            "  dotmatch assay run ${assay_spec}\n"
            "  cp assay_out/assay_report.html assay_report.html\n"
            "  cp assay_out/assay_manifest.json assay_manifest.json\n"
            "  cp assay_out/assay_manifest.summary.tsv assay_manifest.summary.tsv\n"
            "  cp assay_out/sample_qc.tsv sample_qc.tsv\n"
            "  dotmatch --version > versions.yml\n"
            "  \"\"\"\n"
            "}\n"
        ),
        "examples/workflows/nf-core/modules/local/dotmatch/assay_run/meta.yml": (
            "name: dotmatch_assay_run\n"
            "description: Run a DotMatch AssaySpec and emit workflow-facing reports\n"
            "keywords:\n"
            "  - dotmatch\n"
            "  - assayspec\n"
            "output:\n"
            "  - assay_report:\n"
            "  - assay_manifest:\n"
            "  - assay_manifest_summary:\n"
            "  - sample_qc:\n"
            "  - versions:\n"
        ),
        "examples/workflows/nf-core/modules/local/dotmatch/assay_run/tests/main.nf.test": (
            "nextflow_process {\n"
            "  name \"Test Process DOTMATCH_ASSAY_RUN\"\n"
            "  script \"../main.nf\"\n"
            "  process \"DOTMATCH_ASSAY_RUN\"\n"
            "  test(\"runs tiny AssaySpec fixture\") {\n"
            "    when { process { \"\"\" input[0] = Channel.of([ [id:'fixture'], file('examples/workflows/fixtures/crispr_assay.toml'), [file('examples/workflows/fixtures/crispr_library.csv'), file('examples/workflows/fixtures/sample_a.fastq')] ]) \"\"\" } }\n"
            "    then {\n"
            "      assert process.success\n"
            "      assert path(process.out.assay_report[0][1]).text.contains('DotMatch Assay Report')\n"
            "      assert path(process.out.assay_manifest_summary[0][1]).text.contains('primary_report')\n"
            "      assert path(process.out.sample_qc[0][1]).text.contains('assignment_rate')\n"
            "    }\n"
            "  }\n"
            "}\n"
        ),
        "examples/workflows/multiqc/multiqc_config.yaml": (
            "custom_data:\n"
            "  dotmatch_sample_qc:\n"
            "    file_format: tsv\n"
            "    plot_type: \"table\"\n"
            "    fn: \"*sample_qc.tsv\"\n"
            "    headers:\n"
            "      assignment_rate:\n"
            "      ambiguous_rate:\n"
            "      no_match_rate:\n"
            "  dotmatch_assay_manifest:\n"
            "    file_format: tsv\n"
            "    plot_type: \"table\"\n"
            "    fn: \"*assay_manifest.summary.tsv\"\n"
            "    headers:\n"
            "      primary_report:\n"
            "      autopsy_triggered:\n"
            "      warning_count:\n"
        ),
        "examples/workflows/multiqc/data/sample_qc.tsv": (
            "sample_id\tfastq\ttotal_reads\tvalid_extracted_reads\tassigned_reads\texact_reads\t"
            "assignment_rate\tambiguous_rate\tno_match_rate\tcandidates_verified\n"
            "plasmid\treads.fastq.gz\t10\t10\t9\t9\t0.9\t0.0\t0.1\t10\n"
        ),
        "examples/workflows/multiqc/data/assay_manifest.summary.tsv": (
            "schema_version\tmode\tassay_type\tstatus\tnative_version\tautopsy_triggered\twarning_count\t"
            "production_warning_count\tsample_count\tprimary_report\tmanifest\n"
            "1\tcount\tcrispr\tready\tdotmatch 0.1.0\tfalse\t0\t0\t1\tassay_report.html\tassay_manifest.json\n"
        ),
        "examples/workflows/galaxy/README.md": (
            "# Galaxy Wrapper\n\nThese are example wrappers, including an AssaySpec example wrapper, not been published to a Galaxy ToolShed. Use planemo lint.\n"
        ),
        "examples/workflows/galaxy/dotmatch_crispr_count.xml": (
            "<tool id=\"dotmatch_crispr_count\">\n"
            "  <requirements><requirement type=\"package\">dotmatch</requirement></requirements>\n"
            "  <command>dotmatch crispr-count --ambiguous discard --summary '$summary' --sample-qc '$sample_qc'</command>\n"
            "  <outputs><data name=\"counts\"/><data name=\"summary\"/><data name=\"sample_qc\"/></outputs>\n"
            "  <tests><test><param name=\"library\" value=\"crispr_library.csv\"/><param name=\"sample1_fastq\" value=\"sample_a.fastq\"/><param name=\"sample1_label\" value=\"sample_a\"/><param name=\"sample2_fastq\" value=\"sample_b.fastq\"/><param name=\"sample2_label\" value=\"sample_b\"/><output name=\"counts\" file=\"expected_counts.mageck.tsv\"/><output name=\"sample_qc\"><assert_contents><has_text text=\"assignment_rate\"/><has_text text=\"sample_a\"/><has_text text=\"sample_b\"/></assert_contents></output></test></tests>\n"
            "</tool>\n"
        ),
        "examples/workflows/galaxy/dotmatch_assay_run.xml": (
            "<tool id=\"dotmatch_assay_run\">\n"
            "  <requirements><requirement type=\"package\">dotmatch</requirement></requirements>\n"
            "  <command>cat > assay.toml\n dotmatch assay run assay.toml</command>\n"
            "  <inputs><param name=\"library\"/><param name=\"sample1_fastq\"/><param name=\"sample1_label\"/><param name=\"sample2_fastq\"/><param name=\"sample2_label\"/><param name=\"guide_start\"/><param name=\"guide_length\"/><param name=\"k\"/><param name=\"metric\"/><param name=\"ambiguous\"/></inputs>\n"
            "  <outputs><data name=\"assay_report\" from_work_dir=\"assay_out/assay_report.html\"/><data name=\"assay_manifest\" from_work_dir=\"assay_out/assay_manifest.json\"/><data name=\"assay_manifest_summary\" from_work_dir=\"assay_out/assay_manifest.summary.tsv\"/><data name=\"sample_qc\" from_work_dir=\"assay_out/sample_qc.tsv\"/><data name=\"counts\" from_work_dir=\"assay_out/counts.mageck.tsv\"/><data name=\"summary\" from_work_dir=\"assay_out/summary.json\"/></outputs>\n"
            "  <tests><test><param name=\"library\" value=\"crispr_library.csv\"/><param name=\"sample1_fastq\" value=\"sample_a.fastq\"/><param name=\"sample1_label\" value=\"sample_a\"/><param name=\"sample2_fastq\" value=\"sample_b.fastq\"/><param name=\"sample2_label\" value=\"sample_b\"/><output name=\"assay_report\"><assert_contents><has_text text=\"DotMatch Assay Report\"/></assert_contents></output><output name=\"assay_manifest_summary\"><assert_contents><has_text text=\"primary_report\"/></assert_contents></output></test></tests>\n"
            "</tool>\n"
        ),
        "examples/workflows/fixtures/README.md": (
            "# Workflow Test Fixtures\n\n"
            "Tiny CRISPR fixtures exercise unique, ambiguous, unmatched, and invalid DotMatch outcomes.\n"
        ),
        "examples/workflows/fixtures/crispr_library.csv": (
            "id,gRNA.sequence,Gene\n"
            "guide_a,ACGT,GENEA\n"
            "guide_b,ACGA,GENEB\n"
            "guide_c,TTTT,GENEC\n"
        ),
        "examples/workflows/fixtures/sample_a.fastq": (
            "@unique\nACGT\n+\nIIII\n@ambiguous\nACGG\n+\nIIII\n@unmatched\nCCCC\n+\nIIII\n@invalid\nAC\n+\nII\n"
        ),
        "examples/workflows/fixtures/sample_b.fastq": (
            "@unique_b\nTTTT\n+\nIIII\n@unmatched_b\nGGGG\n+\nIIII\n"
        ),
        "examples/workflows/fixtures/expected_counts.mageck.tsv": (
            "sgRNA\tGene\tsample_a\tsample_b\n"
            "guide_a\tGENEA\t1\t0\n"
            "guide_b\tGENEB\t0\t0\n"
            "guide_c\tGENEC\t0\t1\n"
        ),
        "examples/workflows/fixtures/expected_sample_qc.tsv": (
            "sample_id\tfastq\ttotal_reads\tvalid_extracted_reads\tassigned_reads\texact_reads\t"
            "k1_rescued_reads\tk1_sub_reads\tk1_ins_reads\tk1_del_reads\tambiguous_reads\t"
            "no_match_reads\tinvalid_reads\tassignment_rate\texact_rate\trescue_rate\tambiguous_rate\t"
            "no_match_rate\ttargets_observed\tzero_count_targets\tgini_index\ttop_1pct_read_fraction\t"
            "candidates_verified\n"
            "sample_a\tsample_a.fastq\t4\t3\t1\t1\t0\t0\t0\t0\t1\t1\t1\t0.25000000\t"
            "0.25000000\t0.00000000\t0.25000000\t0.25000000\t1\t2\t-0.66666667\t"
            "1.00000000\t3\n"
            "sample_b\tsample_b.fastq\t2\t2\t1\t1\t0\t0\t0\t0\t0\t1\t0\t0.50000000\t"
            "0.50000000\t0.00000000\t0.00000000\t0.50000000\t1\t2\t-0.66666667\t"
            "1.00000000\t1\n"
        ),
        "examples/workflows/fixtures/crispr_assay.toml": (
            "schema_version = 1\n"
            "mode = \"count\"\n"
            "assay_type = \"crispr\"\n"
            "targets = \"examples/workflows/fixtures/crispr_library.csv\"\n"
            "[run]\n"
            "out_dir = \"assay_out\"\n"
            "[extract]\n"
            "start = 0\n"
            "length = 4\n"
            "[assignment]\n"
            "k = 1\n"
            "metric = \"hamming\"\n"
            "[[samples]]\n"
            "id = \"sample_a\"\n"
            "fastq = \"examples/workflows/fixtures/sample_a.fastq\"\n"
        ),
        "examples/workflows/galaxy/test-data/crispr_library.csv": (
            "id,gRNA.sequence,Gene\n"
            "guide_a,ACGT,GENEA\n"
            "guide_b,ACGA,GENEB\n"
            "guide_c,TTTT,GENEC\n"
        ),
        "examples/workflows/galaxy/test-data/sample_a.fastq": (
            "@unique\nACGT\n+\nIIII\n@ambiguous\nACGG\n+\nIIII\n@unmatched\nCCCC\n+\nIIII\n@invalid\nAC\n+\nII\n"
        ),
        "examples/workflows/galaxy/test-data/sample_b.fastq": (
            "@unique_b\nTTTT\n+\nIIII\n@unmatched_b\nGGGG\n+\nIIII\n"
        ),
        "examples/workflows/galaxy/test-data/expected_counts.mageck.tsv": (
            "sgRNA\tGene\tsample_a\tsample_b\n"
            "guide_a\tGENEA\t1\t0\n"
            "guide_b\tGENEB\t0\t0\n"
            "guide_c\tGENEC\t0\t1\n"
        ),
        "examples/workflows/galaxy/test-data/crispr_assay.toml": (
            "schema_version = 1\n"
            "mode = \"count\"\n"
            "assay_type = \"crispr\"\n"
            "targets = \"crispr_library.csv\"\n"
            "[run]\n"
            "out_dir = \"assay_out\"\n"
            "[extract]\n"
            "start = 0\n"
            "length = 4\n"
            "[assignment]\n"
            "k = 1\n"
            "metric = \"hamming\"\n"
            "[[samples]]\n"
            "id = \"sample_a\"\n"
            "fastq = \"sample_a.fastq\"\n"
        ),
    }
    for path, text in files.items():
        full = root / path
        full.parent.mkdir(parents=True, exist_ok=True)
        full.write_text(text, encoding="utf-8")


def test_workflow_examples_ready_accepts_complete_local_examples(tmp_path):
    checker = _load_checker()
    _write_workflow_repo(tmp_path)

    result = checker.audit(tmp_path)

    assert result.failures == []
    assert any("Snakemake" in item for item in result.passed)
    assert any("Nextflow" in item for item in result.passed)
    assert any("nf-core" in item for item in result.passed)
    assert any("MultiQC" in item for item in result.passed)
    assert any("Galaxy" in item for item in result.passed)


def test_workflow_examples_ready_requires_assayspec_workflow_pack(tmp_path):
    checker = _load_checker()
    _write_workflow_repo(tmp_path)
    (tmp_path / "examples" / "workflows" / "galaxy" / "dotmatch_assay_run.xml").unlink()
    (tmp_path / "examples" / "workflows" / "multiqc" / "data" / "assay_manifest.summary.tsv").unlink()

    result = checker.audit(tmp_path)

    assert any("AssaySpec" in failure and "Galaxy" in failure for failure in result.failures)
    assert any("assay_manifest.summary.tsv" in failure for failure in result.failures)


def test_workflow_examples_ready_rejects_external_adoption_claim(tmp_path):
    checker = _load_checker()
    _write_workflow_repo(tmp_path)
    (tmp_path / "examples" / "workflows" / "nf-core" / "README.md").write_text(
        "# nf-core Module\n\nThis is now an upstream nf-core module with external adoption.\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("external-adoption" in failure for failure in result.failures)


def test_workflow_examples_ready_rejects_missing_multiqc_schema_columns(tmp_path):
    checker = _load_checker()
    _write_workflow_repo(tmp_path)
    (tmp_path / "examples" / "workflows" / "multiqc" / "data" / "sample_qc.tsv").write_text(
        "sample_id\tfastq\ttotal_reads\nplasmid\treads.fastq.gz\t10\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("sample_qc.tsv" in failure and "assigned_reads" in failure for failure in result.failures)


def test_workflow_examples_ready_rejects_weak_galaxy_wrapper(tmp_path):
    checker = _load_checker()
    _write_workflow_repo(tmp_path)
    (tmp_path / "examples" / "workflows" / "galaxy" / "dotmatch_crispr_count.xml").write_text(
        "<tool id=\"dotmatch_crispr_count\"><command>dotmatch crispr-count</command></tool>\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("Galaxy" in failure and "--summary" in failure for failure in result.failures)


def test_workflow_examples_ready_requires_snakemake_sample_qc_output(tmp_path):
    checker = _load_checker()
    _write_workflow_repo(tmp_path)
    snakefile = (tmp_path / "examples" / "workflows" / "snakemake" / "Snakefile").read_text(encoding="utf-8")
    (tmp_path / "examples" / "workflows" / "snakemake" / "Snakefile").write_text(
        snakefile.replace(" --sample-qc {output.sample_qc}", "").replace(", sample_qc='sample_qc.tsv'", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("Snakemake" in failure and "sample_qc" in failure for failure in result.failures)


def test_workflow_examples_ready_requires_nextflow_sample_qc_output(tmp_path):
    checker = _load_checker()
    _write_workflow_repo(tmp_path)
    workflow = (tmp_path / "examples" / "workflows" / "nextflow" / "main.nf").read_text(encoding="utf-8")
    (tmp_path / "examples" / "workflows" / "nextflow" / "main.nf").write_text(
        workflow.replace(" --sample-qc sample_qc.tsv", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("Nextflow" in failure and "sample_qc" in failure for failure in result.failures)


def test_workflow_examples_ready_requires_nfcore_sample_qc_output(tmp_path):
    checker = _load_checker()
    _write_workflow_repo(tmp_path)
    module = (
        tmp_path
        / "examples"
        / "workflows"
        / "nf-core"
        / "modules"
        / "local"
        / "dotmatch"
        / "crispr_count"
        / "main.nf"
    ).read_text(encoding="utf-8")
    (
        tmp_path
        / "examples"
        / "workflows"
        / "nf-core"
        / "modules"
        / "local"
        / "dotmatch"
        / "crispr_count"
        / "main.nf"
    ).write_text(
        module.replace("  path 'sample_qc.tsv'\n", "").replace(" --sample-qc sample_qc.tsv", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("nf-core" in failure and "sample_qc" in failure for failure in result.failures)


def test_workflow_examples_ready_requires_galaxy_sample_qc_output(tmp_path):
    checker = _load_checker()
    _write_workflow_repo(tmp_path)
    wrapper_path = tmp_path / "examples" / "workflows" / "galaxy" / "dotmatch_crispr_count.xml"
    wrapper = wrapper_path.read_text(encoding="utf-8")
    wrapper_path.write_text(
        wrapper.replace(" --sample-qc '$sample_qc'", "").replace("<data name=\"sample_qc\"/>", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("Galaxy" in failure and "sample_qc" in failure for failure in result.failures)


def test_workflow_examples_ready_requires_shared_workflow_fixtures(tmp_path):
    checker = _load_checker()
    _write_workflow_repo(tmp_path)
    (tmp_path / "examples" / "workflows" / "fixtures" / "sample_a.fastq").unlink()
    (tmp_path / "examples" / "workflows" / "fixtures" / "expected_sample_qc.tsv").unlink()

    result = checker.audit(tmp_path)

    assert any("workflow test fixture is missing" in failure and "sample_a.fastq" in failure for failure in result.failures)
    assert any("workflow test fixture is missing" in failure and "expected_sample_qc.tsv" in failure for failure in result.failures)


def test_workflow_examples_ready_requires_galaxy_planemo_test_block(tmp_path):
    checker = _load_checker()
    _write_workflow_repo(tmp_path)
    wrapper_path = tmp_path / "examples" / "workflows" / "galaxy" / "dotmatch_crispr_count.xml"
    wrapper = wrapper_path.read_text(encoding="utf-8")
    start = wrapper.index("  <tests>")
    end = wrapper.index("</tests>") + len("</tests>\n")
    wrapper_path.write_text(wrapper[:start] + wrapper[end:], encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("Galaxy wrapper must include a Planemo test" in failure for failure in result.failures)


def test_workflow_examples_ready_requires_nfcore_nf_test_candidate(tmp_path):
    checker = _load_checker()
    _write_workflow_repo(tmp_path)
    (
        tmp_path
        / "examples"
        / "workflows"
        / "nf-core"
        / "modules"
        / "local"
        / "dotmatch"
        / "crispr_count"
        / "tests"
        / "main.nf.test"
    ).unlink()

    result = checker.audit(tmp_path)

    assert any("nf-core module must include an nf-test candidate" in failure for failure in result.failures)
