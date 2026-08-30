from __future__ import annotations

import json
import gzip
import os
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]


def _write_count_spec(tmp_path: Path) -> Path:
    out_dir = tmp_path / "assay_out"
    spec = tmp_path / "crispr.toml"
    spec.write_text(
        f"""
schema_version = 1
mode = "count"
assay_type = "crispr"
targets = "{ROOT / 'examples/workflows/fixtures/crispr_library.csv'}"

[[samples]]
id = "sample_a"
fastq = "{ROOT / 'examples/workflows/fixtures/sample_a.fastq'}"

[[samples]]
id = "sample_b"
fastq = "{ROOT / 'examples/workflows/fixtures/sample_b.fastq'}"

[run]
out_dir = "{out_dir}"
threads = 1

[extract]
start = 0
length = 4

[assignment]
k = 1
metric = "hamming"
ambiguous = "discard"

[reliability]
fail_on_unsafe_targets = false

[outputs]
format = "mageck"
assignments = true
ambiguous = true
unmatched = true
""".lstrip(),
        encoding="utf-8",
    )
    return spec


def _write_demux_spec(tmp_path: Path) -> Path:
    barcodes = tmp_path / "barcodes.tsv"
    reads = tmp_path / "reads.fastq"
    barcodes.write_text("bc0\tACGT\nbc1\tTTTT\nbc2\tAGGA\n", encoding="utf-8")
    reads.write_text(
        "@d0\nACGTAAAA\n+\nIIIIIIII\n"
        "@d1\nTTTGAAAA\n+\nIIIIIIII\n"
        "@d2\nAGGAAAAA\n+\nIIIIIIII\n",
        encoding="utf-8",
    )
    spec = tmp_path / "demux.toml"
    spec.write_text(
        f"""
schema_version = 1
mode = "demux"
assay_type = "inline_barcode"
barcodes = "{barcodes}"
reads = "{reads}"

[run]
out_dir = "{tmp_path / 'demux_out'}"

[extract]
start = 0
length = 4

[assignment]
k = 1
metric = "hamming"

[reliability]
fail_on_unsafe_targets = false

[outputs]
assignments = true
ambiguous = true
unmatched = true
""".lstrip(),
        encoding="utf-8",
    )
    return spec


def _write_pair_spec(tmp_path: Path) -> Path:
    left = tmp_path / "left.tsv"
    right = tmp_path / "right.tsv"
    reads = tmp_path / "pair.fastq"
    left.write_text("L0\tACGT\nL1\tTTTT\n", encoding="utf-8")
    right.write_text("R0\tGGAA\nR1\tCCCC\n", encoding="utf-8")
    reads.write_text(
        "@p0\nACGTGGAA\n+\nIIIIIIII\n"
        "@p1\nTTTTCCCC\n+\nIIIIIIII\n",
        encoding="utf-8",
    )
    spec = tmp_path / "pair.toml"
    spec.write_text(
        f"""
schema_version = 1
mode = "pair-count"
assay_type = "generic"
left_targets = "{left}"
right_targets = "{right}"
reads = "{reads}"

[run]
out_dir = "{tmp_path / 'pair_out'}"

[left]
start = 0
length = 4

[right]
start = 4
length = 4

[assignment]
k = 1
metric = "hamming"

[outputs]
assignments = true
""".lstrip(),
        encoding="utf-8",
    )
    return spec


def _write_paired_pair_spec(tmp_path: Path) -> Path:
    left = tmp_path / "paired_left.tsv"
    right = tmp_path / "paired_right.tsv"
    left_reads = tmp_path / "pair_R1.fastq"
    right_reads = tmp_path / "pair_R2.fastq"
    left.write_text("L0\tACGT\nL1\tTTTT\n", encoding="utf-8")
    right.write_text("R0\tGGAA\nR1\tCCCC\n", encoding="utf-8")
    left_reads.write_text(
        "@p0/1\nACGT\n+\nIIII\n"
        "@p1 1:N:0:1\nTTTT\n+\nIIII\n",
        encoding="utf-8",
    )
    right_reads.write_text(
        "@p0/2\nGGAA\n+\nIIII\n"
        "@p1 2:N:0:1\nCCCC\n+\nIIII\n",
        encoding="utf-8",
    )
    spec = tmp_path / "paired_pair.toml"
    spec.write_text(
        f"""
schema_version = 1
mode = "pair-count"
assay_type = "generic"
left_targets = "{left}"
right_targets = "{right}"
left_reads = "{left_reads}"
right_reads = "{right_reads}"

[run]
out_dir = "{tmp_path / 'paired_pair_out'}"

[left]
start = 0
length = 4

[right]
start = 0
length = 4

[assignment]
k = 1
metric = "hamming"

[outputs]
assignments = true
""".lstrip(),
        encoding="utf-8",
    )
    return spec


def _write_inference_targets(tmp_path: Path) -> Path:
    targets = tmp_path / "targets.tsv"
    targets.write_text("guide_a\tACGT\tGENEA\nguide_b\tTTTT\tGENEB\n", encoding="utf-8")
    return targets


def _write_inference_reads(tmp_path: Path, *, prefix: str = "NN", good: bool = True) -> Path:
    reads = tmp_path / "infer.fastq"
    if good:
        reads.write_text(
            f"@r0\n{prefix}ACGTAAAA\n+\nIIIIIIIIII\n"
            f"@r1\n{prefix}TTTTAAAA\n+\nIIIIIIIIII\n"
            f"@r2\n{prefix}ACGTCCCC\n+\nIIIIIIIIII\n"
            f"@r3\n{prefix}TTTTCCCC\n+\nIIIIIIIIII\n",
            encoding="utf-8",
        )
    else:
        reads.write_text(
            "@r0\nGGGGAAAA\n+\nIIIIIIII\n"
            "@r1\nCCCCAAAA\n+\nIIIIIIII\n"
            "@r2\nAAAACCCC\n+\nIIIIIIII\n",
            encoding="utf-8",
        )
    return reads


def _write_wrong_offset_spec(tmp_path: Path) -> Path:
    targets = _write_inference_targets(tmp_path)
    reads = _write_inference_reads(tmp_path, prefix="NN", good=True)
    spec = tmp_path / "wrong_offset.toml"
    spec.write_text(
        f"""
schema_version = 1
mode = "count"
assay_type = "crispr"
targets = "{targets}"

[[samples]]
id = "shifted"
fastq = "{reads}"

[run]
out_dir = "{tmp_path / 'wrong_offset_out'}"
threads = 1

[extract]
start = 0
length = 4

[assignment]
k = 1
metric = "hamming"
ambiguous = "discard"

[outputs]
format = "mageck"
""".lstrip(),
        encoding="utf-8",
    )
    return spec


def _write_unsafe_count_spec(tmp_path: Path, *, profile: str) -> Path:
    targets = tmp_path / f"unsafe_targets_{profile}.tsv"
    reads = tmp_path / f"unsafe_reads_{profile}.fastq"
    out_dir = tmp_path / f"unsafe_{profile}_out"
    targets.write_text("g0\tACGT\ng1\tACGA\n", encoding="utf-8")
    reads.write_text(
        "@u0\nACGTAAAA\n+\nIIIIIIII\n"
        "@u1\nACGAAAAA\n+\nIIIIIIII\n",
        encoding="utf-8",
    )
    spec = tmp_path / f"unsafe_{profile}.toml"
    spec.write_text(
        f"""
schema_version = 1
mode = "count"
assay_type = "crispr"
targets = "{targets}"

[[samples]]
id = "unsafe"
fastq = "{reads}"

[run]
out_dir = "{out_dir}"
threads = 1

[extract]
start = 0
length = 4

[assignment]
k = 1
metric = "hamming"
ambiguous = "discard"

[reliability]
profile = "{profile}"
fail_on_unsafe_targets = true

[outputs]
format = "mageck"
assignments = true
""".lstrip(),
        encoding="utf-8",
    )
    return spec


def _write_non_acgt_count_spec(tmp_path: Path) -> Path:
    targets = tmp_path / "iupac_targets.tsv"
    reads = tmp_path / "iupac_reads.fastq"
    targets.write_text("g0\tACGN\ng1\tTTTT\n", encoding="utf-8")
    reads.write_text(
        "@n0\nACGNAAAA\n+\nIIIIIIII\n"
        "@n1\nTTTTAAAA\n+\nIIIIIIII\n",
        encoding="utf-8",
    )
    spec = tmp_path / "iupac_count.toml"
    spec.write_text(
        f"""
schema_version = 1
mode = "count"
assay_type = "generic"
targets = "{targets}"

[[samples]]
id = "iupac"
fastq = "{reads}"

[run]
out_dir = "{tmp_path / 'iupac_out'}"
threads = 8

[extract]
start = 0
length = 4

[assignment]
k = 1
metric = "hamming"
ambiguous = "discard"

[reliability]
fail_on_unsafe_targets = false

[outputs]
assignments = true
""".lstrip(),
        encoding="utf-8",
    )
    return spec


def _run_cli(args: list[str], *, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    return subprocess.run(
        [sys.executable, "-m", "dotmatch.cli", *args],
        cwd=ROOT,
        env=merged_env,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def test_load_count_spec_and_compile_deterministic_plan(tmp_path: Path) -> None:
    from dotmatch.assayspec import compile_assay_plan, load_assay_spec

    assay = load_assay_spec(_write_count_spec(tmp_path))
    plan = compile_assay_plan(assay)

    assert [step.name for step in plan.steps] == ["audit", "run", "crispr-qc", "validate"]
    assert plan.steps[0].argv[:3] == ["dotmatch-native", "audit", "--targets"]
    assert plan.steps[1].argv[:2] == ["dotmatch-native", "crispr-count"]
    assert plan.steps[2].argv[:2] == ["dotmatch", "crispr-qc"]
    assert "--ambiguity-policy" in plan.steps[1].argv
    assert "radius" in plan.steps[1].argv
    assert "--sample-qc" in plan.steps[1].argv
    assert "--target-counts-long" in plan.steps[1].argv
    assert "--format" not in plan.steps[1].argv
    assert plan.artifacts["counts"].name == "counts.mageck.tsv"
    assert plan.artifacts["reliability_summary"].name == "reliability_summary.json"
    assert plan.artifacts["reliability_findings"].name == "reliability_findings.tsv"
    assert plan.artifacts["reliability_report"].name == "reliability_report.html"
    assert plan.artifacts["reliability_manifest_summary"].name == "reliability_manifest.summary.tsv"
    assert plan.artifacts["methods"].name == "methods.md"
    assert plan.artifacts["citation_bib"].name == "CITATION.bib"
    assert plan.artifacts["software_versions"].name == "software_versions.yml"


def test_compile_assay_plan_passes_backend_mode(tmp_path: Path) -> None:
    from dotmatch.assayspec import compile_assay_plan, load_assay_spec

    spec = tmp_path / "backend_mode.toml"
    spec.write_text(
        _write_count_spec(tmp_path).read_text(encoding="utf-8")
        + """
[backend]
mode = "gpu-metal-experimental"
""",
        encoding="utf-8",
    )
    plan = compile_assay_plan(load_assay_spec(spec))
    run_argv = plan.steps[1].argv
    assert "--backend" in run_argv
    assert "gpu-metal-experimental" in run_argv


def test_compile_assay_plan_adds_metal_validate_for_gpu_backend(tmp_path: Path) -> None:
    from dotmatch.assayspec import compile_assay_plan, load_assay_spec

    spec = tmp_path / "backend_gpu.toml"
    spec.write_text(
        _write_count_spec(tmp_path).read_text(encoding="utf-8")
        + """
[backend]
mode = "gpu-metal-experimental"
""",
        encoding="utf-8",
    )
    plan = compile_assay_plan(load_assay_spec(spec))
    run_argv = plan.steps[1].argv
    assert "--backend" in run_argv
    assert "gpu-metal-experimental" in run_argv
    assert "--metal-validate" in run_argv


def test_compile_assay_plan_forces_cpu_when_gpu_disallowed(tmp_path: Path) -> None:
    from dotmatch.assayspec import compile_assay_plan, load_assay_spec

    spec = tmp_path / "backend_cpu.toml"
    spec.write_text(
        _write_count_spec(tmp_path).read_text(encoding="utf-8")
        + """
[backend]
allow_gpu = false
""",
        encoding="utf-8",
    )
    plan = compile_assay_plan(load_assay_spec(spec))
    run_argv = plan.steps[1].argv
    assert "--backend" in run_argv
    assert "cpu" in run_argv


def test_assay_check_rejects_invalid_enum(tmp_path: Path) -> None:
    spec = _write_count_spec(tmp_path)
    spec.write_text(spec.read_text(encoding="utf-8").replace('metric = "hamming"', 'metric = "jaccard"'), encoding="utf-8")

    rc = _run_cli(["assay", "check", str(spec)])

    assert rc.returncode == 2
    assert "assignment.metric" in rc.stderr
    assert "hamming" in rc.stderr


def test_assay_plan_prints_native_commands_without_creating_outputs(tmp_path: Path) -> None:
    spec = _write_count_spec(tmp_path)

    rc = _run_cli(["assay", "plan", str(spec)])

    assert rc.returncode == 0, rc.stderr
    assert "dotmatch-native audit --targets" in rc.stdout
    assert "dotmatch-native crispr-count --library" in rc.stdout
    assert "# reliability_report:" in rc.stdout
    assert "# methods:" in rc.stdout
    assert "# citation_bib:" in rc.stdout
    assert "# software_versions:" in rc.stdout
    assert not (tmp_path / "assay_out").exists()


def test_assay_check_writes_preflight_reliability_artifacts(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)
    spec.write_text(spec.read_text(encoding="utf-8").replace("k = 1", "k = 0"), encoding="utf-8")

    rc = _run_cli(["assay", "check", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 0, rc.stderr
    assert "check passed" in rc.stderr
    assert f"{spec.name}: ok" in rc.stdout
    assert "next:" in rc.stderr
    out_dir = tmp_path / "assay_out"
    summary = json.loads((out_dir / "reliability_summary.json").read_text(encoding="utf-8"))
    assert summary["stage"] == "preflight"
    assert summary["overall_status"] == "passed"
    assert summary["profile"] == "production"
    assert summary["backend"]["authority"] == "cpu"
    assert summary["backend"]["gpu_status"] in {"eligible_but_not_used", "not_eligible"}
    assert summary["backend_optimizer"]["authority"] == "cpu"
    assert summary["backend_optimizer"]["candidate_backend"] in {"gpu-metal-experimental", "cpu"}
    assert summary["evidence_boundary"]["status"] == "supported"
    assert "checked public" in summary["evidence_boundary"]["claim_boundary"]
    assert any(finding["finding_id"] == "read_qc_unavailable" for finding in summary["findings"])

    findings = (out_dir / "reliability_findings.tsv").read_text(encoding="utf-8")
    assert "finding_id\tseverity\tstage\tsample_id\tmetric\tobserved\tthreshold\tmessage\trecommended_action\tsource_artifact" in findings
    assert "read_qc_unavailable" in findings

    manifest_summary = (out_dir / "reliability_manifest.summary.tsv").read_text(encoding="utf-8")
    assert "overall_status\tprofile\tfinding_count\tblocked_count\terror_count\twarning_count" in manifest_summary

    report = (out_dir / "reliability_report.html").read_text(encoding="utf-8")
    assert "<title>DotMatch Reliability Report</title>" in report
    assert "Evidence Boundary" in report
    assert "read_qc_unavailable" in report


def test_assay_check_ignores_stale_audit_artifacts(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)
    spec.write_text(spec.read_text(encoding="utf-8").replace("k = 1", "k = 0"), encoding="utf-8")
    stale_audit = tmp_path / "assay_out" / "audit" / "audit_summary.json"
    stale_audit.parent.mkdir(parents=True)
    stale_audit.write_text('{"safe_at_k1": false}\n', encoding="utf-8")

    rc = _run_cli(["assay", "check", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 0, rc.stderr
    reliability = json.loads((tmp_path / "assay_out" / "reliability_summary.json").read_text(encoding="utf-8"))
    assert not any(finding["finding_id"] == "unsafe_targets" for finding in reliability["findings"])


def test_assay_check_ignores_stale_postrun_qc_artifacts(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)
    spec.write_text(spec.read_text(encoding="utf-8").replace("k = 1", "k = 0"), encoding="utf-8")
    out_dir = tmp_path / "assay_out"
    out_dir.mkdir(parents=True)
    (out_dir / "sample_qc.tsv").write_text(
        "sample_id\tassignment_rate\tambiguous_rate\tno_match_rate\ttotal_reads\tinvalid_reads\n"
        "sample_a\t0.10\t0.90\t0.00\t100\t0\n",
        encoding="utf-8",
    )

    rc = _run_cli(["assay", "check", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 0, rc.stderr
    reliability = json.loads((out_dir / "reliability_summary.json").read_text(encoding="utf-8"))
    finding_ids = {finding["finding_id"] for finding in reliability["findings"]}
    assert "assignment_rate_below_min" not in finding_ids
    assert "ambiguous_rate_above_max" not in finding_ids


def test_assay_start_check_only_matches_assay_check(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)
    spec.write_text(spec.read_text(encoding="utf-8").replace("k = 1", "k = 0"), encoding="utf-8")

    check_rc = _run_cli(["assay", "check", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})
    start_rc = _run_cli(
        ["assay", "start", "--check-only", str(spec)],
        env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")},
    )

    assert start_rc.returncode == check_rc.returncode == 0
    assert "check passed" in start_rc.stderr
    assert ": running" not in start_rc.stderr
    assert not (tmp_path / "assay_out" / "counts.mageck.tsv").exists()


def test_assay_start_prints_preflight_verdict_before_continuing_advisory_run(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)

    rc = _run_cli(["assay", "start", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 2
    stderr = rc.stderr
    assert "preflight failed (unsafe_targets); continuing run" in stderr
    assert "reliability (preflight): failed" in stderr
    assert "preflight checks failed" in stderr
    assert "finding: unsafe_targets:" in stderr
    assert ": running" in stderr
    assert "reliability (postrun):" in stderr
    preflight_idx = stderr.index("reliability (preflight):")
    running_idx = stderr.index(": running")
    postrun_idx = stderr.index("reliability (postrun):")
    assert preflight_idx < running_idx < postrun_idx


def test_assay_check_suppresses_audit_warning_when_unsafe_targets_finding_present(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)

    rc = _run_cli(["assay", "check", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 2
    assert "check failed (unsafe_targets)" in rc.stderr
    assert "dotmatch assay: warning:" not in rc.stderr


def test_assay_run_count_reproduces_existing_crispr_fixture(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)

    rc = _run_cli(["assay", "run", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 0, rc.stderr
    out_dir = tmp_path / "assay_out"
    assert (out_dir / "counts.mageck.tsv").read_text(encoding="utf-8") == (
        ROOT / "examples/workflows/fixtures/expected_counts.mageck.tsv"
    ).read_text(encoding="utf-8")
    manifest = json.loads((out_dir / "assay_manifest.json").read_text(encoding="utf-8"))
    assert manifest["mode"] == "count"
    assert manifest["commands"][0]["name"] == "audit"
    assert manifest["commands"][-1]["name"] == "validate"
    assert manifest["commands"][-2]["name"] == "crispr-qc"
    assert manifest["artifacts"]["reliability_summary"].endswith("reliability_summary.json")
    assert manifest["artifacts"]["reliability_findings"].endswith("reliability_findings.tsv")
    assert manifest["artifacts"]["reliability_report"].endswith("reliability_report.html")
    assert manifest["artifacts"]["reliability_manifest_summary"].endswith("reliability_manifest.summary.tsv")
    assert manifest["artifacts"]["methods"].endswith("methods.md")
    assert manifest["artifacts"]["citation_bib"].endswith("CITATION.bib")
    assert manifest["artifacts"]["software_versions"].endswith("software_versions.yml")
    assert (out_dir / "crispr_qc.json").exists()
    assert (out_dir / "crispr_qc.summary.tsv").exists()
    assert (out_dir / "crispr_qc.html").exists()
    methods = (out_dir / "methods.md").read_text(encoding="utf-8")
    assert "DotMatch Methods and Citation" in methods
    assert "Edit radius (`k`): `1`" in methods
    assert "ambiguous reads were not silently counted" in methods
    citation = (out_dir / "CITATION.bib").read_text(encoding="utf-8")
    assert "@software{dotmatch" in citation
    assert "doi = {10.5281/zenodo.22167503}" in citation
    versions = (out_dir / "software_versions.yml").read_text(encoding="utf-8")
    assert "dotmatch_python:" in versions
    assert "dotmatch_native:" in versions
    report_html = (out_dir / "assay_report.html").read_text(encoding="utf-8")
    assert "Methods and Citation" in report_html
    assert "methods.md" in report_html
    crispr_qc = json.loads((out_dir / "crispr_qc.json").read_text(encoding="utf-8"))
    assert crispr_qc["assay"] == "crispr_count_qc"
    assert "low_assignment_rate" in {warning["code"] for warning in crispr_qc["warnings"]}
    assert (out_dir / "assay_report.html").exists()
    summary_lines = (out_dir / "assay_manifest.summary.tsv").read_text(encoding="utf-8").splitlines()
    assert summary_lines[0].split("\t") == [
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
    assert summary_lines[1].split("\t")[1:4] == ["count", "crispr", "ready"]
    reliability = json.loads((out_dir / "reliability_summary.json").read_text(encoding="utf-8"))
    assert reliability["stage"] == "postrun"
    assert reliability["overall_status"] == "failed"
    assert reliability["backend"]["authority"] == "cpu"
    assert reliability["backend"]["gpu_status"] == "eligible_but_not_used"
    assert reliability["evidence_boundary"]["status"] == "supported"
    finding_ids = {finding["finding_id"] for finding in reliability["findings"]}
    assert "assignment_rate_below_min" in finding_ids
    assert "ambiguous_rate_above_max" in finding_ids
    assert "unmatched_rate_above_max" in finding_ids
    assert "coverage_fraction_below_min" in finding_ids
    assert "zero_count_fraction_above_max" in finding_ids
    assert "gini_index_above_max" in finding_ids
    assert "top_1pct_fraction_above_max" in finding_ids
    assert "crispr_qc_guide_collision" in finding_ids
    assert reliability["thresholds"]["min_coverage_fraction"] == 0.90
    assert reliability["thresholds"]["max_gini_index"] == 0.50
    assert (out_dir / "reliability_findings.tsv").exists()
    assert (out_dir / "reliability_report.html").exists()
    assert (out_dir / "reliability_manifest.summary.tsv").exists()
    report = (out_dir / "assay_report.html").read_text(encoding="utf-8")
    assert "<title>DotMatch Assay Report</title>" in report
    assert "Run Status" in report
    assert "Reliability" in report
    assert "Sample QC" in report
    assert "Library Audit" in report
    assert "Native Commands" in report
    assert "assay_manifest.json" in report
    assert "reliability_report.html" in report
    assert "report.html" in report
    assert str(tmp_path) not in report
    assert str(ROOT) not in report
    assert "sample_a.fastq" in report


def test_assay_handoff_writes_review_bundle_with_input_and_output_hashes(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)
    env = {"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")}

    run = _run_cli(["assay", "run", str(spec)], env=env)
    assert run.returncode == 0, run.stderr
    handoff = _run_cli(["assay", "handoff", str(spec)], env=env)

    assert handoff.returncode == 0, handoff.stderr
    bundle = tmp_path / "assay_out" / "handoff"
    record = json.loads((bundle / "handoff_manifest.json").read_text(encoding="utf-8"))
    assert record["bundle_type"] == "dotmatch_assay_handoff"
    assert record["reliability_status"] == "failed"
    assert {item["role"] for item in record["inputs"]} == {
        "targets",
        "sample:sample_a",
        "sample:sample_b",
    }
    assert all(len(item["sha256"]) == 64 for item in record["inputs"])
    assert any(item["role"] == "counts" for item in record["review_files"])
    assert (bundle / "review" / "reliability_report.html").exists()
    assert (bundle / "review" / "counts.mageck.tsv").exists()
    assert "It does not include raw reads" in (bundle / "README_FOR_REVIEW.md").read_text(
        encoding="utf-8"
    )
    assert "review/counts.mageck.tsv" in (bundle / "SHA256SUMS").read_text(encoding="utf-8")


def test_assay_handoff_requires_a_completed_run(tmp_path: Path) -> None:
    spec = _write_count_spec(tmp_path)

    result = _run_cli(["assay", "handoff", str(spec)])

    assert result.returncode == 2
    assert "requires a completed assay run" in result.stderr


def test_assay_handoff_rejects_an_output_path_that_is_a_file(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)
    env = {"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")}
    run = _run_cli(["assay", "run", str(spec)], env=env)
    assert run.returncode == 0, run.stderr
    output_path = tmp_path / "handoff-file"
    output_path.write_text("not a directory\n", encoding="utf-8")

    result = _run_cli(
        ["assay", "handoff", str(spec), "--out-dir", str(output_path)],
        env=env,
    )

    assert result.returncode == 2
    assert "handoff output path must be a directory" in result.stderr


def test_assay_run_manifest_records_configured_reliability_thresholds(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)
    spec.write_text(
        spec.read_text(encoding="utf-8").replace(
            "fail_on_unsafe_targets = false",
            "fail_on_unsafe_targets = false\nmin_assignment_rate = 0.25\nmax_unmatched_rate = 0.75",
        ),
        encoding="utf-8",
    )

    rc = _run_cli(["assay", "run", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 0, rc.stderr
    manifest = json.loads((tmp_path / "assay_out" / "assay_manifest.json").read_text(encoding="utf-8"))
    assert manifest["autopsy_thresholds"]["min_assignment_rate"] == 0.25
    assert manifest["autopsy_thresholds"]["max_unmatched_rate"] == 0.75


def test_assay_run_demux_and_pair_count_specs(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    demux_spec = _write_demux_spec(tmp_path)
    pair_spec = _write_pair_spec(tmp_path)
    paired_pair_spec = _write_paired_pair_spec(tmp_path)
    env = {"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")}

    demux = _run_cli(["assay", "run", str(demux_spec)], env=env)
    pair = _run_cli(["assay", "run", str(pair_spec)], env=env)
    paired_pair = _run_cli(["assay", "run", str(paired_pair_spec)], env=env)

    assert demux.returncode == 0, demux.stderr
    assert (tmp_path / "demux_out" / "demuxed" / "bc0.fastq").exists()
    assert (tmp_path / "demux_out" / "ambiguous.fastq").exists()
    assert pair.returncode == 0, pair.stderr
    assert "L0\tR0\t1" in (tmp_path / "pair_out" / "pair_counts.tsv").read_text(encoding="utf-8")
    assert "L1\tR1\t1" in (tmp_path / "pair_out" / "pair_counts.tsv").read_text(encoding="utf-8")
    pair_reliability = json.loads((tmp_path / "pair_out" / "reliability_summary.json").read_text(encoding="utf-8"))
    assert pair_reliability["evidence_boundary"]["status"] == "smoke"
    assert any(finding["finding_id"] == "evidence_boundary_not_supported" for finding in pair_reliability["findings"])
    assert paired_pair.returncode == 0, paired_pair.stderr
    paired_out = tmp_path / "paired_pair_out"
    assert "L0\tR0\t1" in (paired_out / "pair_counts.tsv").read_text(encoding="utf-8")
    assert "L1\tR1\t1" in (paired_out / "pair_counts.tsv").read_text(encoding="utf-8")
    paired_summary = json.loads((paired_out / "pair_summary.json").read_text(encoding="utf-8"))
    assert paired_summary["input_mode"] == "paired-fastq"
    assert paired_summary["input_sync"] == "canonical-read-id"
    assert paired_summary["total_pairs"] == 2
    assert (paired_out / "pair_assignments.tsv").read_text(encoding="utf-8").splitlines()[1].startswith("p0\t")
    methods = (paired_out / "methods.md").read_text(encoding="utf-8")
    assert "Left FASTQ" in methods
    assert "Right FASTQ" in methods


def test_pair_assayspec_rejects_mixed_fastq_layouts(tmp_path: Path) -> None:
    from dotmatch.assayspec import AssaySpecError, load_assay_spec

    spec = _write_paired_pair_spec(tmp_path)
    spec.write_text(
        spec.read_text(encoding="utf-8").replace(
            f'left_reads = "{tmp_path / "pair_R1.fastq"}"',
            f'reads = "{tmp_path / "pair_R1.fastq"}"\nleft_reads = "{tmp_path / "pair_R1.fastq"}"',
        ),
        encoding="utf-8",
    )

    with pytest.raises(AssaySpecError, match="must use reads or both left_reads and right_reads"):
        load_assay_spec(spec)


def test_demux_gpu_metadata_requires_public_gpu_gate(tmp_path: Path) -> None:
    from dotmatch.assayspec import _backend_summary, load_assay_spec

    assay = load_assay_spec(_write_demux_spec(tmp_path))

    assert _backend_summary(assay)["gpu_status"] == "compute_compatible_no_public_gpu_gate"


def test_backend_optimizer_recommends_gpu_candidate_for_public_crispr(tmp_path: Path) -> None:
    from dotmatch.assayspec import load_assay_spec, optimize_assay_backend

    assay = load_assay_spec(_write_count_spec(tmp_path))
    plan = optimize_assay_backend(assay)

    assert plan["authority"] == "cpu"
    assert plan["selected_backend"] == "cpu"
    assert plan["candidate_backend"] == "gpu-metal-experimental"
    assert plan["recommendation"] == "gpu_candidate_requires_cpu_validation"
    assert plan["expected_speedup_band"] == "1.5-3x"
    assert "public_gpu_gate_validated" in plan["reason_codes"]
    assert "cpu_count_checksum_required" in plan["accuracy_gates"]
    assert plan["cpu_strategy"] == "cpu_hamming_seed_index"
    assert plan["benchmark_prior_count"] == 3
    assert plan["benchmark_confidence"] == "public_prior"
    assert plan["thread_hint"]["recommended_threads"] == 1
    assert "configured_threads_cap" in plan["thread_hint"]["reason_codes"]
    assert "cpu_remains_assignment_authority" in plan["diagnostic_constraints"]
    assert "hamming_seed_index_available" in plan["route_reasons"]
    assert "fixed_length_acgt_targets" in plan["route_reasons"]
    assert "gpu_candidate_public_gate" in plan["route_reasons"]


def test_backend_optimizer_gates_compute_compatible_demux_without_public_gpu_gate(tmp_path: Path) -> None:
    from dotmatch.assayspec import load_assay_spec, optimize_assay_backend

    assay = load_assay_spec(_write_demux_spec(tmp_path))
    plan = optimize_assay_backend(assay)

    assert plan["selected_backend"] == "cpu"
    assert plan["candidate_backend"] == "gpu-metal-experimental"
    assert plan["recommendation"] == "gpu_candidate_gated"
    assert "compute_compatible_no_public_gpu_gate" in plan["reason_codes"]
    assert plan["cpu_strategy"] == "cpu_hamming_seed_index"
    assert plan["benchmark_prior_count"] == 3
    assert plan["benchmark_confidence"] == "nearest_prior"
    assert "hamming_seed_index_available" in plan["route_reasons"]
    assert "gpu_candidate_without_public_gate" in plan["route_reasons"]
    assert plan["thread_hint"]["recommended_threads"] >= 1
    assert "small_target_set" in plan["thread_hint"]["reason_codes"]


def test_backend_optimizer_requires_cpu_for_levenshtein(tmp_path: Path) -> None:
    from dotmatch.assayspec import load_assay_spec, optimize_assay_backend

    spec = _write_count_spec(tmp_path)
    spec.write_text(
        spec.read_text(encoding="utf-8").replace('metric = "hamming"', 'metric = "levenshtein"'),
        encoding="utf-8",
    )
    assay = load_assay_spec(spec)
    plan = optimize_assay_backend(assay)

    assert plan["candidate_backend"] == "cpu"
    assert plan["recommendation"] == "cpu_required"
    assert "metric_not_gpu_supported" in plan["reason_codes"]
    assert plan["cpu_strategy"] == "cpu_levenshtein_indexed"
    assert plan["benchmark_confidence"] == "unsupported_route"
    assert "levenshtein_indexed_cpu" in plan["route_reasons"]
    assert "metric_not_gpu_supported" in plan["route_reasons"]
    assert "gpu_candidate_requires_zero_mismatch_diagnostic" in plan["diagnostic_constraints"]


def test_backend_optimizer_routes_non_acgt_hamming_to_cpu_only(tmp_path: Path) -> None:
    from dotmatch.assayspec import load_assay_spec, optimize_assay_backend

    assay = load_assay_spec(_write_non_acgt_count_spec(tmp_path))
    plan = optimize_assay_backend(assay)

    assert plan["selected_backend"] == "cpu"
    assert plan["candidate_backend"] == "cpu"
    assert plan["recommendation"] == "cpu_required"
    assert plan["cpu_strategy"] == "cpu_hamming_indexed"
    assert plan["benchmark_prior_count"] == 3
    assert plan["benchmark_confidence"] == "unsupported_route"
    assert "target_alphabet_not_gpu_packable" in plan["reason_codes"]
    assert "non_acgt_targets_cpu_only" in plan["route_reasons"]
    assert "gpu_ineligible_cpu_only" in plan["route_reasons"]
    assert plan["thread_hint"]["recommended_threads"] <= plan["thread_hint"]["max_threads"]
    assert "small_target_set" in plan["thread_hint"]["reason_codes"]


def test_assay_optimize_writes_backend_optimization_artifact(tmp_path: Path) -> None:
    spec = _write_count_spec(tmp_path)

    rc = _run_cli(["assay", "optimize", str(spec)])

    assert rc.returncode == 0, rc.stderr
    assert "gpu-metal-experimental" in rc.stdout
    artifact = tmp_path / "assay_out" / "backend_optimization.json"
    data = json.loads(artifact.read_text(encoding="utf-8"))
    assert data["authority"] == "cpu"
    assert data["candidate_backend"] == "gpu-metal-experimental"


def test_assay_init_writes_requested_template(tmp_path: Path) -> None:
    spec = tmp_path / "assay.toml"

    rc = _run_cli(["assay", "init", "--template", "inline-barcode-demux", "--out", str(spec)])

    assert rc.returncode == 0, rc.stderr
    text = spec.read_text(encoding="utf-8")
    assert 'mode = "demux"' in text
    assert 'assay_type = "inline_barcode"' in text


def test_non_assay_cli_delegates_to_native_binary(tmp_path: Path) -> None:
    native = tmp_path / "dotmatch-native"
    native.write_text("#!/bin/sh\necho native:$@\n", encoding="utf-8")
    native.chmod(0o755)

    rc = _run_cli(["dist", "ACGT", "AGGT"], env={"DOTMATCH_NATIVE_CLI": str(native)})

    assert rc.returncode == 0
    assert rc.stdout.strip() == "native:dist ACGT AGGT"


def test_assay_new_scaffolds_multi_sample_crispr_project(tmp_path: Path) -> None:
    targets = _write_inference_targets(tmp_path)
    reads_dir = tmp_path / "fastqs"
    reads_dir.mkdir()
    _write_inference_reads(reads_dir, prefix="NN", good=True).rename(reads_dir / "sample_a.fastq")
    _write_inference_reads(reads_dir, prefix="NN", good=True).rename(reads_dir / "sample_b.fastq")
    project = tmp_path / "crispr_screen"

    rc = _run_cli(
        [
            "assay",
            "new",
            "crispr",
            "--library",
            str(targets),
            "--reads-dir",
            str(reads_dir),
            "--out",
            str(project),
        ]
    )

    assert rc.returncode == 0, rc.stderr
    assert (project / "assay.toml").exists()
    assert (project / "inference_report.json").exists()
    assert (project / "samples.generated.tsv").exists()
    assert (project / "README.md").exists()
    assert (project / "run.sh").exists()
    assert (project / "inputs" / "targets.tsv").exists()
    assert (project / "reads" / "sample_a.fastq").is_file()
    assert (project / "reads" / "sample_b.fastq").is_file()
    text = (project / "assay.toml").read_text(encoding="utf-8")
    assert 'status = "ready"' in text
    assert 'targets = "inputs/targets.tsv"' in text
    assert 'id = "sample_a"' in text
    assert 'id = "sample_b"' in text
    assert 'fastq = "reads/sample_a.fastq"' in text
    report = json.loads((project / "inference_report.json").read_text(encoding="utf-8"))
    assert report["template"] == "crispr"
    assert len(report["samples"]) == 2
    check = _run_cli(["assay", "check", str(project / "assay.toml")])
    assert check.returncode == 0, check.stderr


def test_crispr_new_scaffold_matches_assay_new(tmp_path: Path) -> None:
    targets = _write_inference_targets(tmp_path)
    reads_dir = tmp_path / "fastqs"
    reads_dir.mkdir()
    _write_inference_reads(reads_dir).rename(reads_dir / "plasmid.fastq")
    project = tmp_path / "screen"

    rc = _run_cli(
        [
            "crispr",
            "new",
            "--library",
            str(targets),
            "--reads-dir",
            str(reads_dir),
            "--out",
            str(project),
        ]
    )

    assert rc.returncode == 0, rc.stderr
    assert (project / "assay.toml").exists()


def test_crispr_quickstart_creates_self_contained_reviewable_project(tmp_path: Path) -> None:
    targets = _write_inference_targets(tmp_path)
    reads = _write_inference_reads(tmp_path, prefix="NN", good=True)
    source = tmp_path / "sample.fastq"
    reads.rename(source)
    project = tmp_path / "quickstart"

    rc = _run_cli(
        [
            "crispr",
            "quickstart",
            "--library",
            str(targets),
            "--fastq",
            str(source),
            "--out",
            str(project),
            "--accept-inference",
        ]
    )

    assert rc.returncode == 2, rc.stderr
    assert (project / "assay.toml").exists()
    assert (project / "inference_report.json").exists()
    assert (project / "run.sh").exists()
    staged = project / "reads" / "sample.fastq"
    assert staged.exists()
    assert staged.resolve() != source.resolve()
    assert staged.read_bytes() == source.read_bytes()
    assert "Created reviewable CRISPR project" in rc.stdout
    assert (project / "assay_out" / "reliability_report.html").exists()


def test_crispr_quickstart_keeps_review_only_project_in_draft(tmp_path: Path) -> None:
    targets = _write_inference_targets(tmp_path)
    reads = _write_inference_reads(tmp_path, prefix="NN", good=True)
    source = tmp_path / "sample.fastq"
    reads.rename(source)
    project = tmp_path / "quickstart"

    rc = _run_cli(
        [
            "crispr",
            "quickstart",
            "--library",
            str(targets),
            "--fastq",
            str(source),
            "--out",
            str(project),
        ]
    )

    assert rc.returncode == 0, rc.stderr
    assert 'status = "draft"' in (project / "assay.toml").read_text(encoding="utf-8")
    assert "Draft assay.toml" in (project / "run.sh").read_text(encoding="utf-8")
    assert "quickstart requires explicit --accept-inference" in (
        project / "inference_report.json"
    ).read_text(encoding="utf-8")


def test_crispr_quickstart_rejects_non_fastq_gzip_inputs(tmp_path: Path) -> None:
    targets = _write_inference_targets(tmp_path)
    invalid = tmp_path / "not-a-fastq.gz"
    invalid.write_bytes(b"not FASTQ")

    rc = _run_cli(
        [
            "crispr",
            "quickstart",
            "--library",
            str(targets),
            "--fastq",
            str(invalid),
            "--out",
            str(tmp_path / "quickstart"),
        ]
    )

    assert rc.returncode == 2
    assert "do not look like FASTQ files" in rc.stderr


def test_crispr_quickstart_does_not_delete_existing_staging_path(tmp_path: Path) -> None:
    targets = _write_inference_targets(tmp_path)
    reads = _write_inference_reads(tmp_path, prefix="NN", good=True)
    source = tmp_path / "sample.fastq"
    reads.rename(source)
    project = tmp_path / "quickstart"
    staging_parent = tmp_path / ".quickstart.dotmatch-inputs"
    staging_parent.mkdir()
    sentinel = staging_parent / "do-not-delete.txt"
    sentinel.write_text("preserve", encoding="utf-8")

    rc = _run_cli(
        [
            "crispr",
            "quickstart",
            "--library",
            str(targets),
            "--fastq",
            str(source),
            "--out",
            str(project),
        ]
    )

    assert rc.returncode == 2
    assert "refusing to overwrite existing staging directory" in rc.stderr
    assert sentinel.read_text(encoding="utf-8") == "preserve"


def test_assay_new_refuses_non_empty_project_dir(tmp_path: Path) -> None:
    targets = _write_inference_targets(tmp_path)
    reads_dir = tmp_path / "fastqs"
    reads_dir.mkdir()
    _write_inference_reads(reads_dir).rename(reads_dir / "sample.fastq")
    project = tmp_path / "screen"
    project.mkdir()
    (project / "existing.txt").write_text("x", encoding="utf-8")

    rc = _run_cli(
        [
            "assay",
            "new",
            "crispr",
            "--library",
            str(targets),
            "--reads-dir",
            str(reads_dir),
            "--out",
            str(project),
        ]
    )

    assert rc.returncode == 2
    assert "non-empty" in rc.stderr


def test_assay_new_writes_draft_when_inference_is_low_confidence(tmp_path: Path) -> None:
    targets = _write_inference_targets(tmp_path)
    reads_dir = tmp_path / "fastqs"
    reads_dir.mkdir()
    _write_inference_reads(reads_dir, good=False).rename(reads_dir / "sample.fastq")
    project = tmp_path / "draft_screen"

    rc = _run_cli(
        [
            "assay",
            "new",
            "crispr",
            "--library",
            str(targets),
            "--reads-dir",
            str(reads_dir),
            "--out",
            str(project),
        ]
    )

    assert rc.returncode == 0, rc.stderr
    assert 'status = "draft"' in (project / "assay.toml").read_text(encoding="utf-8")
    assert "Promote To Ready" in (project / "README.md").read_text(encoding="utf-8")


def test_detect_pythonpath_for_scaffold_ignores_relative_env_entries(monkeypatch: pytest.MonkeyPatch) -> None:
    import dotmatch.assayspec as assayspec

    monkeypatch.setenv("PYTHONPATH", "python")
    detected = assayspec._detect_pythonpath_for_scaffold()

    assert Path(detected).is_absolute()
    assert detected.endswith("/python")
    assert ":python" not in detected


def test_scaffold_run_script_skips_pythonpath_for_installed_launcher(tmp_path: Path) -> None:
    import dotmatch.assayspec as assayspec

    project = tmp_path / "screen"
    project.mkdir()
    assayspec._write_scaffold_run_script(
        project,
        status="ready",
        launcher=["/usr/local/bin/dotmatch"],
        native_cli=None,
    )

    run_script = (project / "run.sh").read_text(encoding="utf-8")
    assert "DOTMATCH_LAUNCHER=(" in run_script
    scaffold_header = run_script.split("if ! _dotmatch_ready", maxsplit=1)[0]
    assert "PYTHONPATH=" not in scaffold_header


def test_assay_new_run_script_uses_start(tmp_path: Path) -> None:
    targets = _write_inference_targets(tmp_path)
    reads_dir = tmp_path / "fastqs"
    reads_dir.mkdir()
    _write_inference_reads(reads_dir).rename(reads_dir / "sample.fastq")
    project = tmp_path / "screen"

    rc = _run_cli(
        [
            "assay",
            "new",
            "crispr",
            "--library",
            str(targets),
            "--reads-dir",
            str(reads_dir),
            "--out",
            str(project),
        ]
    )

    assert rc.returncode == 0, rc.stderr
    run_script = (project / "run.sh").read_text(encoding="utf-8")
    readme = (project / "README.md").read_text(encoding="utf-8")
    assert "DOTMATCH_LAUNCHER=(" in run_script
    assert 'assay start assay.toml' in run_script
    assert "PYTHONPATH=" in run_script  # pytest invokes dotmatch via python -m
    if 'status = "draft"' in (project / "assay.toml").read_text(encoding="utf-8"):
        assert "Draft assay.toml" in run_script
        assert "assay_out/assay_fixes.tsv" in run_script
    assert "assay_fixes.tsv" in readme
    assert "reliability_report.html" in readme


def test_assay_start_runs_check_and_prints_reliability_verdict(tmp_path: Path) -> None:
    import dotmatch.assayspec as assayspec

    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)
    spec.write_text(spec.read_text(encoding="utf-8").replace("k = 1", "k = 0"), encoding="utf-8")

    rc = _run_cli(["assay", "start", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    out_dir = tmp_path / "assay_out"
    reliability = json.loads((out_dir / "reliability_summary.json").read_text(encoding="utf-8"))
    assert rc.returncode == assayspec._reliability_exit_code(reliability["overall_status"]), rc.stderr
    assert "preflight passed" in rc.stderr
    assert ": running" in rc.stderr
    assert f"reliability (postrun): {reliability['overall_status']}" in rc.stderr
    assert "report: assay_out/reliability_report.html" in rc.stderr
    assert "next:" in rc.stderr
    assert (out_dir / "assay_fixes.tsv").exists()


def test_assay_start_exit_code_matches_reliability_status() -> None:
    import dotmatch.assayspec as assayspec

    assert assayspec._reliability_exit_code("passed") == 0
    assert assayspec._reliability_exit_code("needs_review") == 1
    assert assayspec._reliability_exit_code("failed") == 2
    assert assayspec._reliability_exit_code("blocked") == 2


def test_assay_start_on_draft_prints_reliability_verdict(tmp_path: Path) -> None:
    spec = _write_count_spec(tmp_path)
    spec.write_text('status = "draft"\n' + spec.read_text(encoding="utf-8"), encoding="utf-8")

    rc = _run_cli(["assay", "start", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 2
    assert "preflight blocked (draft_assayspec)" in rc.stderr
    assert "reliability (preflight): blocked" in rc.stderr
    assert "reliability_report.html" in rc.stderr
    assert "assay_fixes.tsv" in rc.stderr
    assert "next:" in rc.stderr
    assert "Promote status" in rc.stderr
    assert "refusing to run draft" not in rc.stderr
    reliability = json.loads((tmp_path / "assay_out" / "reliability_summary.json").read_text(encoding="utf-8"))
    assert reliability["overall_status"] == "blocked"
    assert any(fix["fix_id"] == "promote_status_ready" for fix in reliability["assay_fixes"])


def test_assay_new_pools_reads_from_multiple_fastqs_for_inference(tmp_path: Path) -> None:
    targets = _write_inference_targets(tmp_path)
    reads_dir = tmp_path / "fastqs"
    reads_dir.mkdir()
    _write_inference_reads(reads_dir, prefix="NN", good=True).rename(reads_dir / "sample_a.fastq")
    _write_inference_reads(reads_dir, prefix="NN", good=True).rename(reads_dir / "sample_b.fastq")
    project = tmp_path / "screen"

    rc = _run_cli(
        [
            "assay",
            "new",
            "crispr",
            "--library",
            str(targets),
            "--reads-dir",
            str(reads_dir),
            "--out",
            str(project),
        ]
    )

    assert rc.returncode == 0, rc.stderr
    report = json.loads((project / "inference_report.json").read_text(encoding="utf-8"))
    assert report["inference_read_sources"] == ["reads/sample_a.fastq", "reads/sample_b.fastq"]
    assert report["chosen"]["sampled_reads"] >= 8


def test_assay_fixes_suggest_concrete_toml_edits(tmp_path: Path) -> None:
    import dotmatch.assayspec as assayspec

    spec = _write_count_spec(tmp_path)
    spec.write_text('status = "draft"\n' + spec.read_text(encoding="utf-8"), encoding="utf-8")
    plan = assayspec.compile_assay_plan(assayspec.load_assay_spec(spec))
    reliability = assayspec._build_reliability_summary(plan, stage="preflight")

    fix_ids = {fix["fix_id"] for fix in reliability["assay_fixes"]}
    assert "promote_status_ready" in fix_ids
    draft = next(finding for finding in reliability["findings"] if finding["finding_id"] == "draft_assayspec")
    assert 'status = "ready"' in draft["recommended_action"]


def test_assay_run_writes_assay_fixes_for_unsafe_targets(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_unsafe_count_spec(tmp_path, profile="production")

    rc = _run_cli(["assay", "run", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 2
    out_dir = tmp_path / "unsafe_production_out"
    fixes = (out_dir / "assay_fixes.tsv").read_text(encoding="utf-8")
    assert "assignment_exact_matching" in fixes
    assert "assignment" in fixes
    assert "\tk\t" in fixes or "assignment\tk" in fixes.replace("\t", " ")
    report = (out_dir / "reliability_report.html").read_text(encoding="utf-8")
    assert "Recommended Assay Fixes" in report


def test_assay_infer_writes_ready_crispr_count_spec_and_report(tmp_path: Path) -> None:
    targets = _write_inference_targets(tmp_path)
    reads = _write_inference_reads(tmp_path)
    spec = tmp_path / "inferred.toml"
    report = tmp_path / "inference_report.json"

    rc = _run_cli(
        [
            "assay",
            "infer",
            "--mode",
            "count",
            "--assay-type",
            "crispr",
            "--targets",
            str(targets),
            "--reads",
            str(reads),
            "--sample-id",
            "sample",
            "--out",
            str(spec),
            "--report",
            str(report),
        ]
    )

    assert rc.returncode == 0, rc.stderr
    text = spec.read_text(encoding="utf-8")
    assert 'status = "ready"' in text
    assert 'start = 2' in text
    assert 'length = 4' in text
    data = json.loads(report.read_text(encoding="utf-8"))
    assert data["status"] == "ready"
    assert data["chosen"]["start"] == 2
    assert data["chosen"]["length"] == 4
    assert data["chosen"]["assignment_rate"] == 1.0
    assert (tmp_path / "inference_candidates.tsv").exists()


def test_assay_infer_accepts_gzipped_fastq(tmp_path: Path) -> None:
    from dotmatch.assayspec import infer_assay_spec

    targets = _write_inference_targets(tmp_path)
    reads = _write_inference_reads(tmp_path)
    gz_reads = tmp_path / "infer.fastq.gz"
    with reads.open("rt", encoding="utf-8") as src, gzip.open(gz_reads, "wt", encoding="utf-8") as dst:
        dst.write(src.read())
    spec = tmp_path / "inferred_gz.toml"
    report = tmp_path / "inference_gz_report.json"

    result = infer_assay_spec(
        mode="count",
        assay_type="crispr",
        targets=targets,
        reads=gz_reads,
        sample_id="sample",
        out=spec,
        report=report,
    )

    assert result["spec"] == spec
    assert 'fastq = "{}"'.format(gz_reads) in spec.read_text(encoding="utf-8")
    data = json.loads(report.read_text(encoding="utf-8"))
    assert data["chosen"]["start"] == 2
    assert data["chosen"]["assignment_rate"] == 1.0


def test_assay_infer_rejects_sample_id_toml_injection(tmp_path: Path) -> None:
    targets = _write_inference_targets(tmp_path)
    reads = _write_inference_reads(tmp_path)
    spec = tmp_path / "injected.toml"
    report = tmp_path / "injected_report.json"

    rc = _run_cli(
        [
            "assay",
            "infer",
            "--mode",
            "count",
            "--assay-type",
            "crispr",
            "--targets",
            str(targets),
            "--reads",
            str(reads),
            "--sample-id",
            'x"\n[run]\nout_dir="/tmp/pwn"',
            "--out",
            str(spec),
            "--report",
            str(report),
        ]
    )

    assert rc.returncode == 2
    assert not spec.exists()


def test_assay_infer_escapes_quoted_paths(tmp_path: Path) -> None:
    try:
        import tomllib
    except ModuleNotFoundError:
        import tomli as tomllib

    targets = tmp_path / 'targets"quoted.tsv'
    reads = tmp_path / 'reads"quoted.fastq'
    targets.write_text("guide_a\tACGT\tGENEA\nguide_b\tTTTT\tGENEB\n", encoding="utf-8")
    reads.write_text("@r0\nNNACGTAAAA\n+\nIIIIIIIIII\n@r1\nNNTTTTAAAA\n+\nIIIIIIIIII\n", encoding="utf-8")
    spec = tmp_path / "quoted.toml"
    report = tmp_path / "quoted_report.json"

    rc = _run_cli(
        [
            "assay",
            "infer",
            "--mode",
            "count",
            "--assay-type",
            "crispr",
            "--targets",
            str(targets),
            "--reads",
            str(reads),
            "--out",
            str(spec),
            "--report",
            str(report),
        ]
    )

    assert rc.returncode == 0, rc.stderr
    data = tomllib.loads(spec.read_text(encoding="utf-8"))
    assert data["targets"] == str(targets)
    assert data["samples"][0]["fastq"] == str(reads)


def test_assay_infer_low_confidence_writes_draft_spec(tmp_path: Path) -> None:
    targets = _write_inference_targets(tmp_path)
    reads = _write_inference_reads(tmp_path, good=False)
    spec = tmp_path / "draft.toml"
    report = tmp_path / "draft_report.json"

    rc = _run_cli(
        [
            "assay",
            "infer",
            "--mode",
            "count",
            "--assay-type",
            "crispr",
            "--targets",
            str(targets),
            "--reads",
            str(reads),
            "--sample-id",
            "sample",
            "--out",
            str(spec),
            "--report",
            str(report),
        ]
    )

    assert rc.returncode == 0, rc.stderr
    assert 'status = "draft"' in spec.read_text(encoding="utf-8")
    data = json.loads(report.read_text(encoding="utf-8"))
    assert data["status"] == "draft"
    assert data["warnings"]


def test_assay_run_refuses_draft_specs(tmp_path: Path) -> None:
    spec = _write_count_spec(tmp_path)
    spec.write_text('status = "draft"\n' + spec.read_text(encoding="utf-8"), encoding="utf-8")

    rc = _run_cli(["assay", "run", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 2
    assert "draft" in rc.stderr
    reliability = json.loads((tmp_path / "assay_out" / "reliability_summary.json").read_text(encoding="utf-8"))
    assert reliability["overall_status"] == "blocked"
    assert any(finding["finding_id"] == "draft_assayspec" for finding in reliability["findings"])


def test_assay_run_allows_draft_when_reliability_policy_allows_it(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)
    spec.write_text(
        'status = "draft"\n'
        + spec.read_text(encoding="utf-8").replace("fail_on_unsafe_targets = false", "fail_on_unsafe_targets = false\nfail_on_draft_inference = false"),
        encoding="utf-8",
    )

    rc = _run_cli(["assay", "run", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 0, rc.stderr
    reliability = json.loads((tmp_path / "assay_out" / "reliability_summary.json").read_text(encoding="utf-8"))
    draft_findings = [finding for finding in reliability["findings"] if finding["finding_id"] == "draft_assayspec"]
    assert draft_findings and draft_findings[0]["severity"] == "error"


def test_assay_run_records_exploratory_draft_without_blocking(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)
    spec.write_text(
        'status = "draft"\n'
        + spec.read_text(encoding="utf-8").replace("fail_on_unsafe_targets = false", 'profile = "exploratory"\nfail_on_unsafe_targets = false'),
        encoding="utf-8",
    )

    rc = _run_cli(["assay", "run", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 0, rc.stderr
    reliability = json.loads((tmp_path / "assay_out" / "reliability_summary.json").read_text(encoding="utf-8"))
    draft_findings = [finding for finding in reliability["findings"] if finding["finding_id"] == "draft_assayspec"]
    assert draft_findings and draft_findings[0]["severity"] == "warning"


def test_assay_run_writes_reliability_artifacts_when_native_command_fails(tmp_path: Path) -> None:
    native = tmp_path / "dotmatch-native"
    native.write_text(
        "#!/bin/sh\n"
        "if [ \"$1\" = '--version' ]; then echo 'dotmatch 0.test'; exit 0; fi\n"
        "if [ \"$1\" = 'audit' ]; then mkdir -p \"$6\"; printf '{\"safe_at_k1\": true}\\n' > \"$6/audit_summary.json\"; printf 'metric\\tvalue\\n' > \"$6/audit_summary.tsv\"; exit 0; fi\n"
        "echo forced failure >&2\n"
        "exit 9\n",
        encoding="utf-8",
    )
    native.chmod(0o755)
    spec = _write_count_spec(tmp_path)

    rc = _run_cli(["assay", "run", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(native)})

    assert rc.returncode == 9
    out_dir = tmp_path / "assay_out"
    reliability = json.loads((out_dir / "reliability_summary.json").read_text(encoding="utf-8"))
    assert reliability["overall_status"] == "failed"
    assert any(finding["finding_id"] == "command_failed" and finding["observed"] == "9" for finding in reliability["findings"])
    assert (out_dir / "reliability_report.html").exists()


def test_assay_run_writes_reliability_artifacts_when_native_cli_missing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import dotmatch.assayspec as assayspec

    spec = _write_count_spec(tmp_path)
    plan = assayspec.compile_assay_plan(assayspec.load_assay_spec(spec))

    def missing_native() -> Path:
        raise FileNotFoundError("missing native")

    monkeypatch.setattr(assayspec, "find_native_cli", missing_native)

    rc = assayspec.run_assay_plan(plan)

    assert rc == 2
    reliability = json.loads((tmp_path / "assay_out" / "reliability_summary.json").read_text(encoding="utf-8"))
    assert reliability["overall_status"] == "blocked"
    assert any(finding["finding_id"] == "native_cli_missing" for finding in reliability["findings"])


def test_assay_rejects_parent_traversal_paths(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    outside = tmp_path / "targets.tsv"
    outside.write_text("guide_a\tACGT\tGENEA\n", encoding="utf-8")
    reads = workspace / "reads.fastq"
    reads.write_text("@r0\nACGTAAAA\n+\nIIIIIIII\n", encoding="utf-8")
    spec = workspace / "escape.toml"
    spec.write_text(
        """
schema_version = 1
mode = "count"
assay_type = "crispr"
targets = "../targets.tsv"

[[samples]]
id = "sample"
fastq = "reads.fastq"

[run]
out_dir = "assay_out"

[extract]
start = 0
length = 4

[assignment]
k = 1
metric = "hamming"
""".lstrip(),
        encoding="utf-8",
    )

    rc = _run_cli(["assay", "check", str(spec)])

    assert rc.returncode == 2
    assert "must stay inside" in rc.stderr
    assert not (workspace / "assay_out").exists()


def test_sample_qc_representation_thresholds_use_fixture_metrics() -> None:
    import dotmatch.assayspec as assayspec

    spec = assayspec.load_assay_spec(ROOT / "examples/workflows/fixtures/crispr_assay.toml")
    plan = assayspec.compile_assay_plan(spec)
    sample_qc = ROOT / "examples/workflows/fixtures/expected_sample_qc.tsv"
    plan.artifacts["sample_qc"] = sample_qc

    reliability = assayspec._build_reliability_summary(plan, stage="postrun", manifest={})
    finding_ids = {finding["finding_id"] for finding in reliability["findings"]}
    sample_a = [finding for finding in reliability["findings"] if finding["sample_id"] == "sample_a"]
    sample_b = [finding for finding in reliability["findings"] if finding["sample_id"] == "sample_b"]

    assert "coverage_fraction_below_min" in finding_ids
    assert "zero_count_fraction_above_max" in finding_ids
    assert {finding["finding_id"] for finding in sample_a} >= {
        "assignment_rate_below_min",
        "ambiguous_rate_above_max",
        "unmatched_rate_above_max",
        "coverage_fraction_below_min",
        "zero_count_fraction_above_max",
    }
    assert {finding["finding_id"] for finding in sample_b} >= {
        "unmatched_rate_above_max",
        "coverage_fraction_below_min",
        "zero_count_fraction_above_max",
        "gini_index_above_max",
        "top_1pct_fraction_above_max",
    }


def test_malformed_sample_qc_is_error_not_zeroed(tmp_path: Path) -> None:
    import dotmatch.assayspec as assayspec

    spec = _write_count_spec(tmp_path)
    plan = assayspec.compile_assay_plan(assayspec.load_assay_spec(spec))
    out_dir = tmp_path / "assay_out"
    out_dir.mkdir()
    plan.artifacts["sample_qc"].write_text(
        "sample_id\tassignment_rate\tambiguous_rate\tno_match_rate\ttotal_reads\tinvalid_reads\n"
        "sample\tbad\t0\t0\t10\t0\n",
        encoding="utf-8",
    )

    reliability = assayspec._build_reliability_summary(plan, stage="postrun", manifest={})

    assert reliability["overall_status"] == "failed"
    assert any(finding["finding_id"] == "sample_qc_malformed" for finding in reliability["findings"])


def test_assay_autopsy_reports_wrong_offset_findings(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_wrong_offset_spec(tmp_path)
    out_dir = tmp_path / "autopsy"

    rc = _run_cli(["assay", "autopsy", str(spec), "--out-dir", str(out_dir)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 0, rc.stderr
    summary = json.loads((out_dir / "autopsy_summary.json").read_text(encoding="utf-8"))
    findings = (out_dir / "findings.tsv").read_text(encoding="utf-8")
    assert summary["findings_count"] >= 1
    assert "wrong_offset" in findings
    assert (out_dir / "top_unmatched.shifted.tsv").exists()


def test_reliability_verdict_uses_project_relative_paths(tmp_path: Path) -> None:
    import dotmatch.assayspec as assayspec

    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)
    spec.write_text(spec.read_text(encoding="utf-8").replace("k = 1", "k = 0"), encoding="utf-8")

    rc = _run_cli(["assay", "check", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 0, rc.stderr
    artifact_lines = [line for line in rc.stderr.splitlines() if line.startswith(("report:", "fixes:", "next:"))]
    assert "report: assay_out/reliability_report.html" in artifact_lines
    assert "fixes: assay_out/assay_fixes.tsv" in artifact_lines
    assert all(str(tmp_path) not in line for line in artifact_lines)


def test_assay_start_does_not_leak_autopsy_audit_path(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    fixture = ROOT / "examples/workflows/fixtures/crispr_assay.toml"

    rc = _run_cli(["assay", "start", str(fixture)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 2
    assert "autopsy/audit" not in rc.stdout
    assert not any(line.strip().endswith("autopsy/audit") for line in rc.stderr.splitlines())


def test_assay_run_auto_triggers_autopsy_on_bad_qc(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_wrong_offset_spec(tmp_path)

    rc = _run_cli(["assay", "run", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 0, rc.stderr
    manifest = json.loads((tmp_path / "wrong_offset_out" / "assay_manifest.json").read_text(encoding="utf-8"))
    assert manifest["autopsy_triggered"] is True
    assert "autopsy" in manifest["autopsy_artifacts"]
    assert (tmp_path / "wrong_offset_out" / "autopsy" / "findings.tsv").exists()
    report = (tmp_path / "wrong_offset_out" / "assay_report.html").read_text(encoding="utf-8")
    assert "Autopsy" in report
    assert "wrong_offset" in report
    reliability = json.loads((tmp_path / "wrong_offset_out" / "reliability_summary.json").read_text(encoding="utf-8"))
    assert any(finding["finding_id"] == "autopsy_wrong_offset" for finding in reliability["findings"])
    assert (tmp_path / "wrong_offset_out" / "assay_fixes.tsv").exists()
    summary = (tmp_path / "wrong_offset_out" / "assay_manifest.summary.tsv").read_text(encoding="utf-8")
    assert "\ttrue\t" in summary


def test_assay_check_reports_failed_preflight_for_unsafe_targets_without_block(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)

    rc = _run_cli(["assay", "check", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 2
    assert "check failed (unsafe_targets)" in rc.stderr
    assert "finding: unsafe_targets:" in rc.stderr


def test_assay_check_blocks_unsafe_targets_with_actionable_verdict(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_unsafe_count_spec(tmp_path, profile="production")

    rc = _run_cli(["assay", "check", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 2
    assert "check blocked (unsafe_targets)" in rc.stderr
    assert "finding: unsafe_targets:" in rc.stderr
    assert "next:" in rc.stderr
    assert "continuing with explicit ambiguity handling" not in rc.stderr


def test_assay_start_blocks_unsafe_targets_before_run(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_unsafe_count_spec(tmp_path, profile="production")

    rc = _run_cli(["assay", "start", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 2
    assert "preflight blocked (unsafe_targets)" in rc.stderr
    assert ": running" not in rc.stderr
    assert "continuing run" not in rc.stderr
    assert "continuing with explicit ambiguity handling" not in rc.stderr
    assert "library safety in the reliability report" not in rc.stderr
    assert "next:" in rc.stderr
    out_dir = tmp_path / "unsafe_production_out"
    assert not (out_dir / "counts.mageck.tsv").exists()
    reliability = json.loads((out_dir / "reliability_summary.json").read_text(encoding="utf-8"))
    assert reliability["overall_status"] == "blocked"
    assert reliability["stage"] == "preflight"
    assert any(finding["finding_id"] == "unsafe_targets" and finding["severity"] == "blocked" for finding in reliability["findings"])


def test_scaffold_run_script_falls_back_to_checkout_pythonpath(tmp_path: Path) -> None:
    import dotmatch.assayspec as assayspec

    project = tmp_path / "screen"
    project.mkdir()
    assayspec._write_scaffold_run_script(
        project,
        status="ready",
        launcher=["/nonexistent/dotmatch"],
        native_cli=None,
    )

    run_script = (project / "run.sh").read_text(encoding="utf-8")
    assert 'for pyroot in "${ROOT}/../../python"' in run_script
    assert 'export PYTHONPATH="${pyroot}' in run_script


def test_scaffold_run_script_unsets_invalid_native_cli(tmp_path: Path) -> None:
    import dotmatch.assayspec as assayspec

    project = tmp_path / "screen"
    project.mkdir()
    assayspec._write_scaffold_run_script(
        project,
        status="ready",
        launcher=["dotmatch"],
        native_cli="/nonexistent/native",
    )

    run_script = (project / "run.sh").read_text(encoding="utf-8")
    assert "unset DOTMATCH_NATIVE_CLI" in run_script


def test_assay_run_production_blocks_unsafe_targets_before_assignment(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_unsafe_count_spec(tmp_path, profile="production")

    rc = _run_cli(["assay", "run", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 2
    assert "unsafe target" in rc.stderr.lower()
    out_dir = tmp_path / "unsafe_production_out"
    assert not (out_dir / "counts.mageck.tsv").exists()
    manifest = json.loads((out_dir / "assay_manifest.json").read_text(encoding="utf-8"))
    assert [command["name"] for command in manifest["commands"]] == ["audit"]
    reliability = json.loads((out_dir / "reliability_summary.json").read_text(encoding="utf-8"))
    assert reliability["overall_status"] == "blocked"
    assert any(finding["finding_id"] == "unsafe_targets" and finding["severity"] == "blocked" for finding in reliability["findings"])


def test_assay_run_exploratory_records_unsafe_targets_without_preflight_block(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_unsafe_count_spec(tmp_path, profile="exploratory")

    rc = _run_cli(["assay", "run", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 0, rc.stderr
    out_dir = tmp_path / "unsafe_exploratory_out"
    assert (out_dir / "counts.mageck.tsv").exists()
    reliability = json.loads((out_dir / "reliability_summary.json").read_text(encoding="utf-8"))
    assert reliability["profile"] == "exploratory"
    assert any(finding["finding_id"] == "unsafe_targets" and finding["severity"] == "warning" for finding in reliability["findings"])


def test_assay_rejects_sample_id_html_injection(tmp_path: Path) -> None:
    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True)
    spec = _write_count_spec(tmp_path)
    spec.write_text(
        spec.read_text(encoding="utf-8").replace('id = "sample_a"', 'id = "sample_<script>alert(1)</script>"'),
        encoding="utf-8",
    )

    rc = _run_cli(["assay", "run", str(spec)], env={"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")})

    assert rc.returncode == 2
    assert "samples[0].id" in rc.stderr


def test_assay_infer_demux_and_pair_reports_are_deterministic(tmp_path: Path) -> None:
    targets = _write_inference_targets(tmp_path)
    reads = _write_inference_reads(tmp_path)
    demux_spec = tmp_path / "demux_inferred.toml"
    demux_report = tmp_path / "demux_report.json"
    pair_spec = tmp_path / "pair_inferred.toml"
    pair_report = tmp_path / "pair_report.json"

    demux = _run_cli(
        [
            "assay",
            "infer",
            "--mode",
            "demux",
            "--assay-type",
            "inline_barcode",
            "--barcodes",
            str(targets),
            "--reads",
            str(reads),
            "--out",
            str(demux_spec),
            "--report",
            str(demux_report),
        ]
    )
    pair = _run_cli(
        [
            "assay",
            "infer",
            "--mode",
            "pair-count",
            "--assay-type",
            "generic",
            "--left-targets",
            str(targets),
            "--right-targets",
            str(targets),
            "--reads",
            str(reads),
            "--out",
            str(pair_spec),
            "--report",
            str(pair_report),
        ]
    )

    assert demux.returncode == 0, demux.stderr
    assert pair.returncode == 0, pair.stderr
    assert json.loads(demux_report.read_text(encoding="utf-8"))["chosen"]["start"] == 2
    pair_data = json.loads(pair_report.read_text(encoding="utf-8"))
    assert pair_data["left"]["chosen"]["start"] == 2
    assert pair_data["right"]["chosen"]["start"] == 2
