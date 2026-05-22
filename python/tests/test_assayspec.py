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
    assert not (tmp_path / "assay_out").exists()


def test_assay_check_writes_preflight_reliability_artifacts(tmp_path: Path) -> None:
    spec = _write_count_spec(tmp_path)

    rc = _run_cli(["assay", "check", str(spec)])

    assert rc.returncode == 0, rc.stderr
    out_dir = tmp_path / "assay_out"
    summary = json.loads((out_dir / "reliability_summary.json").read_text(encoding="utf-8"))
    assert summary["stage"] == "preflight"
    assert summary["overall_status"] == "passed"
    assert summary["profile"] == "production"
    assert summary["backend"]["authority"] == "cpu"
    assert summary["backend"]["gpu_status"] == "eligible_but_not_used"
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
    spec = _write_count_spec(tmp_path)
    stale_audit = tmp_path / "assay_out" / "audit" / "audit_summary.json"
    stale_audit.parent.mkdir(parents=True)
    stale_audit.write_text('{"safe_at_k1": false}\n', encoding="utf-8")

    rc = _run_cli(["assay", "check", str(spec)])

    assert rc.returncode == 0, rc.stderr
    reliability = json.loads((tmp_path / "assay_out" / "reliability_summary.json").read_text(encoding="utf-8"))
    assert not any(finding["finding_id"] == "unsafe_targets" for finding in reliability["findings"])


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
    assert (out_dir / "crispr_qc.json").exists()
    assert (out_dir / "crispr_qc.summary.tsv").exists()
    assert (out_dir / "crispr_qc.html").exists()
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
    env = {"DOTMATCH_NATIVE_CLI": str(ROOT / "dotmatch")}

    demux = _run_cli(["assay", "run", str(demux_spec)], env=env)
    pair = _run_cli(["assay", "run", str(pair_spec)], env=env)

    assert demux.returncode == 0, demux.stderr
    assert (tmp_path / "demux_out" / "demuxed" / "bc0.fastq").exists()
    assert (tmp_path / "demux_out" / "ambiguous.fastq").exists()
    assert pair.returncode == 0, pair.stderr
    assert "L0\tR0\t1" in (tmp_path / "pair_out" / "pair_counts.tsv").read_text(encoding="utf-8")
    assert "L1\tR1\t1" in (tmp_path / "pair_out" / "pair_counts.tsv").read_text(encoding="utf-8")
    pair_reliability = json.loads((tmp_path / "pair_out" / "reliability_summary.json").read_text(encoding="utf-8"))
    assert pair_reliability["evidence_boundary"]["status"] == "smoke"
    assert any(finding["finding_id"] == "evidence_boundary_not_supported" for finding in pair_reliability["findings"])


def test_demux_gpu_metadata_requires_public_gpu_gate(tmp_path: Path) -> None:
    from dotmatch.assayspec import _backend_summary, load_assay_spec

    assay = load_assay_spec(_write_demux_spec(tmp_path))

    assert _backend_summary(assay)["gpu_status"] == "compute_compatible_no_public_gpu_gate"


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
    summary = (tmp_path / "wrong_offset_out" / "assay_manifest.summary.tsv").read_text(encoding="utf-8")
    assert "\ttrue\t" in summary


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
