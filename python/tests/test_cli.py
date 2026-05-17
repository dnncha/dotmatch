import csv
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
LEGACY_ENV = {**os.environ, "DOTMATCH_PYTHON_NO_DELEGATE": "1"}


def _pyproject_version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    assert match is not None
    return match.group(1)


def test_cli_reports_package_version():
    rc = subprocess.run(
        [sys.executable, "-m", "dotmatch.cli", "--version"],
        check=False,
        env=LEGACY_ENV,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert rc.returncode == 0, rc.stderr
    assert rc.stdout.strip() == f"dotmatch {_pyproject_version()}"
    assert rc.stderr == ""


def _write_fixture_files(tmp_path: Path):
    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "target_id\ttarget_seq\tgene\n"
        "guide_1\tACGT\tTP53\n"
        "guide_2\tTTTT\tBRCA1\n"
        "guide_3\tGGGG\tMYC\n",
        encoding="utf-8",
    )
    reads = tmp_path / "reads.fastq"
    reads.write_text(
        "@exact\n"
        "ACGT\n"
        "+\n"
        "IIII\n"
        "@sub\n"
        "ACGA\n"
        "+\n"
        "IIII\n"
        "@none\n"
        "CCCC\n"
        "+\n"
        "IIII\n",
        encoding="utf-8",
    )
    return targets, reads


def test_count_writes_counts_assignments_and_summary(tmp_path):
    targets, reads = _write_fixture_files(tmp_path)
    counts = tmp_path / "counts.tsv"
    assignments = tmp_path / "assignments.tsv"
    summary = tmp_path / "summary.json"

    rc = subprocess.run(
        [
            sys.executable,
            "-m",
            "dotmatch.cli",
            "count",
            "--targets",
            str(targets),
            "--reads",
            str(reads),
            "--target-start",
            "0",
            "--target-length",
            "4",
            "--k",
            "1",
            "--out",
            str(counts),
            "--assignments",
            str(assignments),
            "--summary",
            str(summary),
        ],
        check=False,
        env=LEGACY_ENV,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert rc.returncode == 0, rc.stderr
    count_lines = counts.read_text(encoding="utf-8").splitlines()
    assert count_lines[0].startswith("target_id\ttarget_seq\tgene")
    assert "guide_1\tACGT\tTP53\t1\t1\t0\t0\t0\t2\t0" in count_lines

    assignment_text = assignments.read_text(encoding="utf-8")
    assert "exact\tACGT\tguide_1\tACGT\t0\tunique" in assignment_text
    assert "sub\tACGA\tguide_1\tACGT\t1\tunique" in assignment_text

    summary_data = json.loads(summary.read_text(encoding="utf-8"))
    assert summary_data["total_reads"] == 3
    assert summary_data["assigned_unique"] == 2
    assert summary_data["unmatched"] == 1


def test_count_reads_gzipped_fastq(tmp_path):
    targets, reads = _write_fixture_files(tmp_path)
    gz_reads = tmp_path / "reads.fastq.gz"
    with gzip.open(gz_reads, "wt", encoding="utf-8") as fh:
        fh.write(reads.read_text(encoding="utf-8"))
    counts = tmp_path / "counts.tsv"

    rc = subprocess.run(
        [
            sys.executable,
            "-m",
            "dotmatch.cli",
            "count",
            "--targets",
            str(targets),
            "--reads",
            str(gz_reads),
            "--target-length",
            "4",
            "--out",
            str(counts),
        ],
        check=False,
        env=LEGACY_ENV,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert rc.returncode == 0, rc.stderr
    assert "guide_1\tACGT\tTP53\t1\t1" in counts.read_text(encoding="utf-8")


def test_audit_targets_reports_k1_unsafe_pairs(tmp_path):
    targets = tmp_path / "targets.tsv"
    targets.write_text("a\tACGT\nb\tACGA\nc\tTTTT\n", encoding="utf-8")
    pairs = tmp_path / "pairs.tsv"

    rc = subprocess.run(
        [
            sys.executable,
            "-m",
            "dotmatch.cli",
            "audit-targets",
            "--targets",
            str(targets),
            "--k",
            "1",
            "--out",
            str(pairs),
        ],
        check=False,
        env=LEGACY_ENV,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert rc.returncode == 0, rc.stderr
    summary = json.loads(rc.stdout)
    assert summary["unsafe_for_k"] is True
    assert summary["pairs_within_k"] == 1
    assert "a\tACGT\tb\tACGA\t1" in pairs.read_text(encoding="utf-8")


def _write_barcode_fixture(tmp_path: Path):
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("barcode_id\tbarcode_seq\ns1\tACGT\ns2\tTTTT\n", encoding="utf-8")
    reads = tmp_path / "reads.fastq"
    reads.write_text(
        "@r1\nNACGTAAAA\n+\nIIIIIIIII\n"
        "@r2\nNTTTTAAAA\n+\nIIIIIIIII\n"
        "@r3\nNGGGGAAAA\n+\nIIIIIIIII\n",
        encoding="utf-8",
    )
    return barcodes, reads


def test_barcode_infer_reports_best_offset(tmp_path):
    barcodes, reads = _write_barcode_fixture(tmp_path)
    out = tmp_path / "offset_scan.tsv"
    summary = tmp_path / "infer.json"

    rc = subprocess.run(
        [
            sys.executable,
            "-m",
            "dotmatch.cli",
            "barcode",
            "infer",
            "--barcodes",
            str(barcodes),
            "--reads",
            str(reads),
            "--scan-starts",
            "0:3",
            "--barcode-length",
            "4",
            "--sample-reads",
            "100",
            "--out",
            str(out),
            "--summary",
            str(summary),
        ],
        check=False,
        env=LEGACY_ENV,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert rc.returncode == 0, rc.stderr
    rows = list(csv.DictReader(out.open(encoding="utf-8"), delimiter="\t"))
    assert rows[0]["start"] == "1"
    assert rows[0]["assignment_rate"] == "0.66666667"
    data = json.loads(summary.read_text(encoding="utf-8"))
    assert data["recommended_start"] == 1
    assert data["warnings"] == []


def test_barcode_autopsy_writes_report_and_provenance(tmp_path):
    barcodes, reads = _write_barcode_fixture(tmp_path)
    fake_native = tmp_path / "dotmatch-native"
    fake_native.write_text(
        "#!/bin/sh\n"
        "set -eu\n"
        "cmd=$1; shift\n"
        "if [ \"$cmd\" = audit ]; then\n"
        "  out=''\n"
        "  while [ $# -gt 0 ]; do [ \"$1\" = --out-dir ] && { out=$2; shift 2; } || shift; done\n"
        "  mkdir -p \"$out\"\n"
        "  printf 'metric\\tvalue\\nsafe_at_k1\\tno\\nrisk_pairs_for_k1\\t1\\n' > \"$out/audit_summary.tsv\"\n"
        "  printf '{\"safe_at_k1\": false, \"risk_pairs_for_k1\": 1}\\n' > \"$out/audit_summary.json\"\n"
        "  printf 'target_a\\ttarget_b\\tsequence_a\\tsequence_b\\tdistance\\n' > \"$out/collision_pairs.tsv\"\n"
        "  printf 'target_id\\tsequence\\tsafe_at_k1\\n' > \"$out/target_safety.tsv\"\n"
        "elif [ \"$cmd\" = demux ]; then\n"
        "  out=''; summary=''; assignments=''; ambiguous=''; unmatched=''\n"
        "  while [ $# -gt 0 ]; do\n"
        "    case \"$1\" in\n"
        "      --out-dir) out=$2; shift 2 ;;\n"
        "      --summary) summary=$2; shift 2 ;;\n"
        "      --assignments) assignments=$2; shift 2 ;;\n"
        "      --ambiguous-out) ambiguous=$2; shift 2 ;;\n"
        "      --unmatched-out) unmatched=$2; shift 2 ;;\n"
        "      *) shift ;;\n"
        "    esac\n"
        "  done\n"
        "  mkdir -p \"$out\"\n"
        "  printf '{\"total_reads\": 3, \"assigned_unique\": 2, \"assigned_exact\": 2, \"assigned_corrected\": 0, \"ambiguous\": 0, \"unmatched\": 1, \"invalid\": 0}\\n' > \"$summary\"\n"
        "  printf 'read_id\\tobserved_barcode\\tstatus\\n' > \"$assignments\"\n"
        "  printf '' > \"$ambiguous\"\n"
        "  printf '@r3\\nNGGGGAAAA\\n+\\nIIIIIIIII\\n' > \"$unmatched\"\n"
        "  printf '@r1\\nNACGTAAAA\\n+\\nIIIIIIIII\\n' > \"$out/s1.fastq\"\n"
        "elif [ \"$cmd\" = inspect-unmatched ]; then\n"
        "  out=''\n"
        "  while [ $# -gt 0 ]; do [ \"$1\" = --out ] && { out=$2; shift 2; } || shift; done\n"
        "  printf 'sequence\\tcount\\toffset_hint\\tlow_quality_count\\tadapter_hint\\treverse_complement_nearest\\tnearest_target\\tnearest_distance\\tdiagnosis\\nGGGG\\t1\\t0\\t0\\t\\t\\ts1\\t4\\tno_match\\n' > \"$out\"\n"
        "fi\n",
        encoding="utf-8",
    )
    fake_native.chmod(0o755)
    out_dir = tmp_path / "autopsy"

    rc = subprocess.run(
        [
            sys.executable,
            "-m",
            "dotmatch.cli",
            "barcode",
            "autopsy",
            "--barcodes",
            str(barcodes),
            "--reads",
            str(reads),
            "--scan-starts",
            "0:3",
            "--barcode-length",
            "4",
            "--k-values",
            "0,1",
            "--out-dir",
            str(out_dir),
        ],
        check=False,
        env={**LEGACY_ENV, "DOTMATCH_NATIVE_CLI": str(fake_native)},
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert rc.returncode == 0, rc.stderr
    assert (out_dir / "report.html").exists()
    assert (out_dir / "report.md").exists()
    assert (out_dir / "offset_scan.tsv").exists()
    assert (out_dir / "collision_graph.tsv").exists()
    assert (out_dir / "correction_safety.tsv").exists()
    assert (out_dir / "findings.tsv").exists()
    assert (out_dir / "provenance.json").exists()
    finding_rows = list(csv.DictReader((out_dir / "findings.tsv").open(encoding="utf-8"), delimiter="\t"))
    assert finding_rows
    assert {"finding", "severity", "evidence", "meaning", "next_action"} <= set(finding_rows[0])
    assert {row["finding"] for row in finding_rows} >= {"high_no_match_rate", "unsafe_one_edit_correction"}
    report = (out_dir / "report.md").read_text(encoding="utf-8")
    assert "Barcode Troubleshooting Report" in report
    assert "Speed is reported only after the comparator settings are documented." in report
    assert "Decision Summary" in report
    assert "QC Checklist" in report
    assert "Exact command provenance is recorded in `provenance.json`." in report
    assert "Findings" in report
    assert "What this means" in report
    assert "Next action" in report
    assert "Do not rescue ambiguous reads into either sample without changing the barcode design or assignment policy." in report


def test_validate_compares_indexed_to_scan(tmp_path):
    targets, reads = _write_fixture_files(tmp_path)
    rc = subprocess.run(
        [
            sys.executable,
            "-m",
            "dotmatch.cli",
            "validate",
            "--targets",
            str(targets),
            "--reads",
            str(reads),
            "--target-length",
            "4",
            "--k",
            "1",
        ],
        check=False,
        env=LEGACY_ENV,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert rc.returncode == 0, rc.stderr
    summary = json.loads(rc.stdout)
    assert summary["oracle"] == "native_scan"
    assert summary["checked_reads"] == 3
    assert summary["mismatches"] == 0


def test_crispr_qc_writes_count_qc_report(tmp_path):
    out_json = tmp_path / "crispr_qc.json"
    out_tsv = tmp_path / "crispr_qc.tsv"
    out_html = tmp_path / "crispr_qc.html"

    rc = subprocess.run(
        [
            sys.executable,
            "-m",
            "dotmatch.cli",
            "crispr-qc",
            "--counts",
            str(ROOT / "examples/workflows/fixtures/expected_counts.mageck.tsv"),
            "--sample-qc",
            str(ROOT / "examples/workflows/fixtures/expected_sample_qc.tsv"),
            "--library",
            str(ROOT / "examples/workflows/fixtures/crispr_library.csv"),
            "--out",
            str(out_json),
            "--summary-tsv",
            str(out_tsv),
            "--report",
            str(out_html),
        ],
        check=False,
        env=LEGACY_ENV,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert rc.returncode == 0, rc.stderr
    report = json.loads(out_json.read_text(encoding="utf-8"))
    assert report["schema_version"] == 1
    assert report["guide_count"] == 3
    assert report["sample_count"] == 2
    assert report["library"]["one_edit_collision_pairs"] == 1
    assert report["library"]["collision_radius_audited"] == 1
    assert report["library"]["safe_for_audited_radius"] is False
    assert report["library"]["safe_for_k"] is False
    assert report["samples"]["sample_a"]["assignment_rate"] == 0.33333333
    assert report["samples"]["sample_a"]["invalid_rate"] == 0.25
    assert report["samples"]["sample_a"]["zero_count_fraction"] == 2 / 3
    assert report["samples"]["sample_a"]["gini_index"] >= 0
    assert report["samples"]["sample_a"]["qc_status"] == "review"
    assert report["sample_correlations"][0]["sample_a"] == "sample_a"
    assert report["sample_correlations"][0]["sample_b"] == "sample_b"
    assert report["replicates"][0]["sample_a"] == "sample_a"
    assert report["replicates"][0]["sample_b"] == "sample_b"
    warning_codes = {warning["code"] for warning in report["warnings"]}
    assert "low_assignment_rate" in warning_codes
    assert "high_invalid_rate" in warning_codes
    assert out_tsv.read_text(encoding="utf-8").splitlines()[0].startswith("sample_id\tqc_status\t")
    html = out_html.read_text(encoding="utf-8")
    assert "<title>DotMatch CRISPR QC Report</title>" in html
    assert "Guide Library Audit" in html


def test_crispr_qc_k2_library_safety_is_not_overclaimed(tmp_path):
    out_json = tmp_path / "crispr_qc_k2.json"

    rc = subprocess.run(
        [
            sys.executable,
            "-m",
            "dotmatch.cli",
            "crispr-qc",
            "--counts",
            str(ROOT / "examples/workflows/fixtures/expected_counts.mageck.tsv"),
            "--library",
            str(ROOT / "examples/workflows/fixtures/crispr_library.csv"),
            "--k",
            "2",
            "--out",
            str(out_json),
        ],
        check=False,
        env=LEGACY_ENV,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert rc.returncode == 0, rc.stderr
    report = json.loads(out_json.read_text(encoding="utf-8"))
    assert report["library"]["collision_radius_audited"] == 1
    assert report["library"]["safe_for_k"] is None
    warning_codes = {warning["code"] for warning in report["warnings"]}
    assert "library_collision_audit_radius" in warning_codes


def test_crispr_qc_without_library_does_not_claim_library_safety(tmp_path):
    out_json = tmp_path / "crispr_qc_no_library.json"

    rc = subprocess.run(
        [
            sys.executable,
            "-m",
            "dotmatch.cli",
            "crispr-qc",
            "--counts",
            str(ROOT / "examples/workflows/fixtures/expected_counts.mageck.tsv"),
            "--out",
            str(out_json),
        ],
        check=False,
        env=LEGACY_ENV,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert rc.returncode == 0, rc.stderr
    report = json.loads(out_json.read_text(encoding="utf-8"))
    assert report["status"] == "review"
    assert report["library"]["provided"] is False
    assert report["library"]["safe_for_audited_radius"] is None
    assert report["library"]["safe_for_k"] is None
    warning_codes = {warning["code"] for warning in report["warnings"]}
    assert "sample_qc_not_provided" in warning_codes
    assert "library_not_provided" in warning_codes


def test_crispr_qc_escapes_report_values(tmp_path):
    counts = tmp_path / "counts.tsv"
    sample_qc = tmp_path / "sample_qc.tsv"
    library = tmp_path / "library.csv"
    out_html = tmp_path / "qc.html"
    counts.write_text("sgRNA\tGene\tbad<script>\nguide_a\tGENEA\t1\n", encoding="utf-8")
    sample_qc.write_text(
        "sample_id\tassignment_rate\tambiguous_rate\tno_match_rate\tinvalid_rate\n"
        "bad<script>\t1\t0\t0\t0\n",
        encoding="utf-8",
    )
    library.write_text("id,gRNA.sequence,Gene\nguide_a,ACGT,GENEA\n", encoding="utf-8")

    rc = subprocess.run(
        [
            sys.executable,
            "-m",
            "dotmatch.cli",
            "crispr-qc",
            "--counts",
            str(counts),
            "--sample-qc",
            str(sample_qc),
            "--library",
            str(library),
            "--report",
            str(out_html),
        ],
        check=False,
        env=LEGACY_ENV,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert rc.returncode == 0, rc.stderr
    html = out_html.read_text(encoding="utf-8")
    assert "bad<script>" not in html
    assert "bad&lt;script&gt;" in html


def test_crispr_namespace_qc_matches_crispr_qc_command(tmp_path):
    out_json = tmp_path / "crispr_qc.json"

    rc = subprocess.run(
        [
            sys.executable,
            "-m",
            "dotmatch.cli",
            "crispr",
            "qc",
            "--counts",
            str(ROOT / "examples/workflows/fixtures/expected_counts.mageck.tsv"),
            "--sample-qc",
            str(ROOT / "examples/workflows/fixtures/expected_sample_qc.tsv"),
            "--library",
            str(ROOT / "examples/workflows/fixtures/crispr_library.csv"),
            "--out",
            str(out_json),
        ],
        check=False,
        env=LEGACY_ENV,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert rc.returncode == 0, rc.stderr
    report = json.loads(out_json.read_text(encoding="utf-8"))
    assert report["assay"] == "crispr_count_qc"
    assert report["library"]["one_edit_collision_pairs"] == 1


def test_crispr_namespace_infer_writes_crispr_assayspec(tmp_path):
    targets = tmp_path / "guides.tsv"
    reads = tmp_path / "reads.fastq"
    spec = tmp_path / "assay.toml"
    report = tmp_path / "inference_report.json"
    targets.write_text("guide_a\tACGT\tGENEA\nguide_b\tTTTT\tGENEB\n", encoding="utf-8")
    reads.write_text(
        "@r0\nNNACGTAAAA\n+\nIIIIIIIIII\n"
        "@r1\nNNTTTTAAAA\n+\nIIIIIIIIII\n",
        encoding="utf-8",
    )

    rc = subprocess.run(
        [
            sys.executable,
            "-m",
            "dotmatch.cli",
            "crispr",
            "infer",
            "--library",
            str(targets),
            "--reads",
            str(reads),
            "--out",
            str(spec),
            "--report",
            str(report),
        ],
        check=False,
        env=LEGACY_ENV,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert rc.returncode == 0, rc.stderr
    text = spec.read_text(encoding="utf-8")
    assert 'mode = "count"' in text
    assert 'assay_type = "crispr"' in text
    assert 'start = 2' in text
    assert json.loads(report.read_text(encoding="utf-8"))["status"] == "ready"
