import gzip
import json
import subprocess
import sys
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]


def _pyproject_version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    assert match is not None
    return match.group(1)


def test_cli_reports_package_version():
    rc = subprocess.run(
        [sys.executable, "-m", "dotmatch.cli", "--version"],
        check=False,
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
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert rc.returncode == 0, rc.stderr
    summary = json.loads(rc.stdout)
    assert summary["unsafe_for_k"] is True
    assert summary["pairs_within_k"] == 1
    assert "a\tACGT\tb\tACGA\t1" in pairs.read_text(encoding="utf-8")


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
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert rc.returncode == 0, rc.stderr
    summary = json.loads(rc.stdout)
    assert summary["oracle"] == "native_scan"
    assert summary["checked_reads"] == 3
    assert summary["mismatches"] == 0
