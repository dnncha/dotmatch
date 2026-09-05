import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from editwitness import analyze, load_manifest
from editwitness.cli import main
from editwitness.io import InputError, atomic_write, read_json, verify_result
from editwitness.report import render_report

ROOT = Path(__file__).resolve().parents[1]


def cli(*args, stdin=None):
    env = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
    return subprocess.run([sys.executable, "-m", "editwitness", *map(str, args)],
                          input=stdin, text=True, capture_output=True, env=env, timeout=15)


def test_cli_full_analysis_html_replay_and_compact(tmp_path):
    manifest = ROOT / "examples/demo.json"
    result_path, report = tmp_path / "analysis.json", tmp_path / "report.html"
    result = cli("analyze", manifest, "--output", result_path, "--html", report)
    assert result.returncode == 0, result.stderr
    assert result.stdout == result.stderr == ""
    assert verify_result(result_path).conclusion == "ambiguity_demonstrated"
    verified = cli("verify", result_path, "--manifest", manifest)
    assert json.loads(verified.stdout)["replayed"] is True
    assert report.read_text(encoding="utf-8").startswith("<!doctype html>")
    compact = cli("analyze", manifest, "--compact")
    data = json.loads(compact.stdout)
    assert data["kind"] == "editwitness.summary"
    assert "allele_observations" not in data
    assert "sequence" not in data
    fail_on = cli("analyze", manifest, "--compact", "--fail-on-ambiguity")
    assert fail_on.returncode == 4
    assert json.loads(fail_on.stdout)["conclusion"] == "ambiguity_demonstrated"


def test_demo_stdin_validation_and_informational_commands():
    demo = cli("demo")
    assert demo.returncode == 0
    assert json.loads(demo.stdout)["reference"]["synthetic"]
    validated = cli("validate", "-", stdin=demo.stdout)
    assert validated.returncode == 0
    assert json.loads(validated.stdout)["valid"]
    for args in (("doctor",), ("capabilities",), ("schema", "manifest"), ("demo", "--paired-end")):
        result = cli(*args)
        assert result.returncode == 0, result.stderr
        assert isinstance(json.loads(result.stdout), dict)
    assert "editwitness" in cli("--help").stdout
    assert "0.1.0a1" in cli("--version").stdout


def test_witness_and_scan_commands():
    path = ROOT / "examples/demo.json"
    result = cli("witness", path, "--hypothesis", "hidden_primer_deletion", "--include-sequences")
    assert result.returncode == 0
    assert json.loads(result.stdout)["allele_observations"]
    assert cli("witness", path, "--hypothesis", "intended_biallelic").returncode == 2
    scan = cli("scan", path)
    assert scan.returncode == 0
    assert json.loads(scan.stdout)["enumerated_deletions"] > 0


def test_output_protection_and_inplace_refusal(tmp_path):
    path = tmp_path / "file.json"
    path.write_text((ROOT / "examples/demo.json").read_text())
    original = path.read_bytes()
    assert cli("analyze", path, "--output", path, "--force").returncode == 2
    assert path.read_bytes() == original
    output = tmp_path / "output.json"
    output.write_text("keep me")
    assert cli("analyze", path, "--output", output).returncode == 3
    assert output.read_text() == "keep me"
    assert cli("analyze", path, "--output", output, "--force").returncode == 0
    assert cli("analyze", path, "--output", output, "--html", output, "--force").returncode == 2
    assert cli("analyze", path, "--html", "-").returncode == 2
    assert cli("demo", "-o", tmp_path / "missing/out.json").returncode == 3


def test_tampered_checksum_and_changed_manifest_fail(demo, tmp_path):
    path = tmp_path / "result.json"
    result = analyze(demo)
    path.write_text(result.model_dump_json())
    changed = tmp_path / "changed.json"
    data = demo.model_dump(mode="json")
    data["candidates"][0]["cost_units"] += 1
    changed.write_text(json.dumps(data))
    replay = cli("verify", path, "--manifest", changed)
    assert replay.returncode == 5
    data = json.loads(path.read_text())
    data["plan"]["cost_units"] = 999
    path.write_text(json.dumps(data))
    result = cli("verify", path)
    assert result.returncode == 5
    assert json.loads(result.stderr)["code"] == "INTEGRITY_MISMATCH"


@pytest.mark.parametrize("text", ['{"a":1,"a":2}', '{"a":NaN}', '{"a":Infinity}', '{', '"unterminated', '9'*5000])
def test_malformed_and_duplicate_json(text, tmp_path):
    path = tmp_path / "bad.json"
    path.write_text(text)
    with pytest.raises(InputError):
        read_json(path)
    result = cli("validate", path)
    assert result.returncode == 2
    assert result.stdout == ""
    assert json.loads(result.stderr)["kind"] == "editwitness.error"


def test_bounded_and_binary_input(tmp_path):
    path = tmp_path / "input"
    path.write_bytes(b" " * 20)
    with pytest.raises(InputError, match="byte limit"):
        read_json(path, max_bytes=10)
    path.write_bytes(b"\xff\x00")
    with pytest.raises(InputError):
        read_json(path)


def test_missing_unknown_and_strict_errors_are_json(tmp_path):
    assert cli("validate", tmp_path / "missing").returncode == 3
    result = cli("unknown-command")
    assert result.returncode == 2
    assert json.loads(result.stderr)["code"] == "INVALID_INPUT"
    result = cli("validate", "-", stdin='{"reference":{"sequence":"SECRET_INVALID_DNA"}}')
    assert result.returncode == 2
    assert "SECRET_INVALID_DNA" not in result.stderr


def test_atomic_write_never_silently_replaces(tmp_path):
    path = tmp_path / "output.txt"
    atomic_write(path, "original")
    with pytest.raises(FileExistsError):
        atomic_write(path, "replacement")
    assert path.read_text() == "original"
    atomic_write(path, "replacement", force=True)
    assert path.read_text() == "replacement"
    assert not list(tmp_path.glob(".editwitness-*"))


def test_html_is_offline_script_free_and_escapes_input(demo):
    data = demo.model_dump(mode="json")
    data["reference"]["name"] = '<script>alert("x")</script>'
    from editwitness.models import Manifest
    result = analyze(Manifest.model_validate(data))
    html = render_report(result)
    assert "<script>" not in html
    assert "&lt;script&gt;" in html
    assert "Content-Security-Policy" in html
    assert "src=\"http" not in html
    assert "software-tested" in html.lower()
    assert "not empirically validated" in html
    assert result.result_sha256 in html


def test_verify_rejects_unknown_and_nonobject_results(tmp_path):
    path = tmp_path / "bad-result.json"
    for data in ([], {"kind": "editwitness.summary"}):
        path.write_text(json.dumps(data))
        assert cli("verify", path).returncode == 2


def test_scan_result_can_be_replayed(tmp_path):
    path = tmp_path / "scan.json"
    source = ROOT / "examples/demo.json"
    assert cli("scan", source, "-o", path).returncode == 0
    assert cli("verify", path, "--manifest", source).returncode == 0


@pytest.mark.parametrize("encoding", ["utf-16", "utf-32"])
def test_json_requires_utf8_even_if_json_library_accepts_other_encodings(encoding, tmp_path):
    path = tmp_path / "input.json"
    path.write_bytes('{"a":1}'.encode(encoding))
    with pytest.raises(InputError, match="UTF-8"):
        read_json(path)


def test_cli_init_and_verifying_result_cannot_overwrite_input(demo, tmp_path):
    from editwitness.sequence import reverse_complement
    fasta = tmp_path / "locus.fasta"
    fasta.write_text(f">local\n{demo.reference.sequence}\n")
    alternate = next(x for x in "ACGT" if x != demo.reference.sequence[450])
    result = cli("init", "--fasta", fasta, "--left-primer", demo.reference.sequence[200:220],
                 "--right-primer", reverse_complement(demo.reference.sequence[680:700]),
                 "--edit-position", 450, "--alternate", alternate)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["expected_hypothesis"] == "intended_biallelic"
    output = tmp_path / "analysis.json"
    output.write_text(analyze(demo).model_dump_json())
    assert cli("verify", output, "--output", output, "--force").returncode == 2
