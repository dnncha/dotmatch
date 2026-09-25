"""Exercise the canonical installed/module entrypoint and real first-run engines."""
from __future__ import annotations

import csv
import hashlib
import json
import os
import socket
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

from dotmatch import entrypoint, first_run

ROOT = Path(__file__).resolve().parents[2]


def inputs(tmp_path):
    library = tmp_path / "guides.tsv"
    sequence = "GACTAGCTACGATCGTACGA"
    library.write_text(f"target_id\tsequence\tgene\nguide_a\t{sequence}\tGENE_A\n", encoding="utf-8")
    reads = tmp_path / "sample.fastq"
    reads.write_text("".join(f"@r{i}\n{sequence}\n+\n{'I' * 20}\n" for i in range(20)), encoding="utf-8")
    return library, reads


def quick_args(library, reads, output):
    return ["--library", str(library), "--fastq", str(reads), "--out", str(output), "--max-start", "0", "--max-reads", "20"]


def test_fixture_agrees_with_exhaustive_independent_literal_hamming():
    fixture = json.loads((ROOT / "python/dotmatch/data/first-run.json").read_text())
    assert fixture["synthetic"] is True
    # The package copy must not become a second, drifting scientific fixture.
    for key, filename in [("targets_tsv", "targets.tsv"), ("reads_fastq", "reads.fastq")]:
        assert fixture[key] == (ROOT / "examples/assignment_sensitivity" / filename).read_text()
    targets = list(csv.DictReader(fixture["targets_tsv"].splitlines(), delimiter="\t"))
    sequences = fixture["reads_fastq"].splitlines()[1::4]
    assert len(sequences) == fixture["expected"]["read_count"] == 9
    for mode, k, best in [("exact", 0, False), ("radius_k1", 1, False), ("best_k1", 1, True)]:
        outcomes = Counter({name: 0 for name in ("unique", "ambiguous", "none", "invalid")})
        counts = {row["target_id"]: 0 for row in targets}
        for sequence in sequences:
            if len(sequence) < 20:
                outcomes["invalid"] += 1
                continue
            candidates = [(sum(a != b for a, b in zip(sequence[:20], row["sequence"])), row["target_id"]) for row in targets]
            candidates = [(d, label) for d, label in candidates if d <= k]
            if best and candidates:
                minimum = min(d for d, label in candidates)
                candidates = [(d, label) for d, label in candidates if d == minimum]
            status = "unique" if len(candidates) == 1 else "ambiguous" if candidates else "none"
            outcomes[status] += 1
            if status == "unique":
                counts[candidates[0][1]] += 1
        assert dict(outcomes) == fixture["expected"]["outcomes"][mode]
        assert counts == fixture["expected"]["counts"][mode]


def test_real_demo_offline_with_exact_input_and_output_hashes(tmp_path, monkeypatch):
    def no_network(*args, **kwargs):
        raise AssertionError("first run must not connect to a network")
    monkeypatch.setattr(socket.socket, "connect", no_network)
    output = tmp_path / "demo with spaces"
    result = first_run.run_demo(output)
    assert result["status"] == "expected_results_verified"
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["status"] == "complete" and manifest["synthetic"] is True
    assert "index.html" in manifest["files"]
    for name, metadata in manifest["files"].items():
        content = (output / name).read_bytes()
        assert metadata == {"sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)}
    comparison = json.loads((output / "comparison/report.json").read_text())
    assert comparison["scope"] == "complete"
    sample = comparison["samples"][0]
    assert sample["same_total_different_counts"] is True
    assert sample["changed_guides"] == 3
    assert sample["baseline_compared_total"] == sample["candidate_compared_total"] == 3
    assert (output / "sensitivity/read_changes.tsv").exists()
    assert "SYNTHETIC FIRST RUN" in (output / "index.html").read_text()
    assert str(output) not in (output / "index.html").read_text()


def test_demo_failure_is_not_marked_complete(tmp_path, monkeypatch):
    monkeypatch.setattr(first_run, "run_sensitivity", lambda **kwargs: {"read_count": 0, "outcomes": {}})
    assert first_run.demo_main(["--out-dir", str(tmp_path / "demo")]) == 2
    assert not (tmp_path / "demo/manifest.json").exists()
    assert not (tmp_path / "demo/index.html").exists()


def test_demo_refuses_existing_output_without_modification(tmp_path):
    output = tmp_path / "keep"
    output.mkdir()
    sentinel = output / "important.txt"
    sentinel.write_text("original")
    assert first_run.demo_main(["--out-dir", str(output)]) == 2
    assert sentinel.read_text() == "original"
    assert list(output.iterdir()) == [sentinel]


def test_public_module_demo_and_help(tmp_path):
    env = {**os.environ, "PYTHONPATH": str(ROOT / "python")}
    result = subprocess.run([sys.executable, "-m", "dotmatch", "demo", "--out-dir", str(tmp_path / "run")], cwd=tmp_path, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "Synthetic first run verified" in result.stdout
    assert (tmp_path / "run/manifest.json").exists()
    result = subprocess.run([sys.executable, "-m", "dotmatch", "--help"], env=env, capture_output=True, text=True)
    assert result.returncode == 0
    assert "dotmatch demo" in result.stdout and "dotmatch compare-counts" in result.stdout


@pytest.mark.parametrize("linked", [False, True])
def test_quickstart_original_sources_survive_temporary_staging(tmp_path, linked, capsys):
    library, reads = inputs(tmp_path)
    before = reads.read_bytes()
    project = tmp_path / "project"
    args = quick_args(library, reads, project) + (["--link-reads"] if linked else [])
    assert entrypoint.main(["crispr", "quickstart", *args]) == 0
    staged = project / "reads" / reads.name
    assert staged.is_symlink() is linked
    assert staged.read_bytes() == before == reads.read_bytes()
    assert (staged.resolve() == reads.resolve()) is linked
    assert 'status = "draft"' in (project / "assay.toml").read_text()
    assert not (project / "assay_out").exists()
    report = json.loads((project / "inference_report.json").read_text())
    assert report["samples"][0]["source_fastq"] == str(reads)
    assert report["input_storage"] == ("linked" if linked else "copied")
    assert not list(tmp_path.glob(".project.dotmatch-inputs-*"))
    output = capsys.readouterr().out
    assert "filenames are not biological replicate declarations" in output
    if linked:
        assert "not self-contained" in output


def test_every_glob_must_resolve_before_output_creation(tmp_path, capsys):
    library, reads = inputs(tmp_path)
    project = tmp_path / "project"
    args = quick_args(library, reads, project) + ["--fastq", str(tmp_path / "missing*.fastq.gz")]
    assert first_run.quickstart_main(args) == 2
    assert "matched no files" in capsys.readouterr().err
    assert not project.exists()
    assert not list(tmp_path.glob(".project.dotmatch-inputs-*"))


def test_literal_filename_with_glob_metacharacters(tmp_path):
    library, reads = inputs(tmp_path)
    literal = tmp_path / "sample[one].fastq"
    reads.rename(literal)
    assert first_run.resolve_fastqs([str(literal)]) == [literal]
    assert first_run.quickstart_main(quick_args(library, literal, tmp_path / "project")) == 0


def test_overlapping_patterns_refused_without_double_counting(tmp_path):
    library, reads = inputs(tmp_path)
    with pytest.raises(ValueError, match="duplicate FASTQ input"):
        first_run.resolve_fastqs([str(reads), str(tmp_path / "*.fastq")])


def test_symlink_alias_for_same_input_is_duplicate(tmp_path):
    library, reads = inputs(tmp_path)
    alias = tmp_path / "alias.fastq"
    alias.symlink_to(reads)
    with pytest.raises(ValueError, match="duplicate FASTQ input"):
        first_run.resolve_fastqs([str(reads), str(alias)])


@pytest.mark.parametrize("case_variant", [False, True])
def test_duplicate_basenames_refused_before_creating_project(tmp_path, case_variant):
    library, reads = inputs(tmp_path)
    folder = tmp_path / "second"
    folder.mkdir()
    second = folder / ("SAMPLE.fastq" if case_variant else reads.name)
    second.write_bytes(reads.read_bytes())
    project = tmp_path / "project"
    assert first_run.quickstart_main(quick_args(library, reads, project) + ["--fastq", str(second)]) == 2
    assert not project.exists()


@pytest.mark.parametrize("kind", ["directory", "broken_symlink", "wrong_suffix", "fifo"])
def test_non_regular_or_non_fastq_inputs_rejected(tmp_path, kind):
    library, reads = inputs(tmp_path)
    path = tmp_path / "bad.fastq"
    if kind == "directory":
        path.mkdir()
    elif kind == "broken_symlink":
        path.symlink_to(tmp_path / "missing.fastq")
    elif kind == "wrong_suffix":
        path = tmp_path / "not_reads.gz"
        path.write_bytes(b"not a fastq")
    else:
        if not hasattr(os, "mkfifo"):
            pytest.skip("requires FIFO support")
        os.mkfifo(path)
    assert first_run.quickstart_main(quick_args(library, path, tmp_path / "project")) == 2
    assert not (tmp_path / "project").exists()


@pytest.mark.parametrize("existing", ["empty", "file", "broken_symlink"])
def test_quickstart_no_overwrite_even_empty_or_broken_symlink(tmp_path, existing):
    library, reads = inputs(tmp_path)
    output = tmp_path / "project"
    if existing == "empty":
        output.mkdir()
    elif existing == "file":
        output.write_text("preserve")
    else:
        output.symlink_to(tmp_path / "missing")
    assert first_run.quickstart_main(quick_args(library, reads, output)) == 2
    assert os.path.lexists(output)
    if existing == "file":
        assert output.read_text() == "preserve"


@pytest.mark.parametrize("flag,value", [("--threads", "0"), ("--max-reads", "0"), ("--max-start", "-1")])
def test_invalid_limits_refused_before_creating_output(tmp_path, flag, value):
    library, reads = inputs(tmp_path)
    with pytest.raises(SystemExit) as error:
        first_run.quickstart_main(quick_args(library, reads, tmp_path / "project") + [flag, value])
    assert error.value.code == 2
    assert not (tmp_path / "project").exists()


def test_no_run_overrides_accept_inference(tmp_path):
    library, reads = inputs(tmp_path)
    project = tmp_path / "project"
    assert first_run.quickstart_main(quick_args(library, reads, project) + ["--no-run", "--accept-inference"]) == 0
    assert 'status = "draft"' in (project / "assay.toml").read_text()
    assert not (project / "assay_out").exists()


def test_uncertain_inference_is_not_promoted_automatically(tmp_path, monkeypatch):
    library, reads = inputs(tmp_path)
    original = first_run.scaffold_assay_project
    def uncertain(**kwargs):
        result = original(**kwargs)
        path = kwargs["project_dir"] / "inference_report.json"
        report = json.loads(path.read_text())
        report["status"] = "draft"
        path.write_text(json.dumps(report))
        return result
    def refuse(*args):
        pytest.fail("uncertain inference must not start an analysis")
    monkeypatch.setattr(first_run, "scaffold_assay_project", uncertain)
    monkeypatch.setattr(first_run, "command_assay", refuse)
    assert first_run.quickstart_main(quick_args(library, reads, tmp_path / "project") + ["--accept-inference"]) == 2


def test_ready_inference_delegates_to_existing_assay_engine(tmp_path, monkeypatch):
    library, reads = inputs(tmp_path)
    original = first_run.scaffold_assay_project
    def ready(**kwargs):
        result = original(**kwargs)
        path = kwargs["project_dir"] / "inference_report.json"
        report = json.loads(path.read_text())
        report["status"] = "ready"
        path.write_text(json.dumps(report))
        return result
    calls = []
    monkeypatch.setattr(first_run, "scaffold_assay_project", ready)
    monkeypatch.setattr(first_run, "command_assay", lambda args: calls.append(args) or 0)
    project = tmp_path / "project"
    assert first_run.quickstart_main(quick_args(library, reads, project) + ["--accept-inference"]) == 0
    assert calls == [["start", str(project / "assay.toml")]]


def test_other_commands_keep_the_existing_router(monkeypatch):
    from dotmatch import cli
    calls = []
    monkeypatch.setattr(cli, "main", lambda args: calls.append(args) or 7)
    assert entrypoint.main(["count", "--help"]) == 7
    assert calls == [["count", "--help"]]


def test_compare_counts_alias_uses_existing_implementation(monkeypatch):
    from dotmatch import count_compare
    calls = []
    monkeypatch.setattr(count_compare, "main", lambda args: calls.append(args) or 1)
    assert entrypoint.main(["compare-counts", "--help"]) == 1
    assert calls == [["--help"]]
