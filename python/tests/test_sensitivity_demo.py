"""Native-produced public demo, provenance and no-clobber regression checks."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
import zipfile

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import build_sensitivity_demo as demo


def test_native_bundle_and_all_exported_identities(tmp_path):
    output = tmp_path / "example"
    manifest = demo.build(output, build_native=False)
    assert (manifest["read_count"], manifest["target_count"], manifest["changed_reads"]) == (9, 5, 5)
    assert manifest["kind"] == "synthetic_software_example"
    for name, info in manifest["files"].items():
        content = (output / name).read_bytes()
        assert len(content) == info["bytes"]
        assert hashlib.sha256(content).hexdigest() == info["sha256"]
    summary = json.loads((output / "bundle/summary.json").read_text())
    for name, info in summary["artifacts"].items():
        assert hashlib.sha256((output / "bundle" / name).read_bytes()).hexdigest() == info["sha256"]
    displayed = (output / "report.html").read_text()
    original = (output / "bundle/report.html").read_text()
    assert "Synthetic example — not biological data." in displayed
    assert '<aside role="note"' in displayed and '</aside>' in displayed
    assert displayed != original  # never edit the hash-bound original report
    assert '"source_manifest_checked":true' in displayed
    assert '"source_manifest_checked":false' in original
    with zipfile.ZipFile(output / "dotmatch-review-example.zip") as archive:
        names = archive.namelist()
        assert "dotmatch-review-example.zip" not in names
        assert set(names) == set(manifest["files"]) | {"demo-manifest.json"}
        assert all(not n.startswith("/") and ".." not in Path(n).parts for n in names)
        assert archive.read("report.html") == (output / "report.html").read_bytes()
        assert archive.read("inputs/reads.fastq") == (ROOT / "examples/assignment_sensitivity/reads.fastq").read_bytes()


@pytest.mark.parametrize("kind", ["directory", "file", "symlink"])
def test_existing_output_untouched(tmp_path, kind):
    output = tmp_path / "existing"
    if kind == "directory":
        output.mkdir(); (output / "keep.txt").write_text("keep")
    elif kind == "file":
        output.write_text("keep")
    else:
        output.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError):
        demo.build(output, build_native=False)
    if kind == "directory":
        assert (output / "keep.txt").read_text() == "keep"
    elif kind == "file":
        assert output.read_text() == "keep"
    else:
        assert output.is_symlink()


def test_disagreement_with_independent_example_refuses(tmp_path, monkeypatch):
    import generate_assignment_demo
    real = generate_assignment_demo.generate
    def disagree():
        result = real(); result["read_count"] += 1; return result
    monkeypatch.setattr(generate_assignment_demo, "generate", disagree)
    with pytest.raises(ValueError, match="differs from checked"):
        demo.build(tmp_path / "refused", build_native=False)
    assert not (tmp_path / "refused").exists()


def test_corrupt_producer_artifact_refuses_without_partial_output(tmp_path, monkeypatch):
    import dotmatch.sensitivity
    real = dotmatch.sensitivity.run_sensitivity
    def corrupt(**kwargs):
        result = real(**kwargs)
        path = Path(kwargs["out_dir"]) / "exact.counts.tsv"
        path.write_bytes(path.read_bytes() + b"wrong\n")
        return result
    monkeypatch.setattr(dotmatch.sensitivity, "run_sensitivity", corrupt)
    with pytest.raises(ValueError, match="artifact identity mismatch"):
        demo.build(tmp_path / "refused", build_native=False)
    assert not (tmp_path / "refused").exists()
    assert not list(tmp_path.glob(".review-demo-*"))


def test_site_rebuild_failure_preserves_previous_demo(tmp_path, monkeypatch):
    public = tmp_path / "public"; public.mkdir()
    target = public / "examples/assignment-review"; target.mkdir(parents=True)
    (target / "report.html").write_text("previous")
    monkeypatch.setattr(demo, "ROOT", tmp_path)
    monkeypatch.setattr(demo, "SITE_OUTPUT", target)
    def fail(*args, **kwargs):
        raise ValueError("native disagreement")
    monkeypatch.setattr(demo, "build", fail)
    assert demo.main(["--site"]) == 2
    assert (target / "report.html").read_text() == "previous"


def test_site_rejects_redirected_parent(tmp_path, monkeypatch):
    public = tmp_path / "public"; public.mkdir()
    elsewhere = tmp_path / "private"; elsewhere.mkdir()
    (public / "examples").symlink_to(elsewhere, target_is_directory=True)
    monkeypatch.setattr(demo, "ROOT", tmp_path)
    monkeypatch.setattr(demo, "SITE_OUTPUT", public / "examples/assignment-review")
    assert demo.main(["--site"]) == 2
    assert list(elsewhere.iterdir()) == []
