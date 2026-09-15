#!/usr/bin/env python3
"""Build a portable, native-generated SYNTHETIC reviewer demo; no private inputs.

Only the checked-in nine-read fixture is accepted. The engine writes the real
completed bundle, which is checked against the independently enumerated website
fixture before anything is published. --site regenerates one reserved build
output; --out refuses an existing destination. No package installation or network.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SITE_OUTPUT = ROOT / "public" / "examples" / "assignment-review"
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts"))


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(output: Path, *, build_native: bool = True) -> dict:
    """Publish into a new directory only; preserve existing paths on failure."""
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing existing output: {output}")
    if build_native:
        subprocess.run(["make", "shared"], cwd=ROOT, check=True)
    from dotmatch.sensitivity import run_sensitivity
    from dotmatch.sensitivity_review import build_review_data, render_sensitivity_report
    from generate_assignment_demo import generate

    fixture = ROOT / "examples" / "assignment_sensitivity"
    # The checked-in public fixture is a separate, independently enumerated
    # candidate oracle. Never hand-edit displayed counts to make a demo pass.
    oracle = generate()
    checked = json.loads((ROOT / "public" / "assignment-demo.json").read_text())
    if oracle != checked:
        raise ValueError("Native fixture differs from checked website example")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".review-demo-", dir=output.parent) as temp:
        stage = Path(temp) / "demo"
        stage.mkdir()
        bundle = stage / "bundle"
        run_sensitivity(targets=fixture / "targets.tsv", reads=fixture / "reads.fastq",
                        target_start=0, target_length=20, out_dir=bundle,
                        sample_label="synthetic_nine_read_example", write_read_changes=True)
        summary = json.loads((bundle / "summary.json").read_text())
        for key in ("read_count", "target_count", "changed_reads", "outcomes"):
            if summary[key] != oracle[key]:
                raise ValueError(f"Producer disagrees with independent fixture: {key}")
        for name, digest in oracle["inputs_sha256"].items():
            if _sha(fixture / name) != digest:
                raise ValueError("Fixture input identity changed")
        # Check every producer artifact, not just the two viewer source TSVs.
        for name, info in summary["artifacts"].items():
            path = bundle / name
            if Path(name).name != name or path.is_symlink() or not path.is_file():
                raise ValueError("Unexpected artifact path")
            if info != {"bytes": path.stat().st_size, "sha256": _sha(path)}:
                raise ValueError(f"Producer artifact identity mismatch: {name}")
        model = build_review_data(summary, bundle)
        # A separate labelled display does not modify the manifest-bound report.
        report = render_sensitivity_report(summary, bundle)
        banner = ('<aside role="note" aria-label="Synthetic demonstration" '
                  'style="padding:14px 24px;background:#e6f1ef;border-bottom:1px solid #b4ccc8">'
                  '<strong>Synthetic example — not biological data.</strong> '
                  'Compare Exact with Radius k=1, then open guide_A. To inspect the '
                  'recorded reads, attach bundle/read_changes.tsv from the example download.'
                  '</aside>')
        report = report.replace('<body>', '<body>' + banner, 1)
        (stage / "report.html").write_text(report, encoding="utf-8")
        # HTML reports contain only the viewer projection; raw synthetic inputs
        # are copied explicitly into the example archive, never inferred/uploaded.
        inputs = stage / "inputs"
        inputs.mkdir()
        for name in ("targets.tsv", "reads.fastq"):
            shutil.copyfile(fixture / name, inputs / name)
        (stage / "START-HERE.txt").write_text(
            "DOTMATCH ASSIGNMENT REVIEW — SYNTHETIC EXAMPLE\n\n"
            "Open report.html in a current desktop browser. No installation, server, "
            "account or network is needed. This is an unreleased viewer candidate.\n\n"
            "1. Exact and Radius k=1 both assign three reads, but three guide counts differ.\n"
            "2. Open guide_A: its recorded counts are 1, 0, 1.\n"
            "3. Attach bundle/read_changes.tsv to inspect recorded decisions locally.\n"
            "4. Swap policies or export a figure. More assignments are not evidence of accuracy.\n\n"
            "bundle/report.html is the unchanged producer report; report.html is the separately "
            "labelled, manifest-checked demo. inputs/ contains only the public synthetic fixture.\n"
            "summary.json hashes establish internal consistency, not authentication. No biological "
            "validation or branded Safari acceptance is implied.\n", encoding="utf-8")
        try:
            source = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True)
            dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT, capture_output=True, text=True)
            commit = source.stdout.strip() if source.returncode == 0 else None
            source_dirty = bool(dirty.stdout.strip()) if dirty.returncode == 0 else None
        except OSError:
            commit, source_dirty = None, None
        manifest = {
            "schema_version": "dotmatch.review-demo.v1", "kind": "synthetic_software_example",
            "read_count": summary["read_count"], "target_count": summary["target_count"],
            "changed_reads": summary["changed_reads"],
            "source_commit": commit,
            "tracked_source_dirty": source_dirty,
            "renderer_sha256": model["renderer_sha256"],
            "native_library_sha256": summary["native_library_sha256"],
            "implementation_sha256": summary["implementation_sha256"],
            "files": {str(p.relative_to(stage)): {"sha256": _sha(p), "bytes": p.stat().st_size}
                      for p in sorted(stage.rglob("*")) if p.is_file()},
        }
        (stage / "demo-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        archive = stage / "dotmatch-review-example.zip"
        # Snapshot filenames before opening ZIP: the archive never includes itself.
        files = [p for p in sorted(stage.rglob("*")) if p.is_file()]
        with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED) as z:
            for path in files:
                z.write(path, path.relative_to(stage))
        # Reserve destination exclusively; never replace a prior analysis or file.
        output.mkdir()
        try:
            for path in stage.iterdir():
                shutil.move(str(path), output / path.name)
        except BaseException:
            shutil.rmtree(output)
            raise
    return manifest


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--out", type=Path, help="New example directory; never overwritten")
    choice.add_argument("--site", action="store_true", help="Regenerate the reserved public demo build output")
    args = parser.parse_args(argv)
    try:
        if args.site:
            # This exact ignored directory is generated site output, never an
            # assay path provided by the user. Refuse symlink redirection.
            for parent in (ROOT / "public", SITE_OUTPUT.parent, SITE_OUTPUT):
                if parent.is_symlink():
                    raise ValueError("Refusing a symlink in the public demo path")
            SITE_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(prefix=".review-build-", dir=SITE_OUTPUT.parent) as temp:
                ready = Path(temp) / "ready"
                result = build(ready)
                if SITE_OUTPUT.exists():
                    shutil.rmtree(SITE_OUTPUT)
                ready.rename(SITE_OUTPUT)
        else:
            result = build(args.out)
        print(f"Native synthetic viewer demo: {result['read_count']} reads, "
              f"{result['target_count']} targets, {result['changed_reads']} changed records")
        return 0
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"review demo: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
