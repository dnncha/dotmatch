"""Create release-bound execution provenance and a synthetic demonstration report.

Called only after the distribution gates. The CI run is the authoritative record
of completed jobs. This file records identity, not an independent endorsement.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
from pathlib import Path

import pydantic

from editwitness import __version__, analyze, load_manifest
from editwitness.report import render_report
from editwitness.selftest import self_test


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    out = root / "release-support"
    out.mkdir(exist_ok=True)
    result = analyze(load_manifest(root / "examples/demo.json"))
    software_check = self_test()
    if not software_check["passed"]:
        raise ValueError("cannot prepare release evidence after a failed software self-test")
    sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    repo, run = os.environ.get("GITHUB_REPOSITORY"), os.environ.get("GITHUB_RUN_ID")
    evidence = {
        "format": "editwitness.release-evidence.v1", "package_version": __version__,
        "source_commit": sha,
        "source_inventory_sha256": hashlib.sha256((root / "release-files.json").read_bytes()).hexdigest(),
        "ci_run_url": f"https://github.com/{repo}/actions/runs/{run}" if repo and run else None,
        "python": platform.python_version(), "pydantic": pydantic.__version__,
        "software_self_test": software_check,
        "empirical_biological_validation": False, "pypi_publication": False,
        "caveat": "Execution provenance, not independent scientific review, assay calibration or a signature.",
    }
    (out / "release-evidence.json").write_text(json.dumps(evidence, indent=2)+"\n", encoding="utf-8")
    (out / "editwitness-demo.html").write_text(render_report(result), encoding="utf-8")
    (out / "EVIDENCE_SHA256SUMS").write_text("".join(
        f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n"
        for p in sorted(out.iterdir()) if p.name != "EVIDENCE_SHA256SUMS"
    ), encoding="ascii")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
