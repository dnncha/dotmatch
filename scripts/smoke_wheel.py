"""Install a built wheel outside the source tree and exercise its public CLI.

Dependencies are inherited from the test environment; the package itself is
installed into a temporary venv. This requires no dependency downloads.
"""
from __future__ import annotations

import argparse
import json
import os
import site
import subprocess
import sys
import tempfile
import venv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel", type=Path)
    args = parser.parse_args()
    wheel = args.wheel.resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="editwitness-wheel-") as name:
        root = Path(name)
        environment = root / "venv"
        venv.EnvBuilder(with_pip=True, system_site_packages=True).create(environment)
        python = environment / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
        env = dict(os.environ)
        env.pop("PYTHONPATH", None)
        # A venv created from another venv does not inherit the parent's site-packages.
        # Append those dependency directories after this venv's own packages. Plain
        # paths do not execute parent .pth files; the import-location assertion below
        # prevents accidentally testing a parent/global EditWitness installation.
        purelib = subprocess.run([str(python), "-c", "import sysconfig; print(sysconfig.get_path('purelib'))"],
                                 cwd=root, env=env, text=True, capture_output=True, check=True)
        (Path(purelib.stdout.strip()) / "test_dependencies.pth").write_text(
            "\n".join(site.getsitepackages()) + "\n", encoding="utf-8"
        )
        subprocess.run([str(python), "-m", "pip", "install", "--no-deps", "--force-reinstall", str(wheel)],
                       cwd=root, env=env, check=True)
        def run(*arguments: str) -> subprocess.CompletedProcess[str]:
            return subprocess.run([str(python), "-m", "editwitness", *arguments], cwd=root,
                                  env=env, text=True, capture_output=True, check=True, timeout=60)
        location = subprocess.run([str(python), "-c", "import editwitness; print(editwitness.__file__)"],
                                  cwd=root, env=env, text=True, capture_output=True, check=True)
        if not Path(location.stdout.strip()).resolve().is_relative_to(environment.resolve()):
            raise RuntimeError("smoke test imported a source-tree or global package, not the installed wheel")
        run("demo", "-o", "manifest.json")
        run("analyze", "manifest.json", "-o", "evidence.json", "--html", "report.html")
        replay = json.loads(run("verify", "evidence.json", "--manifest", "manifest.json").stdout)
        if replay.get("replayed") is not True:
            raise RuntimeError("installed wheel failed result replay")
        if not (root / "report.html").read_text(encoding="utf-8").startswith("<!doctype html>"):
            raise RuntimeError("installed wheel did not generate an HTML report")
        witness = json.loads(run("witness", "manifest.json", "--hypothesis", "hidden_primer_deletion",
                                 "--include-sequences").stdout)
        if not witness.get("local_allele_sequences"):
            raise RuntimeError("installed wheel did not provide witness sequence evidence")
        run("compare-models", "manifest.json")
        schema = json.loads(run("schema", "analysis").stdout)
        if schema.get("$id") != "urn:editwitness:schema:analysis:1.1":
            raise RuntimeError("installed wheel has stale result schema")
        print(json.dumps({"installed_wheel": wheel.name, "outside_source_tree": True,
                          "replay": True, "witness_sequences": True, "html": True,
                          "dependency_policy": "preinstalled test-environment dependencies"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
