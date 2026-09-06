"""Install and exercise a built wheel outside the checkout.

Default: a clean virtual environment with dependencies installed normally.
--offline: inherit installed dependencies; do not claim dependency-isolated testing.
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
import zipfile
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    wheels = list((root / "dist").glob("*.whl"))
    if len(wheels) != 1:
        parser.error("build exactly one wheel in a clean dist/ first")
    wheel = wheels[0]
    with zipfile.ZipFile(wheel) as archive:
        names = set(archive.namelist())
        required = {"editwitness/data/demo.json", "editwitness/py.typed",
                    "editwitness/schemas/manifest.schema.json", "editwitness/exact.py"}
        if not required <= names:
            raise ValueError(f"wheel lacks required package files: {sorted(required - names)}")
    with tempfile.TemporaryDirectory(prefix="editwitness-wheel-") as directory:
        work = Path(directory)
        env_path = work / "env"
        venv.EnvBuilder(with_pip=True, system_site_packages=args.offline).create(env_path)
        python = env_path / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
        if args.offline:
            # A nested venv inherits the base interpreter's site packages, not
            # necessarily the invoking venv's dependencies. Add those explicitly
            # for this declared dependency-reuse smoke mode only.
            env_site = env_path / ("Lib/site-packages" if os.name == "nt" else
                                   f"lib/python{sys.version_info.major}.{sys.version_info.minor}/site-packages")
            (env_site / "editwitness-test-dependencies.pth").write_text(
                "\n".join(site.getsitepackages()) + "\n", encoding="utf-8"
            )
        env = dict(os.environ, PYTHONNOUSERSITE="1", PIP_DISABLE_PIP_VERSION_CHECK="1")
        env.pop("PYTHONPATH", None)

        def run(*command: str) -> str:
            response = subprocess.run([str(python), *command], cwd=work, env=env,
                text=True, encoding="utf-8", capture_output=True, timeout=180)
            if response.returncode:
                raise RuntimeError(f"Wheel smoke failed: {command!r}\n{response.stderr}\n{response.stdout}")
            return response.stdout

        install = ["-m", "pip", "install", str(wheel.resolve())]
        if args.offline:
            install += ["--no-deps", "--no-index", "--ignore-installed"]
        run(*install)
        location = Path(run("-c", "import editwitness; print(editwitness.__file__)").strip()).resolve()
        if not location.is_relative_to(env_path.resolve()):
            raise ValueError("smoke test imported a source or global package rather than the installed wheel")
        installed_test = json.loads(run("-m", "editwitness", "self-test"))
        if not installed_test["passed"]:
            raise ValueError("installed package software self-test failed")
        for label, flags in (("full", ()), ("paired", ("--paired-end",))):
            run("-m", "editwitness", "demo", *flags, "-o", f"{label}.json")
            run("-m", "editwitness", "analyze", f"{label}.json", "-o", f"{label}-analysis.json",
                "--html", f"{label}.html")
            integrity = json.loads(run("-m", "editwitness", "verify", f"{label}-analysis.json", "--manifest", f"{label}.json"))
            if not integrity["replayed"]:
                raise ValueError("wheel replay did not run")
        run("-m", "editwitness", "compare-models", "full.json")
        demo = json.loads((work / "full.json").read_text(encoding="utf-8"))
        demo["deletion_scan"] = dict(start_min=195, start_max=201, end_min=215, end_max=221, step=3)
        (work / "grid.json").write_text(json.dumps(demo), encoding="utf-8")
        run("-m", "editwitness", "expand-deletions", "grid.json", "-o", "expanded.json")
        expanded = json.loads(run("-m", "editwitness", "analyze", "expanded.json", "--compact"))
        if expanded["generation"]["added_hypotheses"] != 9:
            raise ValueError("wheel hypothesis-generation smoke failed")
        print(json.dumps({"wheel": wheel.name, "outside_checkout": True,
                          "dependency_isolated": not args.offline,
                          "checks": ["package_data", "full_insert", "paired_end", "html", "replay",
                                     "model_comparison", "hypothesis_generation", "software_self_test"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
