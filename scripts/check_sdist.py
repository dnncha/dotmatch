"""Verify and test the actual source distribution outside the working checkout."""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    sources = list((root / "dist").glob("*.tar.gz"))
    if len(sources) != 1:
        raise ValueError("build exactly one source distribution in a clean dist/ first")
    with tempfile.TemporaryDirectory(prefix="editwitness-sdist-") as directory:
        work = Path(directory)
        with tarfile.open(sources[0], "r:gz") as archive:
            members = archive.getmembers()
            for member in members:
                target = (work / member.name).resolve()
                if not target.is_relative_to(work.resolve()) or not (member.isfile() or member.isdir()):
                    raise ValueError("unsafe source archive member")
            archive.extractall(work, members=members, filter="data")
        roots = [p for p in work.iterdir() if p.is_dir()]
        if len(roots) != 1:
            raise ValueError("source distribution must contain one root")
        source = roots[0]
        env = dict(os.environ, PYTHONPATH=str(source / "src"), PYTHONNOUSERSITE="1")
        for command in ([sys.executable, "scripts/release_manifest.py", "--check"],
                        [sys.executable, "-m", "pytest", "-q"]):
            subprocess.run(command, cwd=source, env=env, check=True, timeout=300)
        # Prove which package the extracted checkout loads.
        location = subprocess.check_output([sys.executable, "-c", "import editwitness; print(editwitness.__file__)"],
                                           cwd=source, env=env, text=True).strip()
        if not Path(location).resolve().is_relative_to(source.resolve()):
            raise ValueError("source-distribution test used the development package")
        print(json.dumps({"sdist": sources[0].name, "outside_checkout": True,
                          "source_inventory_verified": True, "tests_passed": True}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
