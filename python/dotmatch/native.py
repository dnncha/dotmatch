from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Sequence


def native_cli_candidates() -> list[Path]:
    env = os.environ.get("DOTMATCH_NATIVE_CLI")
    if env:
        return [Path(env).expanduser().resolve()]
    here = Path(__file__).resolve()
    candidates = [here.parent / "dotmatch-native"]
    if here.parent.name == "dotmatch" and here.parent.parent.name == "python":
        candidates.append(here.parents[2] / "dotmatch")
    return candidates


def find_native_cli() -> Path:
    for path in native_cli_candidates():
        if path.is_file() and os.access(path, os.X_OK):
            return path
    searched = ", ".join(str(path) for path in native_cli_candidates())
    raise FileNotFoundError(
        "could not find the DotMatch native CLI; searched: "
        f"{searched}. Build it with `make dotmatch`, install a wheel with the "
        "bundled native executable, or set DOTMATCH_NATIVE_CLI=/path/to/dotmatch."
    )


def run_native_cli(argv: Sequence[str]) -> int:
    native = find_native_cli()
    completed = subprocess.run([str(native), *argv], check=False)
    return int(completed.returncode)
