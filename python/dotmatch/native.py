from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Sequence


def native_cli_candidates() -> list[Path]:
    env = os.environ.get("DOTMATCH_NATIVE_CLI")
    candidates = [Path(env)] if env else []
    here = Path(__file__).resolve()
    candidates.extend(
        [
            here.parent / "dotmatch-native",
            here.parents[2] / "dotmatch",
            Path.cwd() / "dotmatch",
        ]
    )
    return candidates


def find_native_cli() -> Path:
    for path in native_cli_candidates():
        if path.exists() and os.access(path, os.X_OK):
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
