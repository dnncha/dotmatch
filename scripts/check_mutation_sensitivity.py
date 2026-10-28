#!/usr/bin/env python3
"""Targeted mutation gate for the scientific invariants.

Each mutant is applied to a throwaway copy of the package, then the scientific
tests must fail. A surviving mutant means the invariant tests do not distinguish
a known-wrong semantic change from correct behaviour.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEST_FILES = ("test_scientific_invariants.py", "test_scientific_fuzz.py")


@dataclass(frozen=True)
class Mutant:
    name: str
    relative_path: str
    before: str
    after: str


MUTANTS = (
    Mutant(
        "radius-becomes-unique",
        "core.py",
        "                    MATCH_AMBIGUOUS,",
        "                    MATCH_UNIQUE,",
    ),
    Mutant(
        "fastq-validation-disabled",
        "fastq_io.py",
        "        if not _valid_printable_ascii(seq):",
        "        if False:",
    ),
    Mutant(
        "counts-parse-through-float",
        "count_io.py",
        "        return int(text)",
        "        return int(float(text))",
    ),
    Mutant(
        "ambiguous-features-counted",
        "feature_matrix.py",
        "                        if status == MATCH_UNIQUE and target_id:\n                            counts[(observation.cell_barcode, target_id)] += 1",
        "                        if status in {MATCH_UNIQUE, MATCH_AMBIGUOUS}:\n                            counts[(observation.cell_barcode, target_id or f\"ambiguous_{result.target_index}\")] += 1",
    ),
)


def _copy_source(destination: Path) -> None:
    shutil.copytree(ROOT / "python" / "dotmatch", destination / "python" / "dotmatch")
    shutil.copytree(ROOT / "python" / "quickdna", destination / "python" / "quickdna")
    (destination / "tests").mkdir()
    for name in TEST_FILES:
        shutil.copy2(ROOT / "python" / "tests" / name, destination / "tests" / name)


def _library_path() -> Path:
    for name in ("libdotmatch.dylib", "libdotmatch.so"):
        candidate = ROOT / name
        if candidate.is_file():
            return candidate
    raise SystemExit("build the shared library first: make shared")


def _run_mutant(mutant: Mutant, library: Path) -> bool:
    with tempfile.TemporaryDirectory(prefix="dotmatch-mutant-") as temporary:
        root = Path(temporary)
        _copy_source(root)
        target = root / "python" / "dotmatch" / mutant.relative_path
        source = target.read_text(encoding="utf-8")
        if mutant.before not in source:
            raise SystemExit(f"mutant {mutant.name}: source pattern not found in {mutant.relative_path}")
        target.write_text(source.replace(mutant.before, mutant.after, 1), encoding="utf-8")
        env = {
            **os.environ,
            "PYTHONPATH": str(root / "python"),
            "DOTMATCH_LIB": str(library),
        }
        completed = subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                str(root / "tests"),
                "-q",
                "-c",
                os.devnull,
                "-p",
                "no:cacheprovider",
            ],
            cwd=root,
            env=env,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if completed.returncode == 0:
            print(f"survived: {mutant.name}", file=sys.stderr)
            print(completed.stdout[-4000:], file=sys.stderr)
            return False
        print(f"killed: {mutant.name}")
        return True


def main() -> int:
    library = _library_path()
    killed = [_run_mutant(mutant, library) for mutant in MUTANTS]
    if not all(killed):
        return 1
    print(f"all {len(MUTANTS)} targeted scientific mutants were killed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
