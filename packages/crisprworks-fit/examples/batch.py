"""Run a validated JSON screen list through the installed Fit command."""

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys


def commands(screens_path, output_dir, workers, seed, rounds, overwrite=False):
    screens = json.loads(screens_path.read_text())
    if not isinstance(screens, list) or not screens:
        raise ValueError("Screen manifest must be a nonempty JSON list")
    seen = set()
    tasks = []
    for screen in screens:
        name = screen["name"]
        if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name):
            raise ValueError("Screen names must be simple file-name components starting with a letter or digit")
        if name in seen:
            raise ValueError(f"Duplicate screen name: {name}")
        seen.add(name)
        counts = (screens_path.parent / screen["count_table"]).resolve()
        design = (screens_path.parent / screen["design_matrix"]).resolve()
        if not counts.is_file() or not design.is_file():
            raise ValueError(f"Missing count table or design for {name}")
        prefix = output_dir / name / "screen"
        if prefix.parent.exists() and not overwrite:
            raise ValueError(f"Output directory exists for {name}; use a new root or --overwrite")
        command = [sys.executable, "-m", "crisprworks_fit", "mle",
                   "-k", str(counts), "-d", str(design), "-n", str(prefix),
                   "--threads", str(workers), "--blas-threads", "1",
                   "--seed", str(seed), "--permutation-round", str(rounds)]
        tasks.append((name, command))
    return tasks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screens", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--permutation-round", type=int, default=10)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.threads < 1 or args.permutation_round < 1 or not 0 <= args.seed < 2**32:
        parser.error("threads/rounds must be positive and seed must fit uint32")
    try:
        tasks = commands(args.screens.resolve(), args.out_dir.resolve(), args.threads,
                         args.seed, args.permutation_round, args.overwrite)
    except (ValueError, KeyError, TypeError, OSError) as error:
        parser.error(str(error))
    for name, command in tasks:
        print(json.dumps({"screen": name, "command": command}), flush=True)
        if not args.dry_run:
            subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
