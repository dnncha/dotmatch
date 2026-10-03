"""Run a validated JSON screen list through the installed Fit command."""

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys


def commands(screens_path, output_dir, workers, seed, rounds, overwrite=False, permutation_pvalues="legacy",
             write_fit_details=False, update_efficiency=False):
    screens = json.loads(screens_path.read_text())
    if not isinstance(screens, list) or not screens:
        raise ValueError("Screen manifest must be a nonempty JSON list")
    seen = set()
    tasks = []
    if permutation_pvalues not in ("legacy", "finite"):
        raise ValueError("Permutation p-values must be legacy or finite")
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
                   "--seed", str(seed), "--permutation-round", str(rounds),
                   "--permutation-pvalues", permutation_pvalues]
        if write_fit_details:
            command.append("--write-fit-details")
        if update_efficiency:
            command.append("--update-efficiency")
        if "control_gene" in screen:
            controls = (screens_path.parent / screen["control_gene"]).resolve()
            if not controls.is_file():
                raise ValueError(f"Missing control-gene file for {name}")
            command.extend(["--control-gene", str(controls), "--norm-method", "control"])
        tasks.append((name, command))
    return tasks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screens", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--permutation-round", type=int, default=10)
    parser.add_argument("--permutation-pvalues", choices=("legacy", "finite"), default="legacy")
    parser.add_argument("--write-fit-details", action="store_true", help="Retain named full-precision effects for calibration")
    parser.add_argument("--update-efficiency", action="store_true", help="Estimate guide efficiency during MLE fitting")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.threads < 1 or args.permutation_round < 1 or not 0 <= args.seed < 2**32:
        parser.error("threads/rounds must be positive and seed must fit uint32")
    try:
        tasks = commands(args.screens.resolve(), args.out_dir.resolve(), args.threads,
                         args.seed, args.permutation_round, args.overwrite, args.permutation_pvalues,
                         args.write_fit_details, args.update_efficiency)
    except (ValueError, KeyError, TypeError, OSError) as error:
        parser.error(str(error))
    for name, command in tasks:
        print(json.dumps({"screen": name, "command": command}), flush=True)
        if not args.dry_run:
            subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
