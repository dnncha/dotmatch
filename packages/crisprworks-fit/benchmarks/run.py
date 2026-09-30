"""Reproducible paired CLI benchmark; reports measurements, not guarantees."""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--genes", type=int, default=1000)
    parser.add_argument("--guides", type=int, default=6)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--count-table", type=Path)
    parser.add_argument("--design-matrix", type=Path)
    parser.add_argument("--update-efficiency", action="store_true")
    parser.add_argument("--genes-varmodeling", type=int, default=0)
    args = parser.parse_args()
    if min(args.genes, args.guides, args.repeats) < 1:
        parser.error("genes, guides and repeats must be positive")
    from crisprworks_fit.cli import synthetic_screen
    from crisprworks_fit.backend import REFERENCE_COMMIT
    if (args.count_table is None) != (args.design_matrix is None):
        parser.error("count-table and design-matrix must be supplied together")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.count_table is None:
        counts, design = synthetic_screen(args.out_dir / "inputs", genes=args.genes, guides=args.guides)
        dataset = "synthetic NB counts; not biological validation"
        genes, guides = args.genes, args.guides
        samples = 6
    else:
        counts, design = args.count_table, args.design_matrix
        lines = [line.split() for line in counts.read_text().splitlines()[1:] if line.strip()]
        genes = len({line[1] for line in lines})
        guides = None
        samples = len(lines[0]) - 2
        dataset = "supplied count table; provenance and validation must be recorded separately"
    report = {
        "dataset": dataset, "genes": genes, "guides_per_gene": guides,
        "samples": samples, "permutation_rounds": 2, "seed": 42,
        "worker_processes": 1, "blas_threads_per_process": 1,
        "reference_commit": REFERENCE_COMMIT,
        "python": platform.python_version(), "platform": platform.platform(),
        "counts_sha256": hashlib.sha256(counts.read_bytes()).hexdigest(),
        "design_sha256": hashlib.sha256(design.read_bytes()).hexdigest(), "runs": [],
    }
    summaries = []
    precisions = []
    for repeat in range(args.repeats):
        # Alternate order to reduce a systematic warm-cache advantage.
        order = ("reference", "accelerated") if repeat % 2 == 0 else ("accelerated", "reference")
        for backend in order:
            prefix = args.out_dir / f"{backend}-{repeat}"
            command = [
                sys.executable, "-m", "crisprworks_fit", "mle", "-k", str(counts),
                "-d", str(design), "-n", str(prefix), "--backend", backend,
                "--threads", "1", "--blas-threads", "1", "--seed", "42", "--permutation-round", "2",
                "--write-fit-details", "--genes-varmodeling", str(args.genes_varmodeling),
            ]
            if args.update_efficiency:
                command.append("--update-efficiency")
            started = time.perf_counter()
            completed = subprocess.run(command, capture_output=True, text=True, timeout=600)
            elapsed = time.perf_counter() - started
            Path(str(prefix) + ".stderr.txt").write_text(completed.stderr)
            if completed.returncode:
                raise RuntimeError(completed.stderr)
            manifest = json.loads(Path(str(prefix) + ".crisprworks.json").read_text())
            summaries.append(Path(str(prefix) + ".gene_summary.txt").read_bytes())
            precisions.append(json.loads(Path(str(prefix) + ".fit-details.json").read_text())["genes"])
            report["runs"].append({
                "backend": backend, "repeat": repeat, "wall_seconds_including_startup": elapsed,
                "workflow_seconds": manifest["elapsed_seconds"], "command": command,
                "manifest": manifest,
            })
            print(f"{backend} repeat {repeat}: {elapsed:.3f}s", flush=True)
    report["identical_gene_summary_bytes"] = all(item == summaries[0] for item in summaries)
    maximum_errors = {}
    for candidate in precisions[1:]:
        if candidate.keys() != precisions[0].keys():
            raise RuntimeError("Gene identifiers changed")
        for gene_id, expected in precisions[0].items():
            if candidate[gene_id]["guides"] != expected["guides"]:
                raise RuntimeError("Guide count changed")
            for field, values in expected.items():
                if field == "guides":
                    continue
                actual = np.asarray(candidate[gene_id][field])
                reference = np.asarray(values)
                np.testing.assert_allclose(actual, reference, rtol=1e-7, atol=1e-8,
                                           err_msg=f"{gene_id}: {field}")
                error = float(np.max(np.abs(actual - reference))) if actual.size else 0
                maximum_errors[field] = max(maximum_errors.get(field, 0), error)
    report["full_precision_max_absolute_errors"] = maximum_errors
    report["full_precision_tolerances"] = {"rtol": 1e-7, "atol": 1e-8}
    medians = {
        backend: statistics.median(run["wall_seconds_including_startup"] for run in report["runs"]
                                   if run["backend"] == backend)
        for backend in ("reference", "accelerated")
    }
    report["median_wall_seconds"] = medians
    report["median_speedup"] = medians["reference"] / medians["accelerated"]
    output = args.out_dir / "benchmark.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in (
        "identical_gene_summary_bytes", "median_wall_seconds", "median_speedup",
    )}, indent=2))
    if not report["identical_gene_summary_bytes"]:
        raise RuntimeError("Gene summary bytes differed; investigate before publishing a speed claim")


if __name__ == "__main__":
    main()
