"""Reproducible synthetic sequence-aware analysis benchmark, not biological accuracy."""
from __future__ import annotations

import json
import platform
import statistics
import time
from pathlib import Path

from editwitness import analyze, expand_deletions, load_manifest
from editwitness._version import __version__


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = expand_deletions(load_manifest(root / "examples/sequence-aware.json"))
    analyze(manifest)
    timings = []
    for _ in range(7):
        start = time.perf_counter()
        result = analyze(manifest)
        timings.append(time.perf_counter()-start)
    print(json.dumps({
        "version": __version__, "python": platform.python_version(), "platform": platform.platform(),
        "model": result.model_version, "manifest_sha256": result.manifest_sha256,
        "reference_bases": len(manifest.reference.sequence), "alleles": len(manifest.alleles),
        "hypotheses": len(manifest.hypotheses), "assays": len(manifest.assays+manifest.candidates),
        "exact_products": sum(len(o.products) for o in result.allele_observations),
        "seconds": timings, "median_seconds": statistics.median(timings),
        "scope": "Seven warmed in-process analysis runs; excludes interpreter startup, expansion and file I/O.",
        "caveat": "Synthetic computational workload; not sensitivity, validation, competitor comparison or a throughput guarantee.",
    }, indent=2))


if __name__ == "__main__":
    main()
