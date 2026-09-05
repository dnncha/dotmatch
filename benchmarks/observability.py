"""Reproducible synthetic workloads; timings are not biological validation."""
from __future__ import annotations

import json
import platform
import statistics
import time
from importlib.resources import files

from editwitness import Manifest, __version__, analyze
from editwitness.scan import scan_deletions


def main() -> None:
    manifest = Manifest.model_validate_json(files("editwitness").joinpath("data/demo.json").read_text())
    data = manifest.model_dump(mode="json")
    data["deletion_scan"] = {"start_min": 0, "start_max": 499, "end_min": 1, "end_max": 900}
    scan_manifest = Manifest.model_validate(data)
    timings: dict[str, list[float]] = {}
    counts = {}
    for name, operation in (("demo_analysis", lambda: analyze(manifest)),
                            ("deletion_grid_450000_pairs", lambda: scan_deletions(scan_manifest))):
        operation()  # Explicit warm-up, not included in timing.
        times = []
        for _ in range(5):
            start = time.perf_counter()
            result = operation()
            times.append(time.perf_counter() - start)
        timings[name] = times
        if name.startswith("deletion"):
            counts["enumerated_valid_deletions"] = result.enumerated_deletions
    print(json.dumps({"version": __version__, "python": platform.python_version(),
                      "platform": platform.platform(), "processor": platform.processor(),
                      "warmups_per_workload": 1, "repetitions": 5, "seconds": timings,
                      "median_seconds": {k: statistics.median(v) for k, v in timings.items()},
                      **counts, "caveat": "Synthetic in-process warm timings on this machine only. No competitor comparison or empirical assay benchmark."}, indent=2))


if __name__ == "__main__":
    main()
