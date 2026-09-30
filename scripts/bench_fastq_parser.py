#!/usr/bin/env python3
"""Synthetic Python parser benchmark; does not measure native assignment."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import platform
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = "python/dotmatch/fastq_io.py"


def worker(source: Path, length: int, records: int) -> dict:
    spec = importlib.util.spec_from_file_location("bench_fastq_source", source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    sequence = ("acgtn" * (length // 5 + 1))[:length]
    quality = ("!I~#5" * (length // 5 + 1))[:length]
    data = (f"@read description\r\n{sequence}\r\n+read\r\n{quality}\r\n" * records)
    expected_sequence = sequence.upper()
    handle = io.StringIO(data, newline="")
    started = time.perf_counter()
    count = 0
    for record in module.iter_fastq_records(handle, "synthetic.fastq"):
        assert (record.read_id, record.seq, record.qual) == ("read", expected_sequence, quality)
        count += 1
    elapsed = time.perf_counter() - started
    assert count == records
    return {
        "length": length, "records": records, "wall_seconds": elapsed,
        "process_peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024),
        "input_sha256": hashlib.sha256(data.encode()).hexdigest(),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "python": platform.python_version(), "platform": platform.platform(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-ref", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--records", type=int, default=20000)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--length", type=int, default=150, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.records <= 0 or args.repetitions <= 0:
        parser.error("records and repetitions must be positive")
    if args.worker:
        import json
        print(json.dumps(worker(args.worker, args.length, args.records)))
        return
    import json
    rows = []
    with tempfile.TemporaryDirectory() as directory:
        baseline = Path(directory) / "baseline.py"
        baseline.write_bytes(subprocess.check_output(["git", "show", f"{args.baseline_ref}:{SOURCE}"], cwd=ROOT))
        for length in (50, 150, 300):
            for repetition in range(args.repetitions):
                variants = [("baseline", baseline), ("candidate", ROOT / SOURCE)]
                for variant, source in variants[::1 if repetition % 2 == 0 else -1]:
                    result = subprocess.check_output([
                        sys.executable, str(Path(__file__).resolve()), "--baseline-ref", args.baseline_ref,
                        "--output", str(args.output), "--worker", str(source), "--length", str(length),
                        "--records", str(args.records),
                    ], text=True)
                    rows.append({"variant": variant, "repetition": repetition + 1, **json.loads(result)})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
