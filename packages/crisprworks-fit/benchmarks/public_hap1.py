"""Fetch a pinned public full-library screen and run the paired benchmark."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from urllib.request import urlopen

COMMIT = "53388adbb4fb0931e5c9dda135502be19e4555f0"
URL = f"https://raw.githubusercontent.com/hart-lab/bagel/{COMMIT}/reads_hap1.txt"
SHA256 = "7638bc6237cf8a3e2302fcd969d014b041819669b37f2f87bc1f7cfbc45ca8a7"
DESIGN = ("Samples\tbaseline\ttreatment\nHAP1_T0\t1\t0\n"
          "HAP1_T18A\t1\t1\nHAP1_T18B\t1\t1\nHAP1_T18C\t1\t1\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--cohort", choices=("full", "four-guide"), default="full",
                        help="Full table or complete four-guide gene labels; selection is recorded")
    parser.add_argument("--local-counts", type=Path,
                        help="Use a previously downloaded file; the same hash is required")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.local_counts:
        data = args.local_counts.read_bytes()
    else:
        with urlopen(URL, timeout=120) as response:
            data = response.read()
    if hashlib.sha256(data).hexdigest() != SHA256:
        raise ValueError("Public count-table hash does not match the pinned dataset")
    rows = data.decode().splitlines()
    guide_counts = Counter(row.split("\t")[1] for row in rows[1:])
    if args.cohort == "four-guide":
        rows = [rows[0]] + [row for row in rows[1:] if guide_counts[row.split("\t")[1]] == 4]
        data = ("\n".join(rows) + "\n").encode()
        guide_counts = {gene: count for gene, count in guide_counts.items() if count == 4}
    counts = args.out_dir / "counts.tsv"
    design = args.out_dir / "design.tsv"
    counts.write_bytes(data)
    design.write_text(DESIGN)
    provenance = {
        "source_url": URL, "source_commit": COMMIT, "source_sha256": SHA256,
        "counts_sha256": hashlib.sha256(data).hexdigest(), "cohort": args.cohort,
        "guides": len(rows) - 1, "gene_labels_including_controls": len(guide_counts),
        "selection": "all labels" if args.cohort == "full" else "labels having exactly four guides; incomplete labels and larger control bins excluded before normalization",
        "design": "one T0 baseline; three T18 replicates; single treatment effect",
        "scope": "numerical parity against MAGeCK2, not validation of biological hits",
        "control_caveat": "Upstream default skips fitting labels with >=40 guides; their permutation tails use the preceding guide-count group",
        "source_license": "MIT; dataset retrieved separately, not redistributed in this package",
    }
    (args.out_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    subprocess.run([
        sys.executable, str(Path(__file__).with_name("run.py")),
        "--out-dir", str(args.out_dir / "paired"), "--count-table", str(counts),
        "--design-matrix", str(design), "--update-efficiency",
        "--genes-varmodeling", "1000", "--repeats", str(args.repeats),
    ], check=True)
    print((args.out_dir / "paired" / "benchmark.json").read_text(), flush=True)


if __name__ == "__main__":
    main()
