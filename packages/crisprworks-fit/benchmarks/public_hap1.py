"""Fetch a pinned public full-library screen and run the paired benchmark."""

import argparse
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
    counts = args.out_dir / "counts.tsv"
    design = args.out_dir / "design.tsv"
    counts.write_bytes(data)
    design.write_text(DESIGN)
    provenance = {
        "source_url": URL, "source_commit": COMMIT, "counts_sha256": SHA256,
        "guides": 71090, "gene_labels_including_controls": 18056,
        "design": "one T0 baseline; three T18 replicates; single treatment effect",
        "scope": "numerical parity against MAGeCK2, not validation of biological hits",
        "control_caveat": "Upstream default skips labels with >=10 guides during permutations",
        "source_license": "MIT; dataset retrieved separately, not redistributed in this package",
    }
    (args.out_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    subprocess.run([
        sys.executable, str(Path(__file__).with_name("run.py")),
        "--out-dir", str(args.out_dir / "paired"), "--count-table", str(counts),
        "--design-matrix", str(design), "--update-efficiency",
        "--genes-varmodeling", "1000", "--repeats", str(args.repeats),
    ], check=True)


if __name__ == "__main__":
    main()
