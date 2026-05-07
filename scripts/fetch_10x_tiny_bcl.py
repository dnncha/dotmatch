#!/usr/bin/env python3
"""Fetch the public 10x Genomics tiny-BCL mkfastq demo dataset.

The run folder is a real public BCL example used by 10x Genomics to document
Cell Ranger mkfastq. It is suitable for parser/workflow validation, but it is
not by itself a state-of-the-art benchmark because it is small and competitor
rows still need to be generated on the same machine.
"""

from __future__ import annotations

import argparse
import csv
import tarfile
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "examples" / "bcl_demux" / "data"

URLS = {
    "archive": "https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz",
    "samplesheet": "https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-samplesheet-1.2.0.csv",
    "simple_csv": "https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv",
    "chromium_i7": "https://cdn.10xgenomics.com/raw/upload/v1655155613/support/in-line%20documents/chromium-shared-sample-indexes-plate.csv",
}

LEGACY_10X_INDEX_SETS = {
    # The public tiny-BCL sample sheet uses the old Cell Ranger SI-P03-C9
    # alias. Cell Ranger expands this kind of 10x index-set name before
    # invoking bcl2fastq. The component oligos below are from the official
    # Chromium i7 multiplex kit CSV downloaded by this script.
    "SI-P03-C9": ["SI-GA-E3", "SI-GA-F3", "SI-GA-G3", "SI-GA-H3"],
}


def download(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        print(dest)
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "DotMatch benchmark fetcher"})
    with urllib.request.urlopen(req) as response, dest.open("wb") as fh:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            fh.write(chunk)
    print(dest)


def load_10x_i7(path: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with path.open(newline="") as fh:
        for row in csv.reader(fh):
            if len(row) >= 5:
                out[row[0]] = row[1:5]
    return out


def write_normalized_tiny_samplesheet(src: Path, index_csv: Path, dest: Path) -> None:
    index_sets = load_10x_i7(index_csv)
    rows: list[dict[str, str]] = []
    in_data = False
    header: list[str] | None = None
    with src.open(newline="") as fh:
        for raw in csv.reader(fh):
            if not raw:
                continue
            first = raw[0].strip()
            if first in ("[Data]", "[BCLConvert_Data]"):
                in_data = True
                header = None
                continue
            if not in_data:
                continue
            if header is None:
                header = raw
                continue
            record = {header[i]: raw[i] if i < len(raw) else "" for i in range(len(header))}
            index_name = record.get("index", "")
            component_sets = LEGACY_10X_INDEX_SETS.get(index_name)
            if component_sets is None:
                rows.append(record)
                continue
            for component in component_sets:
                for seq in index_sets[component]:
                    expanded = dict(record)
                    expanded["index"] = seq
                    expanded["I7_Index_ID"] = component
                    rows.append(expanded)

    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["[Header]"])
        writer.writerow(["IEMFileVersion", "4"])
        writer.writerow(["[Data]"])
        writer.writerow(["Lane", "Sample_ID", "Sample_Name", "index", "Sample_Project"])
        for row in rows:
            writer.writerow([
                row.get("Lane", "1") or "1",
                row.get("Sample_ID", row.get("Sample", "sample")) or "sample",
                row.get("Sample_Name", row.get("Sample", "sample")) or "sample",
                row.get("index", ""),
                row.get("Sample_Project", ""),
            ])
    print(dest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-only", action="store_true", help="download sample sheets only")
    parser.add_argument("--extract", action="store_true", help="extract the run folder after downloading")
    args = parser.parse_args()

    DATA.mkdir(parents=True, exist_ok=True)
    download(URLS["samplesheet"], DATA / "cellranger-tiny-bcl-samplesheet-1.2.0.csv")
    download(URLS["simple_csv"], DATA / "cellranger-tiny-bcl-simple-1.2.0.csv")
    download(URLS["chromium_i7"], DATA / "chromium-shared-sample-indexes-plate.csv")
    write_normalized_tiny_samplesheet(
        DATA / "cellranger-tiny-bcl-samplesheet-1.2.0.csv",
        DATA / "chromium-shared-sample-indexes-plate.csv",
        DATA / "cellranger-tiny-bcl-samplesheet.normalized.csv",
    )
    if args.metadata_only:
        return

    archive = DATA / "cellranger-tiny-bcl-1.2.0.tar.gz"
    download(URLS["archive"], archive)
    if args.extract:
        with tarfile.open(archive, "r:gz") as tf:
            tf.extractall(DATA)
        print(DATA / "cellranger-tiny-bcl-1.2.0")


if __name__ == "__main__":
    main()
