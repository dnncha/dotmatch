#!/usr/bin/env python3
"""Fetch or create the MAGeCK/Yusa CRISPR demo inputs for DotMatch examples."""

from __future__ import annotations

import argparse
import gzip
import shutil
import urllib.request
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "examples" / "crispr_guides" / "data"
YUSA_URL = "https://sourceforge.net/projects/mageck/files/libraries/yusa_library.csv.zip/download"
FASTQ_URLS = {
    "ERR376998": "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR376/ERR376998/ERR376998.fastq.gz",
    "ERR376999": "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR376/ERR376999/ERR376999.fastq.gz",
}


def write_small(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "yusa_library.csv").write_text(
        "id,gRNA.sequence,Gene\n"
        "guide_exact,ACGTACGTACGTACGTACG,GENE1\n"
        "guide_sub,TTTTCCCCAAAAGGGGTTT,GENE2\n"
        "guide_other,GGGGAAAACCCCTTTTGGG,GENE3\n",
        encoding="utf-8",
    )
    reads = (
        "@exact\n"
        "NNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTACGTACG\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        "@one_sub\n"
        "NNNNNNNNNNNNNNNNNNNNNNNTTTTCCCCAAAAGGGGTTA\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        "@no_match\n"
        "NNNNNNNNNNNNNNNNNNNNNNNCCCCCCCCCCCCCCCCCCC\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
    )
    with gzip.open(out_dir / "ERR376998.fastq.gz", "wt", encoding="utf-8") as fh:
        fh.write(reads)
    with gzip.open(out_dir / "ERR376999.fastq.gz", "wt", encoding="utf-8") as fh:
        fh.write(reads)


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as response:
        expected_size = int(response.headers.get("Content-Length") or 0)
        if dest.exists() and expected_size > 0 and dest.stat().st_size == expected_size:
            return
        if dest.exists() and expected_size == 0:
            return
        tmp = dest.with_suffix(dest.suffix + ".tmp")
        with tmp.open("wb") as fh:
            shutil.copyfileobj(response, fh)
        tmp.replace(dest)


def fetch_library(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    library = out_dir / "yusa_library.csv"
    if library.exists() and sum(1 for _ in library.open(encoding="utf-8")) > 100:
        return
    zip_path = out_dir / "yusa_library.csv.zip"
    if zip_path.exists():
        zip_path.unlink()
    download(YUSA_URL, zip_path)
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(out_dir)


def download_fastq_subsample(url: str, dest: Path, records: int) -> None:
    if dest.exists() and count_fastq_records(dest) == records:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as response:
        with gzip.GzipFile(fileobj=response, mode="rb") as inp:
            with gzip.open(dest, "wb") as out:
                for _ in range(records * 4):
                    line = inp.readline()
                    if not line:
                        break
                    out.write(line)


def count_fastq_records(path: Path) -> int:
    if not path.exists():
        return 0
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return sum(1 for _ in fh) // 4


def fetch_subsample(out_dir: Path, records: int) -> None:
    fetch_library(out_dir)
    for sample, url in FASTQ_URLS.items():
        download_fastq_subsample(url, out_dir / f"{sample}.fastq.gz", records)


def fetch_full(out_dir: Path) -> None:
    fetch_library(out_dir)
    for sample, url in FASTQ_URLS.items():
        download(url, out_dir / f"{sample}.fastq.gz")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--small", action="store_true", help="write a tiny local fixture instead of downloading full FASTQ files")
    parser.add_argument("--subsample", type=int, default=0, help="download this many real FASTQ records per sample")
    args = parser.parse_args()
    out_dir = Path(args.out)
    if args.small:
        write_small(out_dir)
    elif args.subsample > 0:
        fetch_subsample(out_dir, args.subsample)
    else:
        fetch_full(out_dir)
    print(out_dir)


if __name__ == "__main__":
    main()
