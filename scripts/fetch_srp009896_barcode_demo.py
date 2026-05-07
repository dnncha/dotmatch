#!/usr/bin/env python3
"""Fetch/subsample the SRP009896 inline-barcode demultiplexing dataset.

SRP009896 is a maize GBS dataset used in public Cutadapt demultiplexing
examples. The public examples describe 5-prime inline barcodes and 96
demultiplexed outputs for runs such as SRR391079-SRR391082.

The FASTQ files are public through ENA. The barcode/sample sheet is not always
available through ENA metadata, so this script accepts either --barcodes-file or
--barcodes-url and records a clear metadata warning when no barcode file is
provided.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "examples" / "barcode_demux" / "data"
ENA_FIELDS = "run_accession,fastq_ftp,fastq_bytes,fastq_md5,read_count,base_count,sample_alias,experiment_title,study_accession"


def ena_metadata(accession: str) -> dict[str, str]:
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={accession}&result=read_run&fields={ENA_FIELDS}&format=tsv&download=false"
    )
    with urllib.request.urlopen(url, timeout=30) as resp:
        text = resp.read().decode("utf-8")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if not rows:
        raise RuntimeError(f"ENA returned no metadata for {accession}")
    return rows[0]


def https_from_ftp(ftp_path: str) -> str:
    if ftp_path.startswith("ftp://"):
        ftp_path = ftp_path[len("ftp://"):]
    return "https://" + ftp_path


def copy_first_fastq_records(remote_url: str, out: Path, records: int) -> int:
    out.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with urllib.request.urlopen(remote_url, timeout=60) as resp:
        with gzip.GzipFile(fileobj=resp, mode="rb") as gz_in:
            with gzip.open(out, "wt") as gz_out:
                while written < records:
                    header = gz_in.readline()
                    if not header:
                        break
                    seq = gz_in.readline()
                    plus = gz_in.readline()
                    qual = gz_in.readline()
                    if not seq or not plus or not qual:
                        raise RuntimeError("remote FASTQ ended mid-record")
                    gz_out.write(header.decode("utf-8", errors="replace"))
                    gz_out.write(seq.decode("utf-8", errors="replace"))
                    gz_out.write(plus.decode("utf-8", errors="replace"))
                    gz_out.write(qual.decode("utf-8", errors="replace"))
                    written += 1
    return written


def download_full(remote_url: str, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".tmp")
    urllib.request.urlretrieve(remote_url, tmp)
    tmp.replace(out)


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def count_barcodes(path: Path | None) -> tuple[int, list[int]]:
    if path is None or not path.exists():
        return 0, []
    count = 0
    lengths: set[int] = set()
    with path.open() as fh:
        first = fh.readline()
        if not first:
            return 0, []
        delim = "\t" if "\t" in first else ","
        fields = [f.strip().lower() for f in first.rstrip("\n").split(delim)]
        has_header = any(f in {"barcode", "sequence", "index", "sample_id", "sample"} for f in fields)
        rows = fh if has_header else [first, *fh]
        for line in rows:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.rstrip("\n").split(delim)]
            candidates = [p for p in parts if p and all(c.upper() in "ACGTN" for c in p)]
            if not candidates:
                continue
            seq = max(candidates, key=len).upper()
            lengths.add(len(seq))
            count += 1
    return count, sorted(lengths)


def install_barcodes(out_dir: Path, barcodes_file: str | None, barcodes_url: str | None) -> str | None:
    if barcodes_file:
        src = Path(barcodes_file)
        dest = out_dir / "barcodes.tsv"
        shutil.copyfile(src, dest)
        return str(dest)
    if barcodes_url:
        dest = out_dir / "barcodes.tsv"
        with urllib.request.urlopen(barcodes_url, timeout=30) as resp:
            dest.write_bytes(resp.read())
        return str(dest)
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--accession", action="append", default=["SRR391079"])
    parser.add_argument("--subsample", type=int, default=100000, help="records per run; use 0 for full FASTQ download")
    parser.add_argument("--barcodes-file")
    parser.add_argument("--barcodes-url")
    parser.add_argument("--barcode-start", type=int, default=0)
    parser.add_argument("--barcode-length", type=int, default=0, help="expected barcode length; 0 infers from barcode sheet")
    parser.add_argument("--require-barcodes", action="store_true")
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    runs: list[dict[str, object]] = []
    for accession in args.accession:
        meta = ena_metadata(accession)
        fastq_ftp = meta.get("fastq_ftp", "")
        if not fastq_ftp:
            raise RuntimeError(f"ENA metadata for {accession} did not include fastq_ftp")
        remote_url = https_from_ftp(fastq_ftp.split(";")[0])
        out_name = f"{accession}.fastq.gz" if args.subsample == 0 else f"{accession}.subsample{args.subsample}.fastq.gz"
        out_path = args.out / out_name
        written = 0
        if not args.metadata_only:
            if args.subsample == 0:
                download_full(remote_url, out_path)
                written = int(meta.get("read_count") or 0)
            else:
                written = copy_first_fastq_records(remote_url, out_path, args.subsample)
        local_md5 = md5_file(out_path) if out_path.exists() else ""
        local_bytes = out_path.stat().st_size if out_path.exists() else 0
        runs.append({
            "accession": accession,
            "remote_fastq": remote_url,
            "local_fastq": str(out_path),
            "subsample_records": args.subsample,
            "written_records": written,
            "local_md5": local_md5,
            "local_bytes": local_bytes,
            "ena": meta,
        })

    barcode_path = install_barcodes(args.out, args.barcodes_file, args.barcodes_url)
    barcode_count, barcode_lengths = count_barcodes(Path(barcode_path) if barcode_path else None)
    if args.require_barcodes and barcode_path is None:
        raise SystemExit("barcode SOTA fetch requires --barcodes-file or --barcodes-url")
    if args.require_barcodes and barcode_count == 0:
        raise SystemExit("barcode SOTA fetch installed a barcode file but no barcodes could be parsed")
    barcode_length = args.barcode_length or (barcode_lengths[0] if len(barcode_lengths) == 1 else 0)
    metadata = {
        "dataset": "SRP009896 maize GBS inline barcode demultiplexing",
        "barcode_position": "5-prime / read start",
        "barcode_start": args.barcode_start,
        "barcode_length": barcode_length,
        "barcode_count": barcode_count,
        "barcode_lengths": barcode_lengths,
        "sources": [
            "https://biobam.atlassian.net/wiki/spaces/OED0324/pages/3525084904/Reads%2BDemultiplexing%2Bwith%2BCutadapt",
            "https://drive.google.com/file/d/1sxiF4ijqp9jHvFrPa3LWsnxtIHmA0rHJ/view?usp=sharing",
            "https://www.ebi.ac.uk/ena/browser/view/SRP009896",
        ],
        "runs": runs,
        "barcodes": barcode_path,
        "barcodes_required_for_benchmark": barcode_path is None,
        "claim_grade_ready": barcode_path is not None and barcode_count > 0,
        "note": (
            "ENA exposes FASTQ files for the SRP009896 runs. Provide --barcodes-file or --barcodes-url "
            "to install the matching sample/barcode sheet before running state-of-the-art demux benchmarks. "
            "The public Cutadapt example links a Google Drive ExampleDataset.zip that includes the FASTQ files "
            "and barcode file; this script avoids downloading that full archive by default."
        ),
    }
    metadata_path = args.out / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(metadata_path)


if __name__ == "__main__":
    main()
