#!/usr/bin/env python3
"""Fetch/subsample the Sanson/Brunello public CRISPR benchmark dataset.

This is the second real CRISPR dataset for DotMatch's SOTA evidence track. It
uses the Sanson et al. Brunello/mod-tracr reads from PRJNA508200, matching the
sample structure used in guide-counter's public README benchmark:

  plasmid: SRR8297997
  RepA:    SRR8297837 + SRR8297836
  RepB:    SRR8297839 + SRR8297838
  RepC:    SRR8297841 + SRR8297840

For subsampled runs the script writes the first N records per biological
sample, spanning source runs in order. For full runs it concatenates the gzip
members into one FASTQ.gz per sample, which is valid gzip and keeps the
benchmark interface one FASTQ per sample.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import urllib.request
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "examples" / "crispr_sanson_brunello" / "data"
ENA_FIELDS = (
    "run_accession,sample_alias,experiment_alias,experiment_title,fastq_ftp,"
    "read_count,base_count,fastq_bytes,fastq_md5"
)
LIB_URL = (
    "https://sourceforge.net/projects/mageck/files/libraries/"
    "broadgpp-brunello-library-corrected.txt.zip/download"
)

SAMPLES: dict[str, list[str]] = {
    "plasmid": ["SRR8297997"],
    "RepA": ["SRR8297837", "SRR8297836"],
    "RepB": ["SRR8297839", "SRR8297838"],
    "RepC": ["SRR8297841", "SRR8297840"],
}


def selected_samples(sample_text: str) -> dict[str, list[str]]:
    names = [part.strip() for part in sample_text.split(",") if part.strip()]
    if not names:
        raise ValueError("at least one sample is required")
    selected: dict[str, list[str]] = {}
    for name in names:
        if name not in SAMPLES:
            raise ValueError(f"unknown sample: {name}")
        selected[name] = SAMPLES[name]
    return selected


def ena_metadata(accession: str) -> dict[str, str]:
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={accession}&result=read_run&fields={ENA_FIELDS}&format=tsv&download=false"
    )
    with urllib.request.urlopen(url, timeout=45) as resp:
        text = resp.read().decode("utf-8")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if not rows:
        raise RuntimeError(f"ENA returned no metadata for {accession}")
    return rows[0]


def https_from_ftp(ftp_path: str) -> str:
    if ftp_path.startswith("ftp://"):
        ftp_path = ftp_path[len("ftp://"):]
    return "https://" + ftp_path


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def count_fastq_gz(path: Path) -> int:
    with gzip.open(path, "rt") as fh:
        return sum(1 for _ in fh) // 4


def cached_fastq_records(path: Path, requested_records: int) -> int:
    if requested_records <= 0 or not path.exists() or path.stat().st_size == 0:
        return 0
    try:
        records = count_fastq_gz(path)
    except (OSError, EOFError, gzip.BadGzipFile, UnicodeDecodeError):
        return 0
    return records if records >= requested_records else 0


def fetch_library(out_dir: Path) -> Path:
    out = out_dir / "broadgpp-brunello-library-corrected.txt"
    if out.exists() and out.stat().st_size > 0:
        return out
    with urllib.request.urlopen(LIB_URL, timeout=90) as resp:
        data = resp.read()
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        wanted = [n for n in zf.namelist() if n.endswith("broadgpp-brunello-library-corrected.txt")]
        if not wanted:
            raise RuntimeError("Brunello library zip did not contain the expected library file")
        out.write_bytes(zf.read(wanted[0]))
    return out


def copy_first_records_from_urls(urls: list[str], out: Path, records: int) -> int:
    written = 0
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt") as gz_out:
        for url in urls:
            if written >= records:
                break
            with urllib.request.urlopen(url, timeout=120) as resp:
                with gzip.GzipFile(fileobj=resp, mode="rb") as gz_in:
                    while written < records:
                        header = gz_in.readline()
                        if not header:
                            break
                        seq = gz_in.readline()
                        plus = gz_in.readline()
                        qual = gz_in.readline()
                        if not seq or not plus or not qual:
                            raise RuntimeError(f"remote FASTQ ended mid-record while reading {url}")
                        gz_out.write(header.decode("utf-8", errors="replace"))
                        gz_out.write(seq.decode("utf-8", errors="replace"))
                        gz_out.write(plus.decode("utf-8", errors="replace"))
                        gz_out.write(qual.decode("utf-8", errors="replace"))
                        written += 1
    return written


def copy_checked_url(url: str, dest, expected_md5: str = "") -> int:
    h = hashlib.md5()
    copied = 0
    with urllib.request.urlopen(url, timeout=120) as src:
        while True:
            chunk = src.read(1024 * 1024)
            if not chunk:
                break
            dest.write(chunk)
            h.update(chunk)
            copied += len(chunk)
    if expected_md5 and h.hexdigest().lower() != expected_md5.lower():
        raise RuntimeError(f"MD5 mismatch while downloading {url}")
    return copied


def concatenate_full(urls: list[str], out: Path, expected_md5s: list[str] | None = None,
                     attempts: int = 3) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".tmp")
    last_error: Exception | None = None
    for _ in range(max(1, attempts)):
        try:
            with tmp.open("wb") as dest:
                for i, url in enumerate(urls):
                    expected_md5 = expected_md5s[i] if expected_md5s is not None and i < len(expected_md5s) else ""
                    copy_checked_url(url, dest, expected_md5)
            tmp.replace(out)
            return
        except Exception as exc:
            last_error = exc
            tmp.unlink(missing_ok=True)
    raise RuntimeError(f"failed to download complete FASTQ for {out}: {last_error}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--subsample", type=int, default=100000,
            help="records per biological sample; use 0 for full FASTQ download/concatenation")
    parser.add_argument("--samples", default=",".join(SAMPLES),
            help="comma-separated biological samples to fetch: plasmid,RepA,RepB,RepC")
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args()
    try:
        samples = selected_samples(args.samples)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    args.out.mkdir(parents=True, exist_ok=True)
    library = fetch_library(args.out)
    sample_rows: list[dict[str, object]] = []
    run_rows: list[dict[str, object]] = []

    for sample, accessions in samples.items():
        metas = [ena_metadata(acc) for acc in accessions]
        urls = [https_from_ftp((m.get("fastq_ftp") or "").split(";")[0]) for m in metas]
        if any(not u for u in urls):
            raise RuntimeError(f"missing FASTQ URL for {sample}")
        out_name = f"{sample}.fastq.gz" if args.subsample == 0 else f"{sample}.subsample{args.subsample}.fastq.gz"
        out_path = args.out / out_name
        written = 0
        local_md5 = ""
        local_size = 0
        if not args.metadata_only:
            if args.subsample == 0:
                expected_records = sum(int(m.get("read_count") or 0) for m in metas)
                written = cached_fastq_records(out_path, expected_records)
                if written == 0:
                    out_path.unlink(missing_ok=True)
                    expected_md5s = [(m.get("fastq_md5") or "").split(";")[0] for m in metas]
                    concatenate_full(urls, out_path, expected_md5s)
                    written = count_fastq_gz(out_path)
                    if expected_records and written != expected_records:
                        out_path.unlink(missing_ok=True)
                        raise RuntimeError(
                            f"downloaded FASTQ for {sample} has {written} records; expected {expected_records}"
                        )
            else:
                written = cached_fastq_records(out_path, args.subsample)
                if written == 0:
                    written = copy_first_records_from_urls(urls, out_path, args.subsample)
            local_md5 = md5_file(out_path)
            local_size = out_path.stat().st_size
        for meta, url in zip(metas, urls):
            run_rows.append({
                "sample": sample,
                "accession": meta.get("run_accession", ""),
                "remote_fastq": url,
                "read_count": meta.get("read_count", ""),
                "fastq_bytes": meta.get("fastq_bytes", ""),
                "fastq_md5": meta.get("fastq_md5", ""),
                "sample_alias": meta.get("sample_alias", ""),
                "experiment_title": meta.get("experiment_title", ""),
            })
        sample_rows.append({
            "sample_id": sample,
            "fastq": str(out_path),
            "source_accessions": accessions,
            "remote_fastqs": urls,
            "subsample_records": args.subsample,
            "written_records": written,
            "local_md5": local_md5,
            "local_bytes": local_size,
            "expected_full_records": sum(int(m.get("read_count") or 0) for m in metas),
        })

    manifest = {
        "dataset_id": "sanson_brunello",
        "name": "Sanson/Brunello mod-tracr public CRISPR benchmark",
        "bioproject": "PRJNA508200",
        "sra_project": "SRP172473",
        "paper": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6303322/",
        "guide_counter_docs": "https://docs.rs/crate/guide-counter/latest/source/README.md",
        "library": str(library),
        "guide_length": 20,
        "target_start": 20,
        "auto_offset": 20,
        "auto_offset_sample": 100000,
        "offset_mode": "multi",
        "offset_min_fraction": 0.005,
        "mageck_trim5": 0,
        "samples": sample_rows,
        "runs": run_rows,
        "notes": [
            "Guide-counter documents this Sanson/Brunello benchmark using plasmid, RepA, RepB, and RepC samples.",
            "The biological RepA/RepB/RepC FASTQs are concatenations of two SRA runs each to match the documented sample sizes.",
            "DotMatch benchmarks use automatic offset detection because the Sanson protocol searches vector flank CACCG and extracts the following 20 nt guide.",
        ],
    }
    manifest_path = args.out / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    sample_tsv = args.out / "samples.tsv"
    with sample_tsv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample_id", "fastq"], delimiter="\t")
        writer.writeheader()
        for sample in sample_rows:
            writer.writerow({"sample_id": sample["sample_id"], "fastq": sample["fastq"]})

    print(manifest_path)


if __name__ == "__main__":
    main()
