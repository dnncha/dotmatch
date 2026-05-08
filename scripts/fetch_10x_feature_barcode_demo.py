#!/usr/bin/env python3
"""Fetch a small public 10x feature-barcode assignment fixture.

The source dataset is the 10x Genomics 1k PBMC TotalSeq-B 3' v3.1
Cell Surface Protein dataset. The full FASTQ tar is several GB, so this script
uses HTTP byte ranges over the uncompressed tar archive and copies only a
prefix of the antibody R2 FASTQ. The feature reference supplies the fixed
R2 pattern used by the benchmark.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import http.client
import json
import re
import sys
import time
import urllib.error
import urllib.request
import zlib
from pathlib import Path
from typing import NamedTuple


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "examples" / "feature_barcode" / "data"
DATASET_ID = "10x_1k_pbmc_totalseq_b_3p"
DATASET_PAGE = "https://www.10xgenomics.com/datasets/1-k-human-pbm-cs-with-total-seq-b-human-tbnk-antibody-cocktail-3-v-3-1-3-1-standard-6-0-0"
FASTQ_TAR_URL = "https://cf.10xgenomics.com/samples/cell-exp/6.0.0/1k_PBMCs_TotalSeq_B_3p/1k_PBMCs_TotalSeq_B_3p_fastqs.tar"
FASTQ_TAR_BYTES = 7500277760
FEATURE_REF_URL = "https://cf.10xgenomics.com/samples/cell-exp/6.0.0/1k_PBMCs_TotalSeq_B_3p/1k_PBMCs_TotalSeq_B_3p_feature_ref.csv"
FEATURE_REF_MD5 = "b168d2c91af23225cfc040bb4aaf189a"
R2_SUFFIX = "antibody_S5_L001_R2_001.fastq.gz"
USER_AGENT = "DotMatch benchmark fetcher/0.3 (+https://github.com/donncha/dotmatch)"
FASTQ_PREFIX_BYTES = 1024 * 1024
TRANSIENT_NETWORK_ERRORS = (
    ConnectionError,
    EOFError,
    TimeoutError,
    http.client.IncompleteRead,
    http.client.RemoteDisconnected,
    urllib.error.URLError,
)


class TarMember(NamedTuple):
    name: str
    data_offset: int
    size: int


def repo_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def request_with_headers(url_or_request: str | urllib.request.Request) -> urllib.request.Request:
    if isinstance(url_or_request, urllib.request.Request):
        for key, value in {"User-Agent": USER_AGENT, "Accept": "*/*"}.items():
            if not url_or_request.has_header(key):
                url_or_request.add_header(key, value)
        return url_or_request
    return urllib.request.Request(url_or_request, headers={"User-Agent": USER_AGENT, "Accept": "*/*"})


def urlopen_with_retries(
    url_or_request: str | urllib.request.Request,
    timeout: int,
    attempts: int = 3,
    sleep_seconds: float = 1.0,
):
    request = request_with_headers(url_or_request)
    last_error: BaseException | None = None
    for attempt in range(1, attempts + 1):
        try:
            return urllib.request.urlopen(request, timeout=timeout)
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code not in {429, 500, 502, 503, 504} or attempt == attempts:
                raise
        except TRANSIENT_NETWORK_ERRORS as exc:
            last_error = exc
            if attempt == attempts:
                raise
        time.sleep(sleep_seconds * attempt)
    assert last_error is not None
    raise last_error


def fetch_range(url: str, start: int, end: int, timeout: int = 60) -> bytes:
    request = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
    with urlopen_with_retries(request, timeout=timeout) as resp:
        try:
            return resp.read()
        except http.client.IncompleteRead as exc:
            if exc.partial:
                return exc.partial
            raise


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _tar_size(header: bytes) -> int:
    raw = header[124:136].split(b"\0", 1)[0].strip() or b"0"
    return int(raw, 8)


def _tar_name(header: bytes) -> str:
    return header[:100].split(b"\0", 1)[0].decode("utf-8", errors="replace")


def _tar_block_span(size: int) -> int:
    return 512 + ((size + 511) // 512) * 512


def find_tar_member(url: str, suffix: str, tar_bytes: int | None = None) -> TarMember:
    offset = 0
    long_name: str | None = None
    limit = FASTQ_TAR_BYTES if tar_bytes is None else tar_bytes
    while offset < limit:
        header = fetch_range(url, offset, offset + 511)
        if len(header) < 512:
            raise RuntimeError(f"truncated tar header while looking for {suffix}")
        if header == b"\0" * 512:
            break
        name = _tar_name(header)
        size = _tar_size(header)
        data_offset = offset + 512
        if name == "././@LongLink":
            payload = fetch_range(url, data_offset, data_offset + size - 1)
            long_name = payload.split(b"\0", 1)[0].decode("utf-8", errors="replace")
        else:
            actual_name = long_name or name
            long_name = None
            if actual_name.endswith(suffix):
                return TarMember(actual_name, data_offset, size)
        offset += _tar_block_span(size)
    raise RuntimeError(f"could not find tar member ending with {suffix}")


def parse_feature_pattern(pattern: str, sequence: str) -> tuple[int, int]:
    match = re.fullmatch(r"\^([Nn]*)\(BC\)([Nn]*)", pattern)
    if not match:
        raise RuntimeError(f"unsupported fixed-position Feature Barcode pattern: {pattern}")
    return len(match.group(1)), len(sequence)


def download_feature_reference(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    with urlopen_with_retries(url, timeout=60) as resp, tmp.open("wb") as fh:
        fh.write(resp.read())
    tmp.replace(dest)


def write_targets_from_feature_ref(feature_ref: Path, targets: Path) -> dict[str, object]:
    with feature_ref.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise RuntimeError(f"feature reference is empty: {feature_ref}")
    first = rows[0]
    read = first.get("read") or ""
    pattern = first.get("pattern") or ""
    feature_type = first.get("feature_type") or ""
    target_start, target_length = parse_feature_pattern(pattern, first.get("sequence") or "")
    targets.parent.mkdir(parents=True, exist_ok=True)
    with targets.open("w", encoding="utf-8") as out:
        out.write("target_id\ttarget_seq\tgene\n")
        for row in rows:
            if row.get("read") != read or row.get("pattern") != pattern:
                raise RuntimeError("feature reference mixes read/pattern values; only one fixed window is supported")
            seq = (row.get("sequence") or "").strip().upper()
            start, length = parse_feature_pattern(pattern, seq)
            if start != target_start or length != target_length:
                raise RuntimeError("feature reference mixes barcode lengths or starts")
            target_id = (row.get("id") or row.get("name") or f"feature_{len(rows)}").strip()
            gene = (row.get("name") or target_id).strip()
            if not target_id or not seq:
                raise RuntimeError("feature reference contains an empty id or sequence")
            out.write(f"{target_id}\t{seq}\t{gene}\n")
    return {
        "feature_count": len(rows),
        "read": read,
        "pattern": pattern,
        "target_start": target_start,
        "target_length": target_length,
        "feature_type": feature_type,
    }


def decompress_concatenated_gzip_prefix(data: bytes) -> bytes:
    remaining = data
    chunks: list[bytes] = []
    while remaining:
        decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
        try:
            chunks.append(decompressor.decompress(remaining))
        except zlib.error:
            break
        if decompressor.unused_data:
            remaining = decompressor.unused_data
            continue
        break
    return b"".join(chunks)


def copy_first_tarred_fastq_records(
    tar_url: str,
    suffix: str,
    dest: Path,
    records: int,
    prefix_bytes: int = FASTQ_PREFIX_BYTES,
    tar_bytes: int | None = None,
) -> int:
    member = find_tar_member(tar_url, suffix) if tar_bytes is None else find_tar_member(tar_url, suffix, tar_bytes=tar_bytes)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    compressed = fetch_range(
        tar_url,
        member.data_offset,
        min(member.data_offset + prefix_bytes - 1, member.data_offset + member.size - 1),
        timeout=120,
    )
    text = decompress_concatenated_gzip_prefix(compressed).decode("utf-8", errors="replace")
    lines = text.splitlines(keepends=True)
    written = 0
    with gzip.open(tmp, "wt", encoding="utf-8") as out:
        for i in range(0, len(lines) - 3, 4):
            if written >= records:
                break
            record = lines[i : i + 4]
            if not record[0].startswith("@") or not record[2].startswith("+"):
                raise RuntimeError("remote FASTQ yielded an invalid record")
            out.writelines(record)
            written += 1
    if written < records:
        raise RuntimeError(f"only copied {written} FASTQ records from prefix; requested {records}")
    tmp.replace(dest)
    return written


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--records", type=int, default=20000)
    parser.add_argument("--prefix-bytes", type=int, default=FASTQ_PREFIX_BYTES)
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    feature_ref = out / "1k_PBMCs_TotalSeq_B_3p_feature_ref.csv"
    targets = out / "feature_barcodes.tsv"
    reads = out / f"1k_PBMCs_TotalSeq_B_3p_antibody_S5_L001_R2.subsample{args.records}.fastq.gz"
    metadata_path = out / "metadata.json"

    download_feature_reference(FEATURE_REF_URL, feature_ref)
    ref_md5 = md5_file(feature_ref)
    if ref_md5 != FEATURE_REF_MD5:
        raise RuntimeError(f"feature reference md5 mismatch: expected {FEATURE_REF_MD5}, observed {ref_md5}")
    feature_metadata = write_targets_from_feature_ref(feature_ref, targets)
    member = find_tar_member(FASTQ_TAR_URL, R2_SUFFIX)
    written = 0
    if not args.metadata_only:
        written = copy_first_tarred_fastq_records(
            FASTQ_TAR_URL,
            R2_SUFFIX,
            reads,
            args.records,
            prefix_bytes=args.prefix_bytes,
        )
    else:
        written = args.records

    metadata = {
        "schema_version": 1,
        "dataset_id": DATASET_ID,
        "dataset_page": DATASET_PAGE,
        "license": "CC BY 4.0",
        "source": "10x Genomics public datasets",
        "fastq_tar_url": FASTQ_TAR_URL,
        "fastq_tar_bytes": FASTQ_TAR_BYTES,
        "feature_ref_url": FEATURE_REF_URL,
        "feature_ref_md5": FEATURE_REF_MD5,
        "feature_ref": repo_path(feature_ref),
        "targets": repo_path(targets),
        "r2_entry": member.name,
        "r2_entry_data_offset": member.data_offset,
        "r2_entry_compressed_bytes": member.size,
        "local_fastq": repo_path(reads),
        "requested_records": args.records,
        "written_records": written,
        "compressed_prefix_bytes": args.prefix_bytes,
        "evidence_ready": bool(written and targets.exists() and feature_ref.exists()),
        **feature_metadata,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(repo_path(metadata_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
