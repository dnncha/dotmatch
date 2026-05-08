#!/usr/bin/env python3
"""Fetch a small public adapter-prefix assignment fixture.

The source fixture comes from the public fast-adapter-trimming repository. The
benchmark uses a fixed R1 window with TruSeq adapter prefixes as known targets;
it is assignment evidence, not adapter trimming evidence.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import http.client
import io
import json
import time
import urllib.error
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "examples" / "oligo_adapter" / "data"
DATASET_ID = "fast_adapter_trimming_truseq_r1"
SOURCE_REPO = "https://github.com/linsalrob/fast-adapter-trimming"
ADAPTER_FASTA_URL = "https://raw.githubusercontent.com/linsalrob/fast-adapter-trimming/main/adapters/truseq.fa"
FASTQ_URL = "https://media.githubusercontent.com/media/linsalrob/fast-adapter-trimming/main/fastq/788707_20180313_S_R1.small.fastq.gz"
FASTQ_NAME = "788707_20180313_S_R1.small.fastq.gz"
TARGET_START = 229
TARGET_LENGTH = 20
USER_AGENT = "DotMatch benchmark fetcher/0.3 (+https://github.com/donncha/dotmatch)"
TRANSIENT_NETWORK_ERRORS = (
    ConnectionError,
    EOFError,
    TimeoutError,
    http.client.IncompleteRead,
    http.client.RemoteDisconnected,
    urllib.error.URLError,
)


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


def fetch_bytes(url: str, timeout: int = 120) -> bytes:
    with urlopen_with_retries(url, timeout=timeout) as resp:
        return resp.read()


def download_file(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    tmp.write_bytes(fetch_bytes(url))
    tmp.replace(dest)


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks).upper()))
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    if not records:
        raise RuntimeError(f"adapter FASTA is empty: {path}")
    return records


def write_targets_from_adapter_fasta(adapter_fasta: Path, targets: Path, target_length: int = TARGET_LENGTH) -> dict[str, object]:
    records = parse_fasta(adapter_fasta)
    seen: dict[str, list[str]] = {}
    for name, seq in records:
        if len(seq) < target_length:
            raise RuntimeError(f"adapter {name} is shorter than target length {target_length}")
        seen.setdefault(seq[:target_length], []).append(name)

    targets.parent.mkdir(parents=True, exist_ok=True)
    with targets.open("w", encoding="utf-8") as out:
        out.write("target_id\ttarget_seq\tgene\n")
        for seq, names in seen.items():
            target_id = names[0]
            out.write(f"{target_id}\t{seq}\t{target_id}\n")

    duplicates = {seq: names for seq, names in seen.items() if len(names) > 1}
    return {
        "adapter_count": len(records),
        "target_count": len(seen),
        "target_length": target_length,
        "duplicate_prefixes": duplicates,
    }


def _source_bytes(source: str | Path) -> bytes:
    if isinstance(source, Path) or not str(source).startswith(("http://", "https://")):
        return Path(source).read_bytes()
    return fetch_bytes(str(source))


def copy_first_fastq_records(source: str | Path, dest: Path, records: int) -> int:
    if records <= 0:
        raise RuntimeError("records must be positive")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    payload = _source_bytes(source)
    written = 0
    with gzip.GzipFile(fileobj=io.BytesIO(payload), mode="rb") as gz, gzip.open(tmp, "wt", encoding="utf-8") as out:
        while written < records:
            header = gz.readline().decode("utf-8")
            if not header:
                break
            seq = gz.readline().decode("utf-8")
            plus = gz.readline().decode("utf-8")
            qual = gz.readline().decode("utf-8")
            if not seq or not plus or not qual or not header.startswith("@") or not plus.startswith("+"):
                raise RuntimeError("remote FASTQ yielded an invalid record")
            out.writelines([header, seq, plus, qual])
            written += 1
    if written < records:
        raise RuntimeError(f"only copied {written} FASTQ records; requested {records}")
    tmp.replace(dest)
    return written


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--records", type=int, default=10000)
    parser.add_argument("--target-start", type=int, default=TARGET_START)
    parser.add_argument("--target-length", type=int, default=TARGET_LENGTH)
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    adapter_fasta = out / "truseq.fa"
    targets = out / "adapter_oligos.tsv"
    source_fastq = out / FASTQ_NAME
    reads = out / f"788707_20180313_S_R1.small.subsample{args.records}.fastq.gz"
    metadata_path = out / "metadata.json"

    download_file(ADAPTER_FASTA_URL, adapter_fasta)
    adapter_metadata = write_targets_from_adapter_fasta(adapter_fasta, targets, target_length=args.target_length)

    written = 0
    if not args.metadata_only:
        download_file(FASTQ_URL, source_fastq)
        written = copy_first_fastq_records(source_fastq, reads, args.records)
    else:
        written = args.records

    metadata = {
        "schema_version": 1,
        "dataset_id": DATASET_ID,
        "source_repo": SOURCE_REPO,
        "source": "linsalrob/fast-adapter-trimming public GitHub repository",
        "license": "MIT",
        "adapter_fasta_url": ADAPTER_FASTA_URL,
        "fastq_url": FASTQ_URL,
        "adapter_fasta": repo_path(adapter_fasta),
        "adapter_fasta_md5": md5_file(adapter_fasta),
        "adapter_fasta_sha256": sha256_file(adapter_fasta),
        "source_fastq": repo_path(source_fastq),
        "source_fastq_md5": md5_file(source_fastq) if source_fastq.exists() else "",
        "source_fastq_sha256": sha256_file(source_fastq) if source_fastq.exists() else "",
        "targets": repo_path(targets),
        "local_fastq": repo_path(reads),
        "requested_records": args.records,
        "written_records": written,
        "target_start": args.target_start,
        "evidence_ready": bool(written and targets.exists() and reads.exists()),
        **adapter_metadata,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(repo_path(metadata_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
