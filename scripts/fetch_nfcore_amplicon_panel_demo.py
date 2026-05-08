#!/usr/bin/env python3
"""Fetch a small public nf-core viralrecon ARTIC amplicon primer fixture.

The source data are nf-core/test-datasets viralrecon Illumina amplicon FASTQs
and the ARTIC V3 SARS-CoV-2 primer FASTA. The fixture selects the full-primer
length group with the most exact R1 prefix assignments and writes that fixed
length target set for DotMatch's count command.
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json
import shutil
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "examples" / "amplicon_panel" / "data"
DATASET_ID = "nfcore_viralrecon_artic_v3_amplicon_sample1"
DATASET_PAGE = "https://github.com/nf-core/test-datasets/tree/viralrecon/illumina/amplicon"
FASTQ_URL = "https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz"
PRIMER_FASTA_URL = "https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/genome/MN908947.3/amplicon/nCoV-2019.artic.V3.primer.fasta"
USER_AGENT = "DotMatch benchmark fetcher/0.3 (+https://github.com/donncha/dotmatch)"


def repo_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def download_url(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "*/*"})
    with urllib.request.urlopen(request, timeout=60) as resp, tmp.open("wb") as out:
        shutil.copyfileobj(resp, out)
    tmp.replace(dest)


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def copy_first_fastq_records(src: Path, dest: Path, records: int) -> int:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    written = 0
    with open_text(src) as inp, gzip.open(tmp, "wt", encoding="utf-8") as out:
        while written < records:
            header = inp.readline()
            if not header:
                break
            seq = inp.readline()
            plus = inp.readline()
            qual = inp.readline()
            if not seq or not plus or not qual:
                raise RuntimeError(f"truncated FASTQ record in {src}")
            out.writelines([header, seq, plus, qual])
            written += 1
    if written <= 0:
        raise RuntimeError(f"no FASTQ records copied from {src}")
    tmp.replace(dest)
    return written


def parse_primer_fasta(path: Path) -> list[tuple[str, str]]:
    primers: list[tuple[str, str]] = []
    name: str | None = None
    seq_parts: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                primers.append((name, "".join(seq_parts).upper()))
            name = line[1:].split()[0]
            seq_parts = []
        else:
            seq_parts.append(line)
    if name is not None:
        primers.append((name, "".join(seq_parts).upper()))
    primers = [(name, seq) for name, seq in primers if seq]
    if not primers:
        raise RuntimeError(f"no primer sequences found in {path}")
    return primers


def read_fastq_sequences(path: Path) -> list[str]:
    sequences: list[str] = []
    with open_text(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip().upper()
            plus = fh.readline()
            qual = fh.readline()
            if not seq or not plus or not qual:
                raise RuntimeError(f"truncated FASTQ record in {path}")
            sequences.append(seq)
    return sequences


def gene_name(primer_id: str) -> str:
    parts = primer_id.rsplit("_", 1)
    return parts[0] if len(parts) == 2 and parts[1] in {"LEFT", "RIGHT"} else primer_id


def select_length_group(primers: list[tuple[str, str]], sequences: list[str]) -> dict[str, object]:
    groups: dict[int, list[tuple[str, str]]] = collections.defaultdict(list)
    for name, seq in primers:
        groups[len(seq)].append((name, seq))
    scored = []
    for length, rows in groups.items():
        target_by_seq = {seq: name for name, seq in rows}
        exact = sum(1 for seq in sequences if seq[:length] in target_by_seq)
        scored.append(
            {
                "target_length": length,
                "feature_count": len(rows),
                "unique_target_sequences": len(target_by_seq),
                "exact_assigned_in_subsample": exact,
                "total_reads_scanned": len(sequences),
                "rows": rows,
            }
        )
    scored.sort(
        key=lambda item: (
            int(item["exact_assigned_in_subsample"]),
            int(item["unique_target_sequences"]),
            int(item["target_length"]),
        ),
        reverse=True,
    )
    selected = scored[0]
    if int(selected["exact_assigned_in_subsample"]) <= 0:
        raise RuntimeError("no ARTIC primer length group produced exact R1 prefix assignments")
    selected["considered_groups"] = [
        {
            key: value
            for key, value in item.items()
            if key
            in {
                "target_length",
                "feature_count",
                "unique_target_sequences",
                "exact_assigned_in_subsample",
                "total_reads_scanned",
            }
        }
        for item in scored
    ]
    return selected


def write_observed_targets_from_primer_fasta(primer_fasta: Path, reads: Path, targets: Path) -> dict[str, object]:
    primers = parse_primer_fasta(primer_fasta)
    sequences = read_fastq_sequences(reads)
    selected = select_length_group(primers, sequences)
    rows = selected["rows"]
    targets.parent.mkdir(parents=True, exist_ok=True)
    with targets.open("w", encoding="utf-8") as out:
        out.write("target_id\ttarget_seq\tgene\n")
        for primer_id, seq in rows:
            out.write(f"{primer_id}\t{seq}\t{gene_name(primer_id)}\n")
    return {
        key: value
        for key, value in selected.items()
        if key != "rows"
    } | {
        "selected_target_length": selected["target_length"],
        "target_start": 0,
        "dropped_feature_count": len(primers) - len(rows),
        "feature_type": "ARTIC V3 primer start",
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--records", type=int, default=20000)
    args = parser.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    source_fastq = out / "nfcore_viralrecon_sample1_R1.fastq.gz"
    reads = out / f"nfcore_viralrecon_sample1_R1.subsample{args.records}.fastq.gz"
    primer_fasta = out / "nCoV-2019.artic.V3.primer.fasta"
    metadata_path = out / "metadata.json"

    download_url(FASTQ_URL, source_fastq)
    download_url(PRIMER_FASTA_URL, primer_fasta)
    written = copy_first_fastq_records(source_fastq, reads, args.records)
    temporary_targets = out / "artic_v3_primers.selected.tsv"
    primer_metadata = write_observed_targets_from_primer_fasta(primer_fasta, reads, temporary_targets)
    targets = out / f"artic_v3_primers_len{primer_metadata['target_length']}.tsv"
    temporary_targets.replace(targets)

    metadata = {
        "schema_version": 1,
        "dataset_id": DATASET_ID,
        "dataset_page": DATASET_PAGE,
        "license": "See nf-core/test-datasets repository license for fixture reuse terms.",
        "source": "nf-core/test-datasets viralrecon public test data",
        "fastq_url": FASTQ_URL,
        "primer_fasta_url": PRIMER_FASTA_URL,
        "source_fastq": repo_path(source_fastq),
        "local_fastq": repo_path(reads),
        "primer_fasta": repo_path(primer_fasta),
        "targets": repo_path(targets),
        "requested_records": args.records,
        "written_records": written,
        "target_start": 0,
        "target_length": primer_metadata["target_length"],
        "evidence_ready": bool(written and targets.exists() and reads.exists() and primer_fasta.exists()),
        **primer_metadata,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(repo_path(metadata_path))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
