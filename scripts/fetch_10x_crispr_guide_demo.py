#!/usr/bin/env python3
"""Fetch a small public 10x CRISPR Guide Capture assignment fixture.

The source dataset is the 10x Genomics 1k A375 GEM-X Single Cell 5' CRISPR
Screening dataset. The full FASTQ tar is several GB, so this script uses HTTP
byte ranges over the uncompressed tar archive and copies only a prefix of the
CRISPR R2 FASTQ. The feature reference contains multiple CRISPR guide patterns;
this fixture selects the observed fixed-window group with the most exact
assignments in the copied R2 prefix.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from fetch_10x_feature_barcode_demo import (  # noqa: E402
    FASTQ_PREFIX_BYTES,
    copy_first_tarred_fastq_records,
    download_feature_reference,
    find_tar_member,
    md5_file,
)


DEFAULT_OUT = ROOT / "examples" / "perturb_seq" / "data"
DATASET_ID = "10x_1k_a375_gemx_5p_crispr"
DATASET_PAGE = "https://www.10xgenomics.com/datasets/1k-CRISPR-5p-gemx"
FASTQ_TAR_URL = "https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/1k_CRISPR_5p_gemx_Multiplex/1k_CRISPR_5p_gemx_Multiplex_fastqs.tar"
FASTQ_TAR_BYTES = 4961945600
FEATURE_REF_URL = "https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/1k_CRISPR_5p_gemx_Multiplex/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv"
FEATURE_REF_MD5 = "89088af5191f421a3a21d59fa60d39ed"
R2_SUFFIX = "crispr_S1_L001_R2_001.fastq.gz"


def repo_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def parse_crispr_pattern(pattern: str) -> str:
    cleaned = pattern.strip()
    if cleaned.startswith("^"):
        cleaned = cleaned[1:]
    match = re.fullmatch(r"([ACGTNacgtn]+)\(BC\)", cleaned)
    if not match:
        raise RuntimeError(f"unsupported CRISPR Guide Capture pattern: {pattern}")
    return match.group(1).upper()


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


def read_crispr_rows(feature_ref: Path) -> list[dict[str, str]]:
    with feature_ref.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    crispr_rows = []
    for row in rows:
        if row.get("feature_type") != "CRISPR Guide Capture":
            continue
        seq = (row.get("sequence") or "").strip().upper()
        if not seq:
            continue
        parse_crispr_pattern(row.get("pattern") or "")
        row = dict(row)
        row["sequence"] = seq
        crispr_rows.append(row)
    if not crispr_rows:
        raise RuntimeError(f"no CRISPR Guide Capture rows found in {feature_ref}")
    return crispr_rows


def grouped_crispr_rows(rows: list[dict[str, str]]) -> list[list[dict[str, str]]]:
    groups: dict[tuple[str, int], list[dict[str, str]]] = collections.defaultdict(list)
    for row in rows:
        seq = row["sequence"]
        if len(seq) < 12:
            continue
        groups[(row["pattern"], len(seq))].append(row)
    if not groups:
        raise RuntimeError("no fixed-length CRISPR Guide Capture groups found")
    return list(groups.values())


def score_group(rows: list[dict[str, str]], sequences: list[str]) -> dict[str, object]:
    pattern = rows[0]["pattern"]
    anchor = parse_crispr_pattern(pattern)
    target_length = len(rows[0]["sequence"])
    targets = {row["sequence"] for row in rows}
    starts: collections.Counter[int] = collections.Counter()
    for seq in sequences:
        anchor_start = seq.find(anchor)
        if anchor_start < 0:
            continue
        target_start = anchor_start + len(anchor)
        starts[target_start] += 1
    target_start, anchor_observed = starts.most_common(1)[0] if starts else (0, 0)
    exact_assigned = 0
    if starts:
        for seq in sequences:
            observed = seq[target_start : target_start + target_length]
            if len(observed) == target_length and observed in targets:
                exact_assigned += 1
    return {
        "pattern": pattern,
        "selected_pattern": pattern,
        "anchor": anchor,
        "target_start": target_start,
        "target_length": target_length,
        "feature_count": len(rows),
        "anchor_observed_reads": anchor_observed,
        "exact_assigned_in_subsample": exact_assigned,
        "total_reads_scanned": len(sequences),
        "rows": rows,
    }


def select_observed_group(rows: list[dict[str, str]], sequences: list[str]) -> dict[str, object]:
    scored = [score_group(group, sequences) for group in grouped_crispr_rows(rows)]
    scored.sort(
        key=lambda item: (
            int(item["exact_assigned_in_subsample"]),
            int(item["anchor_observed_reads"]),
            int(item["feature_count"]),
            int(item["target_length"]),
        ),
        reverse=True,
    )
    selected = scored[0]
    if int(selected["exact_assigned_in_subsample"]) <= 0:
        raise RuntimeError("no CRISPR guide group produced exact assignments in the FASTQ prefix")
    selected["considered_groups"] = [
        {
            key: value
            for key, value in item.items()
            if key
            in {
                "pattern",
                "target_start",
                "target_length",
                "feature_count",
                "anchor_observed_reads",
                "exact_assigned_in_subsample",
                "total_reads_scanned",
            }
        }
        for item in scored
    ]
    return selected


def write_observed_targets_from_feature_ref(feature_ref: Path, reads: Path, targets: Path) -> dict[str, object]:
    rows = read_crispr_rows(feature_ref)
    sequences = read_fastq_sequences(reads)
    selection = select_observed_group(rows, sequences)
    selected_rows = selection["rows"]
    targets.parent.mkdir(parents=True, exist_ok=True)
    with targets.open("w", encoding="utf-8") as out:
        out.write("target_id\ttarget_seq\tgene\n")
        for row in selected_rows:
            target_id = (row.get("id") or row.get("name") or "").strip()
            gene = (row.get("target_gene_name") or row.get("name") or target_id).strip()
            out.write(f"{target_id}\t{row['sequence']}\t{gene}\n")
    return {
        key: value
        for key, value in selection.items()
        if key != "rows"
    } | {
        "dropped_feature_count": len(rows) - len(selected_rows),
        "read": selected_rows[0].get("read", ""),
        "feature_type": "CRISPR Guide Capture",
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--records", type=int, default=20000)
    parser.add_argument("--prefix-bytes", type=int, default=FASTQ_PREFIX_BYTES)
    args = parser.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    feature_ref = out / "1k_CRISPR_5p_gemx_count_feature_reference.csv"
    targets = out / "crispr_guides.tsv"
    reads = out / f"1k_CRISPR_5p_gemx_crispr_S1_L001_R2.subsample{args.records}.fastq.gz"
    metadata_path = out / "metadata.json"

    download_feature_reference(FEATURE_REF_URL, feature_ref)
    ref_md5 = md5_file(feature_ref)
    if ref_md5 != FEATURE_REF_MD5:
        raise RuntimeError(f"feature reference md5 mismatch: expected {FEATURE_REF_MD5}, observed {ref_md5}")
    member = find_tar_member(FASTQ_TAR_URL, R2_SUFFIX, tar_bytes=FASTQ_TAR_BYTES)
    written = copy_first_tarred_fastq_records(
        FASTQ_TAR_URL,
        R2_SUFFIX,
        reads,
        args.records,
        prefix_bytes=args.prefix_bytes,
        tar_bytes=FASTQ_TAR_BYTES,
    )
    feature_metadata = write_observed_targets_from_feature_ref(feature_ref, reads, targets)
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
