#!/usr/bin/env python3
"""Run the bounded, accessioned GSE146194 direct-guide-capture case study.

The public workflow streams a fixed prefix of SRR11214031, extracts the 32
published UPR guide-barcode targets from the paper supplement, chooses a fixed
window on a discovery prefix, and evaluates held-out reads with DotMatch and
independent exact/Hamming-radius oracles.  It does not perform cell/UMI
quantification or perturbation-effect analysis.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import json
import os
import platform
import re
import resource
import shlex
import subprocess
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterable, Iterator, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROTOCOL = ROOT / "examples" / "perturb_seq_gse146194" / "protocol.json"
DEFAULT_FIXTURE = ROOT / "examples" / "perturb_seq_gse146194" / "fixture"
DEFAULT_WORK = ROOT / "examples" / "perturb_seq_gse146194" / "work"
DEFAULT_RESULTS = ROOT / "examples" / "perturb_seq_gse146194" / "results.json"
DEFAULT_PROVENANCE = ROOT / "examples" / "perturb_seq_gse146194" / "provenance.json"
DEFAULT_CSV = ROOT / "benchmarks" / "raw" / "perturb_seq_gse146194.csv"
DEFAULT_REPORT = ROOT / "docs" / "benchmarks" / "perturb_seq_gse146194" / "README.md"
DEFAULT_HTML = ROOT / "docs" / "benchmarks" / "perturb_seq_gse146194" / "report.html"
ASSIGNMENT_FIELDS = [
    "read_id",
    "observed_seq",
    "target_id",
    "target_seq",
    "distance",
    "status",
    "match_count",
    "second_best_distance",
    "correction",
]


def access_and_reuse_record() -> dict[str, object]:
    """Return the source-access record without asserting rights we do not grant."""
    return {
        "sequence_data": {
            "access": "unrestricted public NCBI SRA record mirrored by ENA; no credentials required",
            "policy": "NCBI places no restrictions on molecular-data use or distribution but cannot transfer any rights claimed by submitters",
            "policy_url": "https://www.ncbi.nlm.nih.gov/home/about/policies/",
        },
        "publisher_supplement": {
            "access": "direct publisher download linked to the primary article",
            "license_status": "source terms apply; this workflow asserts no redistribution license",
            "policy_url": "https://support.springernature.com/en/support/solutions/articles/6000210902-supplementary-information",
        },
        "repository_redistribution": {
            "raw_reads": False,
            "publisher_workbook": False,
            "committed_material": "aggregate results, hashes, protocol, provenance, and a synthetic contract fixture",
        },
    }


@dataclass(frozen=True)
class Target:
    target_id: str
    sequence: str
    gene: str = ""


@dataclass(frozen=True)
class FastqRecord:
    read_id: str
    header: str
    sequence: str
    plus: str
    quality: str

    def canonical_bytes(self) -> bytes:
        return (
            f"{self.header}\n{self.sequence}\n{self.plus}\n{self.quality}\n".encode(
                "ascii"
            )
        )


class CountingReader:
    def __init__(self, source: BinaryIO, maximum_bytes: int) -> None:
        self.source = source
        self.maximum_bytes = maximum_bytes
        self.bytes_read = 0

    def read(self, size: int = -1) -> bytes:
        data = self.source.read(size)
        self.bytes_read += len(data)
        if self.bytes_read > self.maximum_bytes:
            raise RuntimeError(
                f"compressed stream exceeded the {self.maximum_bytes}-byte retrieval cap"
            )
        return data

    def readable(self) -> bool:
        return True


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def public_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def public_command(command: Sequence[str]) -> str:
    displayed: list[str] = []
    for index, argument in enumerate(command):
        path = Path(argument)
        if path.is_absolute():
            try:
                relative = path.resolve().relative_to(ROOT)
            except ValueError:
                displayed.append(argument)
            else:
                prefix = "./" if index == 0 else ""
                displayed.append(prefix + relative.as_posix())
        else:
            displayed.append(argument)
    return shlex.join(displayed)


def stable_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"expected a JSON object in {path}")
    return value


def request(url: str, timeout: int):
    return urllib.request.urlopen(
        urllib.request.Request(
            url,
            headers={
                "User-Agent": "DotMatch-GSE146194-case-study/1.0 (+https://github.com/dnncha/dotmatch)"
            },
        ),
        timeout=timeout,
    )


def download_verified(
    url: str,
    destination: Path,
    expected_sha256: str,
    maximum_bytes: int,
    timeout: int,
) -> dict[str, object]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file() and sha256_path(destination) == expected_sha256:
        return {
            "url": url,
            "path": public_path(destination),
            "bytes": destination.stat().st_size,
            "sha256": expected_sha256,
            "cache": "verified_reuse",
        }
    temporary = destination.with_suffix(destination.suffix + ".partial")
    digest = hashlib.sha256()
    downloaded = 0
    with request(url, timeout) as response, temporary.open("wb") as output:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            downloaded += len(chunk)
            if downloaded > maximum_bytes:
                raise RuntimeError(
                    f"download exceeded the {maximum_bytes}-byte cap: {url}"
                )
            digest.update(chunk)
            output.write(chunk)
    observed = digest.hexdigest()
    if observed != expected_sha256:
        raise RuntimeError(
            f"SHA-256 mismatch for {url}: expected {expected_sha256}, observed {observed}"
        )
    temporary.replace(destination)
    return {
        "url": url,
        "path": public_path(destination),
        "bytes": downloaded,
        "sha256": observed,
        "cache": "downloaded",
    }


def excel_column_index(reference: str) -> int:
    match = re.match(r"([A-Z]+)", reference.upper())
    if match is None:
        raise RuntimeError(f"invalid XLSX cell reference: {reference}")
    value = 0
    for char in match.group(1):
        value = value * 26 + ord(char) - ord("A") + 1
    return value - 1


def shared_strings(archive: zipfile.ZipFile) -> list[str]:
    try:
        data = archive.read("xl/sharedStrings.xml")
    except KeyError:
        return []
    root = ET.fromstring(data)
    values: list[str] = []
    for item in root.findall("{*}si"):
        values.append("".join(node.text or "" for node in item.findall(".//{*}t")))
    return values


def worksheet_path(archive: zipfile.ZipFile, wanted: str) -> str:
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    targets = {
        item.attrib["Id"]: item.attrib["Target"]
        for item in relationships.findall("{*}Relationship")
    }
    relationship_key = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
    for sheet in workbook.findall(".//{*}sheet"):
        if sheet.attrib.get("name") != wanted:
            continue
        target = targets[sheet.attrib[relationship_key]].lstrip("/")
        return target if target.startswith("xl/") else f"xl/{target}"
    raise RuntimeError(f"worksheet not found: {wanted}")


def worksheet_rows(path: Path, wanted: str) -> list[list[str]]:
    with zipfile.ZipFile(path) as archive:
        strings = shared_strings(archive)
        sheet = ET.fromstring(archive.read(worksheet_path(archive, wanted)))
    rows: list[list[str]] = []
    for row in sheet.findall(".//{*}row"):
        values: dict[int, str] = {}
        for cell in row.findall("{*}c"):
            reference = cell.attrib.get("r", "")
            column = excel_column_index(reference)
            cell_type = cell.attrib.get("t", "")
            if cell_type == "inlineStr":
                value = "".join(node.text or "" for node in cell.findall(".//{*}t"))
            else:
                node = cell.find("{*}v")
                value = "" if node is None else (node.text or "")
                if cell_type == "s" and value:
                    value = strings[int(value)]
            values[column] = value.strip()
        if values:
            width = max(values) + 1
            rows.append([values.get(index, "") for index in range(width)])
    return rows


def extract_targets(workbook: Path, protocol: dict, output: Path) -> list[Target]:
    config = protocol["inputs"]["guide_library"]
    rows = worksheet_rows(workbook, str(config["worksheet"]))
    id_name = str(config["id_column"])
    sequence_name = str(config["sequence_column"])
    header_index = -1
    id_column = sequence_column = gene_column = -1
    for index, row in enumerate(rows):
        if id_name in row and sequence_name in row:
            header_index = index
            id_column = row.index(id_name)
            sequence_column = row.index(sequence_name)
            gene_column = row.index("Target gene name") if "Target gene name" in row else -1
            break
    if header_index < 0:
        raise RuntimeError("published guide-library columns were not found in the workbook")
    target_length = int(config["target_length"])
    targets: list[Target] = []
    for row in rows[header_index + 1 :]:
        if max(id_column, sequence_column) >= len(row):
            continue
        target_id = row[id_column].strip()
        sequence = row[sequence_column].strip().upper()
        if not target_id and not sequence:
            continue
        if not target_id or not re.fullmatch(r"[ACGT]+", sequence):
            raise RuntimeError(f"invalid guide-library row: {row}")
        if len(sequence) != target_length:
            raise RuntimeError(
                f"guide {target_id} has length {len(sequence)}, expected {target_length}"
            )
        gene = row[gene_column].strip() if 0 <= gene_column < len(row) else ""
        targets.append(Target(target_id, sequence, gene))
    expected = int(protocol["dataset"]["guide_count_expected"])
    if len(targets) != expected:
        raise RuntimeError(f"expected {expected} guides, extracted {len(targets)}")
    if len({target.target_id for target in targets}) != len(targets):
        raise RuntimeError("guide identifiers are not unique")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["target_id", "target_seq", "gene"])
        for target in targets:
            writer.writerow([target.target_id, target.sequence, target.gene])
    return targets


def parse_fastq_record(handle: BinaryIO) -> FastqRecord | None:
    lines = [handle.readline() for _ in range(4)]
    if not lines[0]:
        return None
    if not all(lines):
        raise RuntimeError("truncated FASTQ record in streamed prefix")
    try:
        header, sequence, plus, quality = [line.rstrip(b"\r\n").decode("ascii") for line in lines]
    except UnicodeDecodeError as exc:
        raise RuntimeError("FASTQ stream is not ASCII") from exc
    if not header.startswith("@") or not plus.startswith("+"):
        raise RuntimeError("invalid FASTQ record in streamed prefix")
    sequence = sequence.upper()
    if len(sequence) != len(quality):
        raise RuntimeError("FASTQ sequence and quality lengths differ")
    return FastqRecord(header[1:].split()[0], header, sequence, plus, quality)


def deterministic_gzip_writer(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return raw, compressed


def stream_prefix(protocol: dict, evaluation_path: Path) -> tuple[list[FastqRecord], dict]:
    config = protocol["inputs"]["fastq"]
    plan = protocol["analysis_plan"]
    total_records = int(plan["total_prefix_records"])
    discovery_count = int(plan["discovery_records"])
    timeout = int(protocol["resource_bounds"]["network_timeout_seconds"])
    network_cap = int(protocol["resource_bounds"]["maximum_local_work_bytes"])
    discovery: list[FastqRecord] = []
    prefix_digest = hashlib.sha256()
    evaluation_digest = hashlib.sha256()
    response = request(str(config["url"]), timeout)
    counter = CountingReader(response, network_cap)
    compressed = gzip.GzipFile(fileobj=counter, mode="rb")
    raw_output, gzip_output = deterministic_gzip_writer(evaluation_path.with_suffix(".partial"))
    observed = 0
    try:
        while observed < total_records:
            record = parse_fastq_record(compressed)
            if record is None:
                break
            canonical = record.canonical_bytes()
            prefix_digest.update(canonical)
            if observed < discovery_count:
                discovery.append(record)
            else:
                gzip_output.write(canonical)
                evaluation_digest.update(canonical)
            observed += 1
    finally:
        gzip_output.close()
        raw_output.close()
        compressed.close()
        response.close()
    if observed != total_records:
        raise RuntimeError(f"requested {total_records} FASTQ records but received {observed}")
    temporary = evaluation_path.with_suffix(".partial")
    temporary.replace(evaluation_path)
    return discovery, {
        "source_url": str(config["url"]),
        "selected_records": observed,
        "discovery_records": discovery_count,
        "evaluation_records": observed - discovery_count,
        "compressed_network_bytes": counter.bytes_read,
        "selected_prefix_uncompressed_sha256": prefix_digest.hexdigest(),
        "evaluation_uncompressed_sha256": evaluation_digest.hexdigest(),
        "evaluation_gzip_sha256": sha256_path(evaluation_path),
        "full_archive_md5_status": "registry_value_recorded_not_reverified",
        "full_archive_reported_md5": str(config["archive_reported_md5"]),
    }


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def discover_window(
    records: Sequence[FastqRecord], targets: Sequence[Target], protocol: dict
) -> tuple[list[Target], dict]:
    config = protocol["analysis_plan"]["window_discovery"]
    length = int(config["sequence_length"])
    if not records:
        raise RuntimeError("window discovery received no records")
    maximum_start = min(len(record.sequence) for record in records) - length
    if maximum_start < 0:
        raise RuntimeError("discovery reads are shorter than the guide-barcode length")
    candidates: list[tuple[int, int, int, str, int, list[Target]]] = []
    for orientation_index, orientation in enumerate(config["orientations"]):
        oriented = [
            Target(
                target.target_id,
                target.sequence
                if orientation == "forward"
                else reverse_complement(target.sequence),
                target.gene,
            )
            for target in targets
        ]
        by_sequence: dict[str, list[str]] = {}
        for target in oriented:
            by_sequence.setdefault(target.sequence, []).append(target.target_id)
        for start in range(maximum_start + 1):
            assigned = 0
            observed_targets: set[str] = set()
            for record in records:
                window = record.sequence[start : start + length]
                hits = by_sequence.get(window, [])
                if len(hits) == 1:
                    assigned += 1
                    observed_targets.add(hits[0])
            candidates.append(
                (
                    assigned,
                    len(observed_targets),
                    -orientation_index,
                    str(orientation),
                    start,
                    oriented,
                )
            )
    assigned, distinct, _orientation_order, orientation, start, oriented_targets = max(
        candidates, key=lambda row: (row[0], row[1], row[2], -row[4])
    )
    fraction = assigned / len(records)
    if distinct < int(config["minimum_distinct_exact_targets"]):
        raise RuntimeError(
            f"window discovery observed {distinct} distinct exact targets; "
            f"minimum is {config['minimum_distinct_exact_targets']}"
        )
    if fraction < float(config["minimum_exact_assignment_fraction"]):
        raise RuntimeError(
            f"window discovery exact-assignment fraction {fraction:.6f} is below "
            f"{config['minimum_exact_assignment_fraction']}"
        )
    ranked = sorted(candidates, key=lambda row: (row[0], row[1]), reverse=True)[:5]
    return oriented_targets, {
        "orientation": orientation,
        "target_start": start,
        "target_length": length,
        "records": len(records),
        "exact_unique_assignments": assigned,
        "exact_assignment_fraction": fraction,
        "distinct_exact_targets": distinct,
        "top_candidates": [
            {
                "orientation": row[3],
                "target_start": row[4],
                "exact_unique_assignments": row[0],
                "distinct_exact_targets": row[1],
            }
            for row in ranked
        ],
    }


def write_targets(path: Path, targets: Sequence[Target]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["target_id", "target_seq", "gene"])
        for target in targets:
            writer.writerow([target.target_id, target.sequence, target.gene])


def iter_fastq(path: Path) -> Iterator[FastqRecord]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as handle:
        while True:
            record = parse_fastq_record(handle)
            if record is None:
                return
            yield record


def hamming(left: str, right: str) -> int:
    if len(left) != len(right):
        return max(len(left), len(right))
    return sum(a != b for a, b in zip(left, right))


def oracle_row(
    record: FastqRecord, targets: Sequence[Target], start: int, length: int, k: int
) -> dict[str, str]:
    if start + length > len(record.sequence):
        return {
            "read_id": record.read_id,
            "observed_seq": "",
            "target_id": "",
            "target_seq": "",
            "distance": "-1",
            "status": "invalid",
            "match_count": "0",
            "second_best_distance": "-1",
            "correction": "invalid",
        }
    observed = record.sequence[start : start + length]
    compatible = sorted(
        (
            (hamming(observed, target.sequence), index, target)
            for index, target in enumerate(targets)
            if hamming(observed, target.sequence) <= k
        ),
        key=lambda value: (value[0], value[1]),
    )
    if not compatible:
        return {
            "read_id": record.read_id,
            "observed_seq": observed,
            "target_id": "",
            "target_seq": "",
            "distance": "-1",
            "status": "none",
            "match_count": "0",
            "second_best_distance": "-1",
            "correction": "none",
        }
    best_distance, _index, best = compatible[0]
    status = "unique" if len(compatible) == 1 else "ambiguous"
    correction = "exact" if best_distance == 0 else "substitution"
    if status == "ambiguous":
        correction = "ambiguous"
    return {
        "read_id": record.read_id,
        "observed_seq": observed,
        "target_id": best.target_id,
        "target_seq": best.sequence,
        "distance": str(best_distance),
        "status": status,
        "match_count": str(len(compatible)),
        "second_best_distance": str(compatible[1][0] if len(compatible) > 1 else -1),
        "correction": correction,
    }


def write_oracle(
    reads: Path,
    targets: Sequence[Target],
    start: int,
    length: int,
    k: int,
    output: Path,
) -> tuple[dict, float]:
    begin = time.perf_counter()
    rows = [oracle_row(record, targets, start, length, k) for record in iter_fastq(reads)]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=ASSIGNMENT_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    return summarize_assignments(rows), time.perf_counter() - begin


def read_assignment_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        if "distance" not in row and "best_distance" in row:
            row["distance"] = row["best_distance"]
    return rows


def summarize_assignments(rows: Iterable[dict[str, str]]) -> dict:
    rows = list(rows)
    counts = Counter(row["status"] for row in rows)
    assigned = [row for row in rows if row["status"] == "unique"]
    exact = sum(int(row["distance"]) == 0 for row in assigned)
    corrected = sum(int(row["distance"]) > 0 for row in assigned)
    per_target = Counter(row["target_id"] for row in assigned)
    unmatched_windows = Counter(
        row["observed_seq"] for row in rows if row["status"] == "none" and row["observed_seq"]
    )
    total = len(rows)
    return {
        "total_reads": total,
        "assigned_unique": counts["unique"],
        "assigned_exact": exact,
        "assigned_corrected": corrected,
        "ambiguous": counts["ambiguous"],
        "unmatched": counts["none"],
        "invalid": counts["invalid"],
        "assignment_rate": counts["unique"] / total if total else 0.0,
        "ambiguous_rate": counts["ambiguous"] / total if total else 0.0,
        "unmatched_rate": counts["none"] / total if total else 0.0,
        "invalid_rate": counts["invalid"] / total if total else 0.0,
        "distinct_assigned_targets": len(per_target),
        "per_target_unique_counts": dict(sorted(per_target.items())),
        "top_unmatched_windows": [
            {"sequence": sequence, "count": count}
            for sequence, count in unmatched_windows.most_common(10)
        ],
    }


def compare_assignments(dotmatch: Path, oracle: Path) -> dict:
    left = read_assignment_rows(dotmatch)
    right = read_assignment_rows(oracle)
    fields = ["read_id", "status", "target_id", "distance"]
    mismatches: list[dict[str, object]] = []
    for index in range(max(len(left), len(right))):
        if index >= len(left) or index >= len(right):
            mismatches.append({"row": index + 1, "reason": "row_count"})
            continue
        differences = {
            field: {"dotmatch": left[index].get(field, ""), "oracle": right[index].get(field, "")}
            for field in fields
            if left[index].get(field, "") != right[index].get(field, "")
        }
        if differences and len(mismatches) < 20:
            mismatches.append({"row": index + 1, "differences": differences})
    mismatch_count = sum(
        1
        for index in range(min(len(left), len(right)))
        if any(left[index].get(field, "") != right[index].get(field, "") for field in fields)
    ) + abs(len(left) - len(right))
    checked = max(len(left), len(right))
    return {
        "checked_records": checked,
        "validation_mismatches": mismatch_count,
        "agreement_rate": (checked - mismatch_count) / checked if checked else 0.0,
        "comparison_fields": fields,
        "first_mismatches": mismatches,
    }


def run_command(command: list[str]) -> tuple[float, str, str]:
    begin = time.perf_counter()
    process = subprocess.run(
        command,
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    seconds = time.perf_counter() - begin
    if process.returncode != 0:
        raise RuntimeError(
            f"command failed ({process.returncode}): {' '.join(command)}\n{process.stderr}"
        )
    return seconds, process.stdout.strip(), process.stderr.strip()


def dotmatch_version(dotmatch: Path) -> str:
    resolved = command_path(str(dotmatch))
    _seconds, stdout, stderr = run_command([str(resolved), "--version"])
    return stdout or stderr


def run_dotmatch(
    dotmatch: Path,
    targets: Path,
    reads: Path,
    start: int,
    length: int,
    k: int,
    work: Path,
) -> tuple[dict, dict]:
    prefix = work / f"dotmatch.k{k}"
    command = [
        str(dotmatch),
        "count",
        "--targets",
        str(targets),
        "--reads",
        str(reads),
        "--sample-label",
        "GSE146194_SRR11214031_heldout_prefix",
        "--target-start",
        str(start),
        "--target-length",
        str(length),
        "--k",
        str(k),
        "--metric",
        "hamming",
        "--ambiguity-policy",
        "radius",
        "--format",
        "dotmatch",
        "--out",
        str(prefix) + ".counts.tsv",
        "--summary",
        str(prefix) + ".summary.json",
        "--assignments",
        str(prefix) + ".assignments.tsv",
        "--ambiguous",
        "report",
        "--sample-qc",
        str(prefix) + ".sample_qc.tsv",
    ]
    seconds, stdout, stderr = run_command(command)
    assignments = Path(str(prefix) + ".assignments.tsv")
    return summarize_assignments(read_assignment_rows(assignments)), {
        "seconds": seconds,
        "command": public_command(command),
        "stdout": stdout,
        "stderr": stderr,
        "assignments": public_path(assignments),
        "summary": public_path(Path(str(prefix) + ".summary.json")),
        "counts": public_path(Path(str(prefix) + ".counts.tsv")),
        "sample_qc": public_path(Path(str(prefix) + ".sample_qc.tsv")),
    }


def library_audit(targets: Sequence[Target]) -> dict:
    pairs: list[dict[str, object]] = []
    minimum = None
    for left_index, left in enumerate(targets):
        for right in targets[left_index + 1 :]:
            distance = hamming(left.sequence, right.sequence)
            minimum = distance if minimum is None else min(minimum, distance)
            if distance <= 2:
                pairs.append(
                    {
                        "left": left.target_id,
                        "right": right.target_id,
                        "distance": distance,
                    }
                )
    return {
        "minimum_pairwise_hamming_distance": minimum,
        "pairs_within_hamming_1": sum(pair["distance"] <= 1 for pair in pairs),
        "pairs_within_hamming_2": len(pairs),
        "pairs_within_hamming_2_detail": pairs,
    }


def peak_rss_mib(who: int) -> float:
    value = float(resource.getrusage(who).ru_maxrss)
    return value / (1024 * 1024) if sys.platform == "darwin" else value / 1024


def run_analysis(
    dotmatch: Path,
    targets: Sequence[Target],
    targets_path: Path,
    reads: Path,
    start: int,
    length: int,
    work: Path,
) -> list[dict]:
    dotmatch = command_path(str(dotmatch))
    runs: list[dict] = []
    for k, comparator in [(0, "exact_slice_hash"), (1, "exhaustive_hamming_radius")]:
        dotmatch_summary, command = run_dotmatch(
            dotmatch, targets_path, reads, start, length, k, work
        )
        oracle_path = work / f"oracle.k{k}.assignments.tsv"
        oracle_summary, oracle_seconds = write_oracle(
            reads, targets, start, length, k, oracle_path
        )
        agreement = compare_assignments(command_path(command["assignments"]), oracle_path)
        runs.append(
            {
                "k": k,
                "metric": "hamming",
                "ambiguity_policy": "radius",
                "dotmatch": {**dotmatch_summary, **command},
                "comparator": {
                    "name": comparator,
                    "seconds": oracle_seconds,
                    "assignments": public_path(oracle_path),
                    **oracle_summary,
                },
                "agreement": agreement,
            }
        )
    return runs


def command_path(path_text: str) -> Path:
    path = Path(path_text)
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def write_csv(path: Path, results: dict) -> None:
    fields = [
        "case_study_id",
        "source_accession",
        "tool",
        "comparator",
        "k",
        "metric",
        "ambiguity_policy",
        "target_start",
        "target_length",
        "orientation",
        "reads",
        "targets",
        "assigned_unique",
        "assigned_exact",
        "assigned_corrected",
        "ambiguous",
        "unmatched",
        "invalid",
        "distinct_assigned_targets",
        "validation_mismatches",
        "agreement_rate",
        "seconds",
        "peak_rss_mib",
        "command",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for run in results["runs"]:
            stats = run["dotmatch"]
            writer.writerow(
                {
                    "case_study_id": results["case_study_id"],
                    "source_accession": "SRR11214031",
                    "tool": "dotmatch",
                    "comparator": run["comparator"]["name"],
                    "k": run["k"],
                    "metric": run["metric"],
                    "ambiguity_policy": run["ambiguity_policy"],
                    "target_start": results["window_discovery"]["target_start"],
                    "target_length": results["window_discovery"]["target_length"],
                    "orientation": results["window_discovery"]["orientation"],
                    "reads": stats["total_reads"],
                    "targets": results["guide_library"]["target_count"],
                    "assigned_unique": stats["assigned_unique"],
                    "assigned_exact": stats["assigned_exact"],
                    "assigned_corrected": stats["assigned_corrected"],
                    "ambiguous": stats["ambiguous"],
                    "unmatched": stats["unmatched"],
                    "invalid": stats["invalid"],
                    "distinct_assigned_targets": stats["distinct_assigned_targets"],
                    "validation_mismatches": run["agreement"]["validation_mismatches"],
                    "agreement_rate": f"{run['agreement']['agreement_rate']:.8f}",
                    "seconds": f"{stats['seconds']:.6f}",
                    "peak_rss_mib": f"{results['resources']['peak_rss_mib']:.3f}",
                    "command": stats["command"],
                }
            )


def percent(value: float) -> str:
    return f"{value * 100:.2f}%"


def report_markdown(results: dict, protocol: dict) -> str:
    window = results["window_discovery"]
    lines = [
        "# GSE146194 Direct-Guide-Capture Perturb-seq Case Study",
        "",
        "This bounded public-data case study tests DotMatch at its evidence boundary: per-read fixed-window assignment to a known multi-guide barcode library. It uses the UPR GBC sample from Replogle et al. and holds window-discovery reads out of the reported evaluation.",
        "",
        "## Result",
        "",
        f"The publisher supplement yielded `{results['guide_library']['target_count']}` guide barcodes. A frozen discovery rule selected `{window['orientation']}` orientation at zero-based start `{window['target_start']}` for an `{window['target_length']}`-base window, then `{results['input_slice']['evaluation_records']:,}` held-out reads were evaluated.",
        "",
        "| Rule | Unique | Exact | Corrected | Ambiguous | Unmatched | Invalid | Distinct guides | Oracle mismatches | Agreement |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for run in results["runs"]:
        stats = run["dotmatch"]
        lines.append(
            f"| Hamming k={run['k']} radius | {stats['assigned_unique']} | {stats['assigned_exact']} | {stats['assigned_corrected']} | {stats['ambiguous']} | {stats['unmatched']} | {stats['invalid']} | {stats['distinct_assigned_targets']} | {run['agreement']['validation_mismatches']} | {percent(run['agreement']['agreement_rate'])} |"
        )
    lines.extend(
        [
            "",
            "The unmatched and ambiguous columns are part of the result, not discarded failures. The small deterministic fixture separately requires exact, corrected, ambiguous, unmatched, and invalid outcomes. If this public slice contains zero ambiguous reads, that is reported as zero rather than manufactured.",
            "",
            "## Scientific question and criterion",
            "",
            f"**Question.** {protocol['scientific_question']}",
            "",
            "**Ground truth.** For k=0, a transparent exact-slice hash classifies target multiplicity. For k=1, an independent exhaustive Hamming-radius implementation checks every guide. The frozen pass criterion is zero held-out per-read differences in status, target identifier, and distance.",
            "",
            "The authors' published guide caller is not used as a per-read comparator: it combines guide-aligned reads with Cell Ranger-corrected cell barcodes and UMIs, then makes cell-level threshold calls. Those units and rules are not matched to this bounded per-read assignment question.",
            "",
            "## Dataset and provenance",
            "",
            "- GEO series/sample: [GSE146194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146194) / [GSM4367979](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4367979)",
            "- SRA study/experiment/run: [SRP251252](https://www.ncbi.nlm.nih.gov/sra/?term=SRP251252) / [SRX7826824](https://www.ncbi.nlm.nih.gov/sra/SRX7826824) / [SRR11214031](https://www.ncbi.nlm.nih.gov/sra/SRR11214031)",
            "- BioProject/BioSample: [PRJNA609688](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA609688) / [SAMN14258014](https://www.ncbi.nlm.nih.gov/biosample/SAMN14258014)",
            "- Primary paper: Replogle JM et al., *Nature Biotechnology* 38, 954-961 (2020), [doi:10.1038/s41587-020-0470-y](https://doi.org/10.1038/s41587-020-0470-y), [PMID 32231336](https://pubmed.ncbi.nlm.nih.gov/32231336/)",
            "- Access and reuse: [NCBI's molecular-data policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) places no NCBI restrictions on use or distribution but does not transfer any rights asserted by submitters. The [publisher supplement](https://support.springernature.com/en/support/solutions/articles/6000210902-supplementary-information) is retrieved from its article link; this workflow asserts no redistribution license and does not commit the workbook.",
            "",
            "Guide workbook SHA-256:",
            "",
            "```text",
            str(results["guide_library"]["workbook_sha256"]),
            "```",
            "",
            "Target table SHA-256:",
            "",
            "```text",
            str(results["guide_library"]["targets_sha256"]),
            "```",
            "",
            "Evaluation FASTQ SHA-256:",
            "",
            "```text",
            str(results["input_slice"]["evaluation_gzip_sha256"]),
            "```",
            "",
            "Selected-prefix uncompressed SHA-256:",
            "",
            "```text",
            str(results["input_slice"]["selected_prefix_uncompressed_sha256"]),
            "```",
            "",
            "The raw read archive is not committed. The workflow records ENA's full-file MD5 and byte count, streams only the frozen prefix, and verifies the derived prefix and evaluation files with SHA-256. The full archive MD5 is not claimed as locally reverified in bounded mode.",
            "",
            "## Reproduce",
            "",
            "```bash",
            "make dotmatch",
            "python3 scripts/run_perturb_seq_gse146194.py public --dotmatch ./dotmatch",
            "python3 scripts/check_perturb_seq_gse146194.py --require-public",
            "```",
            "",
            "For the no-network deterministic harness:",
            "",
            "```bash",
            "make perturb-seq-case-study-fixture-gate",
            "```",
            "",
            "Expected public outputs are `examples/perturb_seq_gse146194/results.json`, `examples/perturb_seq_gse146194/provenance.json`, `benchmarks/raw/perturb_seq_gse146194.csv`, and this report. Large reads and per-read work products stay under the ignored `examples/perturb_seq_gse146194/work/` directory.",
            "",
            "## Methods",
            "",
            f"The workflow extracts `{results['guide_library']['target_count']}` 18-base GBCs from Supplementary Table 2, scans every valid start and both target orientations on the first `{window['records']}` reads, and applies deterministic score and tie-break rules from `protocol.json`. Those discovery reads are excluded. DotMatch uses Hamming distance with conservative radius ambiguity at k=0 and k=1. Independent matched-rule oracles produce per-read status, target, and distance for the held-out prefix. Commands, versions, hashes, and resources are recorded in machine-readable artifacts.",
            "",
            "## QC interpretation",
            "",
        ]
    )
    for run in results["runs"]:
        stats = run["dotmatch"]
        lines.append(
            f"- k={run['k']}: assignment `{percent(stats['assignment_rate'])}`, ambiguous `{percent(stats['ambiguous_rate'])}`, unmatched `{percent(stats['unmatched_rate'])}`, invalid `{percent(stats['invalid_rate'])}`; `{stats['distinct_assigned_targets']}` distinct guides received unique reads."
        )
    audit = results["guide_library"]["audit"]
    lines.extend(
        [
            f"- Guide-library audit: minimum pairwise Hamming distance `{audit['minimum_pairwise_hamming_distance']}`; `{audit['pairs_within_hamming_1']}` pairs within distance 1 and `{audit['pairs_within_hamming_2']}` pairs within distance 2.",
            f"- Resource record for this bounded run: `{results['resources']['wall_seconds']:.3f}` wall seconds, `{results['resources']['peak_rss_mib']:.1f}` MiB peak RSS, and `{results['input_slice']['compressed_network_bytes']}` compressed bytes read from the FASTQ stream.",
            "",
            "## What this proves",
            "",
            "This proves reproducible, matched-rule per-read assignment and explicit QC behavior for the checked held-out SRR11214031 prefix and the published 32-guide GBC list. It closes the earlier single-guide extraction gap with a real multi-guide direct-capture dataset.",
            "",
            "## What this does not prove",
            "",
            "It does not prove guide-per-cell accuracy, UMI deduplication, Cell Ranger parity, expression quantification, perturbation effects, biological validity, full-run prevalence, or a speed advantage. The bounded prefix and fixed-position assumption can fail to represent other runs or protocols.",
            "",
            "## Next step",
            "",
            "A core facility or Perturb-seq workflow maintainer can reuse the manifest and held-out protocol on a complete guide-enrichment run, then add matched cell-barcode/UMI aggregation and compare guide-per-cell calls with the authors' published `cell_identities.csv` under explicitly matched thresholds.",
            "",
        ]
    )
    return "\n".join(lines)


def report_html(markdown_results: dict, protocol: dict) -> str:
    rows = []
    for run in markdown_results["runs"]:
        stats = run["dotmatch"]
        rows.append(
            "<tr>"
            f"<td>Hamming k={run['k']} radius</td><td>{stats['assigned_unique']}</td>"
            f"<td>{stats['assigned_exact']}</td><td>{stats['assigned_corrected']}</td>"
            f"<td>{stats['ambiguous']}</td><td>{stats['unmatched']}</td>"
            f"<td>{stats['invalid']}</td><td>{stats['distinct_assigned_targets']}</td>"
            f"<td>{run['agreement']['validation_mismatches']}</td>"
            f"<td>{html.escape(percent(run['agreement']['agreement_rate']))}</td></tr>"
        )
    window = markdown_results["window_discovery"]
    return f"""<!doctype html>
<html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">
<title>GSE146194 Direct-Guide-Capture Perturb-seq Case Study</title>
<style>body{{font:16px/1.55 system-ui,-apple-system,sans-serif;color:#17212b;background:#f7f8fa;margin:0}}main{{box-sizing:border-box;max-width:1000px;margin:auto;padding:48px 24px 72px}}h1{{font-size:2.35rem;line-height:1.12}}h2{{margin-top:2.2rem}}.lead,.boundary{{background:white;border:1px solid #d8dee6;border-radius:12px;padding:20px}}.table-wrap{{overflow-x:auto;border-radius:8px;-webkit-overflow-scrolling:touch}}table{{border-collapse:collapse;width:100%;min-width:880px;background:white}}th,td{{padding:10px;border-bottom:1px solid #d8dee6;text-align:right}}th:first-child,td:first-child{{text-align:left}}code{{overflow-wrap:anywhere}}a{{color:#075fa8}}.boundary{{border-left:5px solid #8a5a00}}@media(max-width:600px){{main{{padding:32px 24px 56px}}h1{{font-size:2rem}}}}</style></head>
<body><main><h1>GSE146194 direct-guide-capture Perturb-seq</h1>
<p class=\"lead\">A bounded, accessioned multi-guide case study for deterministic per-read assignment. The first {window['records']:,} records select a frozen fixed window; {markdown_results['input_slice']['evaluation_records']:,} held-out records carry the reported result.</p>
<h2>Matched-rule agreement</h2><div class="table-wrap" role="region" aria-label="Matched-rule agreement table" tabindex="0"><table><thead><tr><th>Rule</th><th>Unique</th><th>Exact</th><th>Corrected</th><th>Ambiguous</th><th>Unmatched</th><th>Invalid</th><th>Guides</th><th>Mismatches</th><th>Agreement</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div>
<h2>Question and method</h2><p>{html.escape(protocol['scientific_question'])}</p><p>The published 32-guide GBC table is extracted from Supplementary Table 2. DotMatch Hamming-radius calls at k=0 and k=1 are checked record-for-record against independent exact-slice and exhaustive Hamming oracles.</p>
<h2>Dataset identity</h2><p><a href=\"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146194\">GSE146194</a> · <a href=\"https://www.ncbi.nlm.nih.gov/sra/SRR11214031\">SRR11214031</a> · <a href=\"https://doi.org/10.1038/s41587-020-0470-y\">primary paper</a></p><p><a href=\"https://www.ncbi.nlm.nih.gov/home/about/policies/\">NCBI's molecular-data policy</a> places no NCBI restrictions on use or distribution but does not transfer submitter rights. The publisher workbook is retrieved from the article and is not redistributed or relicensed here.</p>
<p>Window: {window['orientation']} orientation, zero-based start {window['target_start']}, length {window['target_length']}. Evaluation SHA-256: <code>{markdown_results['input_slice']['evaluation_gzip_sha256']}</code>.</p>
<h2>Evidence boundary</h2><div class=\"boundary\"><strong>Proves:</strong> matched-rule per-read assignment and explicit QC on the checked held-out prefix. <strong>Does not prove:</strong> cell/UMI calls, expression quantification, Cell Ranger parity, perturbation effects, biological validity, full-run prevalence, or speed superiority.</div>
<h2>Next step</h2><p>Run the same manifest on a complete guide-enrichment library, then add explicitly matched cell-barcode/UMI aggregation and compare guide-per-cell calls with the authors' published cell identities.</p>
</main></body></html>"""


def write_reports(results: dict, protocol: dict, markdown_path: Path, html_path: Path) -> None:
    markdown_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    markdown_path.write_text(report_markdown(results, protocol), encoding="utf-8")
    html_path.write_text(report_html(results, protocol), encoding="utf-8")


def public_run(args: argparse.Namespace) -> int:
    started = time.perf_counter()
    protocol_path = Path(args.protocol)
    protocol = load_json(protocol_path)
    if protocol.get("protocol_status") != "frozen_before_read_inspection":
        raise RuntimeError("protocol must be frozen before public analysis")
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    workbook_path = work / "guide_library.xlsx"
    workbook = protocol["inputs"]["guide_library"]
    workbook_record = download_verified(
        str(workbook["url"]),
        workbook_path,
        str(workbook["sha256"]),
        int(protocol["resource_bounds"]["maximum_local_work_bytes"]),
        int(protocol["resource_bounds"]["network_timeout_seconds"]),
    )
    extracted_targets_path = work / "targets.published.tsv"
    published_targets = extract_targets(workbook_path, protocol, extracted_targets_path)
    evaluation_fastq = work / "evaluation.fastq.gz"
    discovery_records, input_slice = stream_prefix(protocol, evaluation_fastq)
    oriented_targets, window = discover_window(discovery_records, published_targets, protocol)
    targets_path = work / "targets.tsv"
    write_targets(targets_path, oriented_targets)
    runs = run_analysis(
        Path(args.dotmatch),
        oriented_targets,
        targets_path,
        evaluation_fastq,
        int(window["target_start"]),
        int(window["target_length"]),
        work,
    )
    if any(run["agreement"]["validation_mismatches"] for run in runs):
        raise RuntimeError("DotMatch disagreed with a matched independent oracle")
    results = {
        "schema_version": 1,
        "case_study_id": protocol["case_study_id"],
        "protocol": public_path(protocol_path),
        "dataset_accessions": protocol["dataset"]["accessions"],
        "access_and_reuse": access_and_reuse_record(),
        "guide_library": {
            "target_count": len(oriented_targets),
            "workbook_sha256": workbook_record["sha256"],
            "published_targets_sha256": sha256_path(extracted_targets_path),
            "targets_sha256": sha256_path(targets_path),
            "orientation": window["orientation"],
            "audit": library_audit(oriented_targets),
        },
        "input_slice": input_slice,
        "window_discovery": window,
        "runs": runs,
        "resources": {
            "wall_seconds": time.perf_counter() - started,
            "peak_rss_mib": max(
                peak_rss_mib(resource.RUSAGE_SELF),
                peak_rss_mib(resource.RUSAGE_CHILDREN),
            ),
            "host": platform.platform(),
            "python": platform.python_version(),
        },
    }
    results_path = Path(args.results)
    provenance = {
        "schema_version": 1,
        "case_study_id": protocol["case_study_id"],
        "protocol_sha256": sha256_path(protocol_path),
        "dotmatch_version": dotmatch_version(Path(args.dotmatch)),
        "access_and_reuse": access_and_reuse_record(),
        "source_artifacts": {
            "guide_workbook": workbook_record,
            "fastq": protocol["inputs"]["fastq"],
        },
        "derived_artifacts": {
            "published_targets": {
                "path": public_path(extracted_targets_path),
                "sha256": sha256_path(extracted_targets_path),
            },
            "oriented_targets": {
                "path": public_path(targets_path),
                "sha256": sha256_path(targets_path),
            },
            "evaluation_fastq": {
                "path": public_path(evaluation_fastq),
                "sha256": sha256_path(evaluation_fastq),
            },
        },
        "commands": [run["dotmatch"]["command"] for run in runs],
        "software": results["resources"],
        "evidence_boundary": {
            "proves": "matched-rule per-read fixed-window guide-barcode assignment on the held-out bounded SRR11214031 prefix",
            "does_not_prove": protocol["exclusions"],
        },
    }
    stable_json(results_path, results)
    stable_json(Path(args.provenance), provenance)
    write_csv(Path(args.csv), results)
    write_reports(results, protocol, Path(args.report), Path(args.html))
    print(f"PUBLIC CASE STUDY: PASS ({input_slice['evaluation_records']} held-out reads)")
    return 0


def load_targets(path: Path) -> list[Target]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return [
            Target(row["target_id"], row["target_seq"].upper(), row.get("gene", ""))
            for row in rows
        ]


def fixture_run(args: argparse.Namespace) -> int:
    fixture = Path(args.fixture)
    targets_path = fixture / "targets.tsv"
    reads = fixture / "reads.fastq"
    expected = load_json(fixture / "expected.json")
    targets = load_targets(targets_path)
    work = Path(args.work) / "fixture"
    work.mkdir(parents=True, exist_ok=True)
    runs = run_analysis(
        Path(args.dotmatch),
        targets,
        targets_path,
        reads,
        int(expected["target_start"]),
        int(expected["target_length"]),
        work,
    )
    failures: list[str] = []
    for run in runs:
        wanted = expected["runs"][f"k{run['k']}"]
        observed = run["dotmatch"]
        for field, wanted_value in wanted.items():
            observed_value = (
                run["agreement"][field]
                if field == "validation_mismatches"
                else observed[field]
            )
            if observed_value != wanted_value:
                failures.append(
                    f"k={run['k']} {field}: expected {wanted_value}, observed {observed_value}"
                )
    fixture_results = {
        "schema_version": 1,
        "note": expected["note"],
        "runs": runs,
        "failures": failures,
    }
    stable_json(work / "results.json", fixture_results)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        print("PERTURB-SEQ CASE-STUDY FIXTURE: FAIL")
        return 1
    print("PERTURB-SEQ CASE-STUDY FIXTURE: PASS")
    return 0


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    subparsers = root.add_subparsers(dest="command", required=True)
    public = subparsers.add_parser("public", help="run the bounded public-data workflow")
    public.add_argument("--protocol", default=str(DEFAULT_PROTOCOL))
    public.add_argument("--work", default=str(DEFAULT_WORK))
    public.add_argument("--dotmatch", default=str(ROOT / "dotmatch"))
    public.add_argument("--results", default=str(DEFAULT_RESULTS))
    public.add_argument("--provenance", default=str(DEFAULT_PROVENANCE))
    public.add_argument("--csv", default=str(DEFAULT_CSV))
    public.add_argument("--report", default=str(DEFAULT_REPORT))
    public.add_argument("--html", default=str(DEFAULT_HTML))
    public.set_defaults(func=public_run)
    fixture = subparsers.add_parser("fixture", help="run the no-network deterministic fixture")
    fixture.add_argument("--fixture", default=str(DEFAULT_FIXTURE))
    fixture.add_argument("--work", default=str(DEFAULT_WORK))
    fixture.add_argument("--dotmatch", default=str(ROOT / "dotmatch"))
    fixture.set_defaults(func=fixture_run)
    return root


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        return int(args.func(args))
    except (OSError, RuntimeError, ValueError, zipfile.BadZipFile) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
