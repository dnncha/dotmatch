#!/usr/bin/env python3
"""Benchmark DotMatch classic-BCL demultiplexing.

The default workload is a generated classic per-cycle BCL run folder. It is a
smoke benchmark for parser/output plumbing, not a state-of-the-art claim.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import hashlib
import shutil
import struct
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw" / "bcl_demux.csv"
WORK = ROOT / "benchmarks" / "work" / "bcl_demux"


def parse_time_rss(path: Path) -> int:
    if not path.exists():
        return 0
    for line in path.read_text(errors="replace").splitlines():
        stripped = line.strip()
        if sys.platform == "darwin" and "maximum resident set size" in stripped:
            try:
                return int(stripped.split()[0]) // 1024
            except (IndexError, ValueError):
                return 0
        if "Maximum resident set size" in stripped:
            try:
                return int(stripped.rsplit(":", 1)[1].strip())
            except (IndexError, ValueError):
                return 0
    return 0


def timed_command(cmd: list[str], time_out: Path) -> list[str]:
    if not Path("/usr/bin/time").exists():
        return cmd
    if sys.platform == "darwin":
        return ["/usr/bin/time", "-l", "-o", str(time_out), *cmd]
    return ["/usr/bin/time", "-v", "-o", str(time_out), *cmd]


def run(cmd: list[str], allow_missing: bool = False, log_prefix: str = "command") -> tuple[float, int, int, str, str]:
    log_dir = WORK / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    safe_prefix = "".join(c if c.isalnum() or c in "._-" else "_" for c in log_prefix)
    stdout_log = log_dir / f"{safe_prefix}.stdout.log"
    stderr_log = log_dir / f"{safe_prefix}.stderr.log"
    with tempfile.NamedTemporaryFile(prefix="dotmatch-bcl-time-", delete=False) as tmp:
        time_path = Path(tmp.name)
    start = time.perf_counter()
    try:
        with stdout_log.open("wb") as out, stderr_log.open("wb") as err:
            rc = subprocess.run(timed_command(cmd, time_path), cwd=ROOT, check=False, stdout=out, stderr=err).returncode
    except FileNotFoundError:
        if allow_missing:
            return 0.0, 127, 0, str(stdout_log), str(stderr_log)
        raise
    seconds = time.perf_counter() - start
    peak_rss = parse_time_rss(time_path)
    time_path.unlink(missing_ok=True)
    return seconds, rc, peak_rss, str(stdout_log), str(stderr_log)


def tool_version(exe: str, args: list[str]) -> str:
    found = shutil.which(exe)
    if found is None:
        return "not_installed"
    try:
        p = subprocess.run([found, *args], check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except Exception:
        return "unknown"
    lines = [line.strip() for line in p.stdout.splitlines() if line.strip()]
    for line in lines:
        lower = line.lower()
        if "version" in lower or " v" in lower or lower.startswith("v"):
            return line.replace(",", ";")
    for line in lines:
        return line.replace(",", ";")
    return "unknown"


def encode_base(base: str, q: int = 40) -> int:
    return (q << 2) | {"A": 0, "C": 1, "G": 2, "T": 3}[base]


def write_bcl(path: Path, bases: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as fh:
        fh.write(struct.pack("<I", len(bases)))
        fh.write(bytes(encode_base(b) for b in bases))


def make_classic_bcl_fixture(root: Path, clusters: int) -> tuple[Path, Path, int, int]:
    shutil.rmtree(root, ignore_errors=True)
    base = root / "Data" / "Intensities" / "BaseCalls" / "L001"
    for cycle in range(1, 9):
        (base / f"C{cycle}.1").mkdir(parents=True, exist_ok=True)
    (root / "RunInfo.xml").write_text("""<?xml version="1.0"?>
<RunInfo>
  <Run Id="bench" Number="1">
    <Flowcell>DOTMATCH</Flowcell>
    <Reads>
      <Read Number="1" NumCycles="4" IsIndexedRead="N"/>
      <Read Number="2" NumCycles="4" IsIndexedRead="Y"/>
    </Reads>
    <FlowcellLayout LaneCount="1" SurfaceCount="1" SwathCount="1" TileCount="1"/>
  </Run>
</RunInfo>
""")
    sheet = root / "SampleSheet.csv"
    barcodes = [
        ("s1", "Sample One", "ACGT"),
        ("s2", "Sample Two", "AGGT"),
        ("s3", "Sample Three", "ACGA"),
        ("s4", "Sample Four", "TTTT"),
    ]
    with sheet.open("w") as fh:
        fh.write("[Header]\nIEMFileVersion,4\n[Data]\nSample_ID,Sample_Name,index\n")
        for row in barcodes:
            fh.write(",".join(row) + "\n")
    read_payloads = ["AAAA", "CCCC", "GGGG", "TTTT"]
    read_cycles = [[] for _ in range(4)]
    index_cycles = [[] for _ in range(4)]
    pf = []
    for i in range(clusters):
        read = read_payloads[i % len(read_payloads)]
        index = barcodes[i % len(barcodes)][2]
        if i % 17 == 0:
            index = index[:-1] + ("A" if index[-1] != "A" else "C")
        if i % 43 == 0:
            index = "GGGG"
        for j in range(4):
            read_cycles[j].append(read[j])
            index_cycles[j].append(index[j])
        pf.append(0 if i % 101 == 0 else 1)
    for j, bases in enumerate(read_cycles, start=1):
        write_bcl(base / f"C{j}.1" / "s_1_1101.bcl.gz", bases)
    for j, bases in enumerate(index_cycles, start=5):
        write_bcl(base / f"C{j}.1" / "s_1_1101.bcl.gz", bases)
    with (base / "s_1_1101.filter").open("wb") as fh:
        fh.write(struct.pack("<II", 0, clusters))
        fh.write(bytes(pf))
    return root, sheet, len(barcodes), 8


def summary_stats(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    data = json.loads(path.read_text())
    return {
        "total_clusters": str(data.get("total_clusters", "")),
        "assigned_reads": str(data.get("assigned_reads", "")),
        "undetermined_reads": str(data.get("undetermined_reads", "")),
        "filtered_clusters": str(data.get("filtered_clusters", "")),
        "tiles": str(data.get("tiles", "")),
        "requested_threads": str(data.get("requested_threads", "")),
        "effective_threads": str(data.get("effective_threads", "")),
        "gzip_level": str(data.get("gzip_level", "")),
    }


def int_stat(stats: dict[str, str], key: str, default: int = 0) -> int:
    try:
        value = stats.get(key, "")
        if value == "":
            return default
        return int(value)
    except (TypeError, ValueError):
        return default


def run_info_cycles(run_folder: Path) -> int:
    path = run_folder / "RunInfo.xml"
    if not path.exists():
        return 0
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError:
        return 0
    total = 0
    for elem in root.iter():
        if elem.tag.endswith("Read"):
            try:
                total += int(elem.attrib.get("NumCycles", "0"))
            except ValueError:
                return 0
    return total


def sample_sheet_sample_count(path: Path) -> int:
    if not path.exists():
        return 0
    in_data = False
    header: list[str] | None = None
    sample_col = -1
    lane_col = -1
    seen: set[tuple[str, str]] = set()
    with path.open(newline="") as fh:
        for row in csv.reader(fh):
            if not row:
                continue
            first = row[0].strip()
            if first in ("[Data]", "[BCLConvert_Data]"):
                in_data = True
                header = None
                continue
            if first.startswith("["):
                in_data = False
                header = None
                continue
            if not in_data:
                continue
            if header is None:
                header = [x.strip() for x in row]
                candidates = ["Sample_ID", "SampleID", "Sample_Name", "Sample", "Sample_Project"]
                sample_col = next((header.index(c) for c in candidates if c in header), -1)
                lane_col = header.index("Lane") if "Lane" in header else -1
                continue
            if sample_col < 0 or sample_col >= len(row):
                continue
            sample = row[sample_col].strip()
            if not sample:
                continue
            lane = row[lane_col].strip() if lane_col >= 0 and lane_col < len(row) else ""
            seen.add((lane, sample))
    return len(seen)


def count_fastq_records(path: Path) -> int:
    opener = gzip.open if "".join(path.suffixes[-2:]) == ".fastq.gz" or path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return sum(1 for _ in fh) // 4


def output_stats(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    assigned = 0
    undetermined = 0
    fastqs = sorted(path.rglob("*.fastq.gz"))
    primary_fastqs = [p for p in fastqs if "_R1_" in p.name]
    if primary_fastqs:
        fastqs = primary_fastqs
    for fastq in fastqs:
        n = count_fastq_records(fastq)
        if "Undetermined" in fastq.name:
            undetermined += n
        else:
            assigned += n
    return {"assigned_reads": str(assigned), "undetermined_reads": str(undetermined)}


def output_fingerprint(path: Path) -> dict[str, str]:
    if not path.exists():
        return {"output_bytes": "", "output_sha256": "", "fastq_content_sha256": ""}
    total = 0
    artifact_digest = hashlib.sha256()
    fastq_digest = hashlib.sha256()
    saw_fastq = False
    for file in sorted(p for p in path.rglob("*") if p.is_file()):
        rel = file.relative_to(path).as_posix().encode()
        artifact_digest.update(rel)
        artifact_digest.update(b"\0")
        total += file.stat().st_size
        with file.open("rb") as fh:
            while True:
                chunk = fh.read(1024 * 1024)
                if not chunk:
                    break
                artifact_digest.update(chunk)
        if file.name.endswith(".fastq.gz") or file.name.endswith(".fq.gz"):
            saw_fastq = True
            fastq_digest.update(rel)
            fastq_digest.update(b"\0")
            with gzip.open(file, "rb") as fh:
                while True:
                    chunk = fh.read(1024 * 1024)
                    if not chunk:
                        break
                    fastq_digest.update(chunk)
    return {
        "output_bytes": str(total),
        "output_sha256": artifact_digest.hexdigest(),
        "fastq_content_sha256": fastq_digest.hexdigest() if saw_fastq else "",
    }


def validate_stats(dotmatch_out: Path, truth_out: Path) -> dict[str, str]:
    cmd = [str(ROOT / "dotmatch"), "bcl-validate", "--dotmatch-out", str(dotmatch_out), "--truth-out", str(truth_out)]
    p = subprocess.run(cmd, cwd=ROOT, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    try:
        data = json.loads(p.stdout)
    except json.JSONDecodeError:
        data = None
    if p.returncode == 0 and data is not None:
        mismatches = int(data.get("missing_fastq_files", 0)) + int(data.get("mismatched_fastq_files", 0))
        return {
            "validation_mismatches": str(mismatches),
            "validation_exit_code": str(p.returncode),
            "validation_mode": "strict_fastq",
        }
    dot = output_stats(dotmatch_out)
    truth = output_stats(truth_out)
    mismatches = 0
    for key in ("assigned_reads", "undetermined_reads"):
        if dot.get(key, "") != truth.get(key, ""):
            mismatches += 1
    return {
        "validation_mismatches": str(mismatches),
        "validation_exit_code": "0" if mismatches == 0 else "1",
        "validation_mode": "count_totals",
    }


def row(tool: str, version: str, workflow: str, fmt: str, clusters: int, cycles: int, samples: int,
        seconds: float, rc: int, peak_rss: int, cmd: list[str], stats: dict[str, str] | None = None) -> dict[str, str]:
    out = {
        "tool": tool,
        "version": version,
        "workflow": workflow,
        "format": fmt,
        "clusters": str(clusters),
        "cycles": str(cycles),
        "samples": str(samples),
        "seconds": f"{seconds:.6f}",
        "clusters_per_sec": f"{clusters / seconds:.1f}" if seconds > 0 and rc == 0 else "0.0",
        "peak_rss_kb": str(peak_rss),
        "total_clusters": str(clusters),
        "requested_threads": "",
        "effective_threads": "",
        "gzip_level": "",
        "output_bytes": "",
        "output_sha256": "",
        "fastq_content_sha256": "",
        "assigned_reads": "",
        "undetermined_reads": "",
        "filtered_clusters": "",
        "tiles": "",
        "validation_mismatches": "",
        "validation_exit_code": "",
        "validation_mode": "",
        "exit_code": str(rc),
        "command": " ".join(cmd),
        "stdout_log": "",
        "stderr_log": "",
    }
    if stats:
        out.update(stats)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clusters", type=int, default=int(os.environ.get("DOTMATCH_BCL_CLUSTERS", "20000")))
    parser.add_argument("--run-folder")
    parser.add_argument("--sample-sheet")
    parser.add_argument("--workflow-name", default="")
    parser.add_argument("--out", default=str(RAW))
    parser.add_argument("--detect-competitors", action="store_true")
    parser.add_argument("--run-installed-competitors", action="store_true")
    parser.add_argument("--threads", type=int, default=int(os.environ.get("DOTMATCH_BCL_THREADS", "1")))
    parser.add_argument("--gzip-level", type=int, default=int(os.environ.get("DOTMATCH_BCL_GZIP_LEVEL", "1")))
    parser.add_argument("--emit-index-fastqs", action="store_true")
    args = parser.parse_args()

    subprocess.run(["make", "dotmatch"], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    WORK.mkdir(parents=True, exist_ok=True)
    if args.run_folder and args.sample_sheet:
        run_folder = Path(args.run_folder).resolve()
        sample_sheet = Path(args.sample_sheet).resolve()
        workflow = args.workflow_name or "real_classic_bcl_user_supplied"
        samples = sample_sheet_sample_count(sample_sheet)
        cycles = run_info_cycles(run_folder)
    else:
        run_folder, sample_sheet, samples, cycles = make_classic_bcl_fixture(WORK / "classic_fixture", args.clusters)
        workflow = "synthetic_classic_bcl_fixture"

    out_dir = WORK / "dotmatch_bcl_out"
    shutil.rmtree(out_dir, ignore_errors=True)
    summary = WORK / "dotmatch_bcl_summary.json"
    cmd = [
        str(ROOT / "dotmatch"),
        "bcl-demux",
        "--run-folder", str(run_folder),
        "--sample-sheet", str(sample_sheet),
        "--out-dir", str(out_dir),
        "--barcode-mismatches", "1",
        "--threads", str(args.threads),
        "--gzip-level", str(args.gzip_level),
        "--summary", str(summary),
    ]
    if args.emit_index_fastqs:
        cmd.append("--emit-index-fastqs")
    seconds, rc, peak_rss, stdout_log, stderr_log = run(cmd, log_prefix="dotmatch_bcl_demux")
    dotmatch_stats = summary_stats(summary)
    dotmatch_stats.update(output_fingerprint(out_dir))
    dotmatch_stats.update({"stdout_log": stdout_log, "stderr_log": stderr_log})
    effective_clusters = int_stat(dotmatch_stats, "total_clusters", args.clusters)
    rows = [row("dotmatch_bcl_demux", "local", workflow, "classic_bcl", effective_clusters, cycles, samples,
                seconds, rc, peak_rss, cmd, dotmatch_stats)]

    if args.run_installed_competitors and shutil.which("bcl-convert") is not None:
        competitor_out = WORK / "bcl_convert_out"
        shutil.rmtree(competitor_out, ignore_errors=True)
        competitor_out.mkdir(parents=True, exist_ok=True)
        competitor_cmd = [
            "bcl-convert",
            "--force", "true",
            "--bcl-input-directory", str(run_folder),
            "--output-directory", str(competitor_out),
            "--sample-sheet", str(sample_sheet),
        ]
        seconds, comp_rc, comp_rss, stdout_log, stderr_log = run(competitor_cmd, allow_missing=True, log_prefix="bcl_convert")
        comp_stats = output_stats(competitor_out)
        comp_stats.update(output_fingerprint(competitor_out))
        comp_stats.update({"stdout_log": stdout_log, "stderr_log": stderr_log})
        if comp_rc == 0 and rc == 0:
            comp_stats.update(validate_stats(out_dir, competitor_out))
        rows.append(row("bcl-convert", tool_version("bcl-convert", ["--version"]), workflow, "classic_bcl_or_cbcl",
                        effective_clusters, cycles, samples, seconds, comp_rc, comp_rss, competitor_cmd, comp_stats))

    if args.run_installed_competitors and shutil.which("bcl2fastq") is not None:
        competitor_out = WORK / "bcl2fastq_out"
        shutil.rmtree(competitor_out, ignore_errors=True)
        competitor_out.mkdir(parents=True, exist_ok=True)
        competitor_cmd = [
            "bcl2fastq",
            "--runfolder-dir", str(run_folder),
            "--output-dir", str(competitor_out),
            "--sample-sheet", str(sample_sheet),
            "--ignore-missing-filter",
            "--ignore-missing-positions",
            "--ignore-missing-controls",
            "--fastq-compression-level", str(args.gzip_level),
        ]
        seconds, comp_rc, comp_rss, stdout_log, stderr_log = run(competitor_cmd, allow_missing=True, log_prefix="bcl2fastq")
        comp_stats = output_stats(competitor_out)
        comp_stats.update(output_fingerprint(competitor_out))
        comp_stats.update({"stdout_log": stdout_log, "stderr_log": stderr_log, "gzip_level": str(args.gzip_level)})
        if comp_rc == 0 and rc == 0:
            comp_stats.update(validate_stats(out_dir, competitor_out))
        rows.append(row("bcl2fastq", tool_version("bcl2fastq", ["--version"]), workflow, "classic_bcl",
                        effective_clusters, cycles, samples, seconds, comp_rc, comp_rss, competitor_cmd, comp_stats))

    if args.detect_competitors:
        for exe, version_args in [("bcl-convert", ["--version"]), ("bcl2fastq", ["--version"]), ("cuda-demux", ["--version"])]:
            if any(r["tool"] == exe for r in rows):
                continue
            rows.append(row(exe, tool_version(exe, version_args), workflow, "classic_bcl_or_cbcl",
                            effective_clusters, cycles, samples, 0.0, 127 if shutil.which(exe) is None else 125,
                            0, [exe, "not-run"]))

    out_csv = Path(args.out)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "tool", "version", "workflow", "format", "clusters", "cycles", "samples", "seconds",
        "clusters_per_sec", "peak_rss_kb", "total_clusters", "requested_threads", "effective_threads", "gzip_level",
        "output_bytes", "output_sha256", "fastq_content_sha256", "assigned_reads", "undetermined_reads",
        "filtered_clusters", "tiles", "validation_mismatches", "validation_exit_code", "validation_mode",
        "exit_code", "stdout_log", "stderr_log", "command",
    ]
    with out_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(out_csv)


if __name__ == "__main__":
    main()
