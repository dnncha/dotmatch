"""One-pass comparison of exact, radius-one and nearest-distance assignments.

The comparison reuses the native Hamming index. It is a sensitivity analysis,
not a policy-selection rule or an estimate of biological assignment accuracy.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import io
import json
import os
import re
import shutil
import sys
import tempfile
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence

from .core import (
    MATCH_AMBIGUOUS,
    MATCH_INVALID,
    MATCH_NONE,
    MATCH_UNIQUE,
    Matcher,
    status_name,
)
from .target_io import read_target_table

MODES = ("exact", "radius_k1", "best_k1")
STATES = ("unique", "ambiguous", "none", "invalid")
PAIRS = (("exact", "radius_k1"), ("exact", "best_k1"), ("radius_k1", "best_k1"))


@dataclass(frozen=True)
class PolicyCall:
    status: int
    target_index: int = -1


@dataclass(frozen=True)
class PolicyComparison:
    exact: PolicyCall
    radius_k1: PolicyCall
    best_k1: PolicyCall
    candidates_within_one: int

    @property
    def changed(self) -> bool:
        return not (self.exact == self.radius_k1 == self.best_k1)


def compare_hamming_policies(
    matcher: Matcher, windows: Sequence[str | None]
) -> list[PolicyComparison]:
    """Compare three policies with one native query per distinct valid window.

    ``None`` denotes a failed window extraction. A tied best match stays tied;
    a radius policy also abstains when a farther target is within the radius.
    Exact status can be derived from the same result because Hamming distance
    is nonnegative. No exact-match candidate is skipped by this query.
    """
    distinct = list(dict.fromkeys(window for window in windows if window is not None))
    native = matcher.assign_hamming(distinct, k=1, policy="best") if distinct else []
    calls: dict[str | None, PolicyComparison] = {}
    invalid = PolicyCall(MATCH_INVALID)
    calls[None] = PolicyComparison(invalid, invalid, invalid, 0)
    for window, result in zip(distinct, native):
        best = PolicyCall(
            result.status, result.target_index if result.status == MATCH_UNIQUE else -1
        )
        radius = PolicyCall(MATCH_AMBIGUOUS) if result.match_count > 1 else best
        exact = (
            best
            if result.best_distance == 0 or result.status == MATCH_INVALID
            else PolicyCall(MATCH_NONE)
        )
        calls[window] = PolicyComparison(exact, radius, best, result.match_count)
    return [calls[window] for window in windows]


class _HashReader(io.RawIOBase):
    def __init__(self, raw, digest):
        super().__init__()
        self.raw, self.digest = raw, digest

    def readable(self):
        return True

    def readinto(self, buffer):
        count = self.raw.readinto(buffer)
        if count:
            self.digest.update(memoryview(buffer)[:count])
        return count


def _fastq(handle, source: Path) -> Iterator[tuple[str, str]]:
    ordinal = 0
    while True:
        header = handle.readline()
        if not header:
            return
        ordinal += 1
        seq, plus, qual = handle.readline(), handle.readline(), handle.readline()
        if not seq or not plus or not qual:
            raise ValueError(f"truncated FASTQ record {ordinal} in {source}")
        header, seq, plus, qual = (
            text.rstrip("\r\n") for text in (header, seq, plus, qual)
        )
        if (
            not header.startswith("@")
            or not header[1:].split()
            or not plus.startswith("+")
        ):
            raise ValueError(f"invalid FASTQ record {ordinal} in {source}")
        if (
            len(seq) != len(qual)
            or not seq.isascii()
            or any(ch.isspace() for ch in seq)
        ):
            raise ValueError(
                f"invalid sequence/quality lengths or sequence symbols at record {ordinal} in {source}"
            )
        if any(ord(ch) < 33 or ord(ch) > 126 for ch in qual):
            raise ValueError(
                f"invalid Phred+33 quality at record {ordinal} in {source}"
            )
        yield header[1:].split()[0], seq.upper()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _writer(path: Path, fields: Sequence[str]):
    handle = path.open("w", encoding="utf-8", newline="")
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(fields)
    return handle, writer


def _report(summary: dict, changes: list[dict]) -> str:
    escape = lambda value: html.escape(str(value))
    rows = "".join(
        "<tr><th>"
        + mode
        + "</th>"
        + "".join(f"<td>{summary['outcomes'][mode][state]:,}</td>" for state in STATES)
        + "</tr>"
        for mode in MODES
    )
    guide_rows = "".join(
        f"<tr><th>{escape(row['target_id'])}</th><td>{row['exact']}</td><td>{row['radius_k1']}</td><td>{row['best_k1']}</td></tr>"
        for row in changes[:50]
    )
    if not guide_rows:
        guide_rows = '<tr><td colspan="4">No per-guide count differences.</td></tr>'
    return f"""<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>DotMatch | Assignment sensitivity</title><style>
:root{{color-scheme:light}}*{{box-sizing:border-box}}body{{margin:0;color:#14271f;background:#fff;font:16px/1.65 system-ui,sans-serif}}main{{max-width:1040px;margin:auto;padding:40px 24px}}h1{{font-size:clamp(2rem,5vw,3.4rem);line-height:1.1;letter-spacing:-.045em}}h2{{margin-top:40px;font-size:1.4rem}}.eyebrow{{color:#356756;font-size:.8rem;letter-spacing:.12em;text-transform:uppercase}}.lede{{font-size:1.2rem;max-width:780px}}.scroll{{overflow-x:auto}}table{{border-collapse:collapse;min-width:580px;width:100%;text-align:left}}th,td{{border-bottom:1px solid #dce7e2;padding:13px}}thead{{background:#f3f7f5}}.note{{border-left:3px solid #39715c;padding:16px 20px;background:#f3f7f5}}code{{overflow-wrap:anywhere}}footer{{margin-top:44px;font-size:.86rem;color:#4d6258}}a{{color:#14583d}}:focus-visible{{outline:3px solid #357863;outline-offset:3px}}
</style><main><p class="eyebrow">DotMatch · {escape(summary['sample_label'])} · local analysis</p>
<h1>What changes when you allow one mismatch?</h1><p class="lede">{summary['read_count']:,} reads. {summary['target_count']:,} targets. Three policies applied to the same fixed windows. No policy is selected automatically.</p>
<p class="note"><strong>{summary['changed_reads']:,} reads change assignment outcome or target between policies.</strong> A larger assigned fraction does not establish greater biological accuracy.</p>
<h2>Every read is accounted for</h2><div class="scroll" tabindex="0" role="region" aria-label="Policy outcomes"><table><thead><tr><th>Policy</th><th>Unique</th><th>Ambiguous</th><th>Unmatched</th><th>Invalid window</th></tr></thead><tbody>{rows}</tbody></table></div>
<p><strong>Exact:</strong> zero substitutions. <strong>Radius k=1:</strong> one and only one target within one substitution. <strong>Best k=1:</strong> one nearest target, allowing one substitution; ties stay ambiguous.</p>
<h2>Which guide counts changed?</h2><p>Up to 50 guides with the largest absolute count changes are shown. <a href="guide_deltas.tsv">guide_deltas.tsv</a> contains every supplied target, including unchanged and zero-count targets.</p><div class="scroll" tabindex="0" role="region" aria-label="Changed guide counts"><table><thead><tr><th>Target ID</th><th>Exact</th><th>Radius k=1</th><th>Best k=1</th></tr></thead><tbody>{guide_rows}</tbody></table></div>
<h2>Reproduce and inspect</h2><p><a href="summary.json">Machine-readable summary and checksums</a> · <a href="transitions.tsv">Read-state transitions</a> · <a href="sample_qc.tsv">Sample QC</a></p>
<p>Raw count tables: <a href="exact.counts.tsv">exact</a> · <a href="radius_k1.counts.tsv">radius k=1</a> · <a href="best_k1.counts.tsv">best k=1</a>.</p>
<footer>Hamming substitutions only; same supplied library, orientation and read window. N and other ASCII symbols compare literally, not as wildcards. No offsets, quality thresholds, indels, cell calls or downstream screen statistics are fitted here. Unique is an assignment outcome, not proof of biological origin. This report contains target identifiers but not raw read sequences.</footer></main></html>"""


def run_sensitivity(
    *,
    targets: str | Path,
    reads: str | Path,
    target_start: int,
    target_length: int,
    out_dir: str | Path,
    sample_label: str = "sample",
    batch_size: int = 4096,
    write_read_changes: bool = False,
) -> dict:
    """Write a no-overwrite review bundle; reject incomplete/invalid input runs."""
    if target_start < 0 or target_length <= 0:
        raise ValueError(
            "target_start must be non-negative and target_length must be positive"
        )
    if not 1 <= batch_size <= 65536:
        raise ValueError("batch_size must be between 1 and 65536")
    if not re.fullmatch(
        r"[A-Za-z0-9_.-]{1,120}", sample_label
    ) or sample_label.lower() in {"sgrna", "gene"}:
        raise ValueError(
            "sample_label must be 1–120 letters, digits, dots, underscores or hyphens, and not a reserved column name"
        )
    library_path, fastq_path, output = Path(targets), Path(reads), Path(out_dir)
    library_digest = _sha256(library_path)
    library = read_target_table(library_path)
    if _sha256(library_path) != library_digest:
        raise ValueError("target library changed during input validation")
    if any(len(row.sequence) != target_length for row in library):
        raise ValueError(
            "all target sequences must match target_length for this fixed-window Hamming comparison"
        )
    if not fastq_path.is_file():
        raise ValueError(f"FASTQ is not a regular file: {fastq_path}")
    # Reserve an empty directory exclusively; never replace an earlier run.
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir(exist_ok=False)
    try:
        stage = Path(tempfile.mkdtemp(prefix=".pending-", dir=output))
    except BaseException:
        output.rmdir()
        raise
    published: list[Path] = []
    started = time.perf_counter()
    counts = {mode: [0] * len(library) for mode in MODES}
    outcomes = {mode: Counter({state: 0 for state in STATES}) for mode in MODES}
    transitions = {pair: Counter() for pair in PAIRS}
    read_count = changed_reads = distinct_windows_queried = 0
    read_digest = hashlib.sha256()
    changes_handle = None
    try:
        change_writer = None
        if write_read_changes:
            changes_handle, change_writer = _writer(
                stage / "read_changes.tsv",
                ["record_index", "read_id"]
                + [
                    f"{mode}_{kind}"
                    for mode in MODES
                    for kind in ("status", "target_id")
                ],
            )
        with (
            Matcher([row.sequence for row in library]) as matcher,
            fastq_path.open("rb") as raw,
        ):
            before = os.fstat(raw.fileno())
            buffered = io.BufferedReader(_HashReader(raw, read_digest))
            decoded = (
                gzip.GzipFile(fileobj=buffered, mode="rb")
                if fastq_path.name.lower().endswith(".gz")
                else buffered
            )
            with io.TextIOWrapper(decoded, encoding="utf-8", newline="") as handle:
                iterator = iter(_fastq(handle, fastq_path))
                while True:
                    batch = []
                    for _ in range(batch_size):
                        record = next(iterator, None)
                        if record is None:
                            break
                        batch.append(record)
                    if not batch:
                        break
                    end = target_start + target_length
                    windows = [
                        sequence[target_start:end] if len(sequence) >= end else None
                        for _, sequence in batch
                    ]
                    distinct_windows_queried += len(
                        set(window for window in windows if window is not None)
                    )
                    results = compare_hamming_policies(matcher, windows)
                    for (read_id, _sequence), result in zip(batch, results):
                        read_count += 1
                        for mode in MODES:
                            call = getattr(result, mode)
                            outcomes[mode][status_name(call.status)] += 1
                            if call.status == MATCH_UNIQUE:
                                counts[mode][call.target_index] += 1
                        for left, right in PAIRS:
                            transitions[(left, right)][
                                (
                                    status_name(getattr(result, left).status),
                                    status_name(getattr(result, right).status),
                                )
                            ] += 1
                        if result.changed:
                            changed_reads += 1
                            if change_writer is not None:
                                values = [read_count, read_id]
                                for mode in MODES:
                                    call = getattr(result, mode)
                                    values.extend(
                                        [
                                            status_name(call.status),
                                            (
                                                library[call.target_index].target_id
                                                if call.status == MATCH_UNIQUE
                                                else ""
                                            ),
                                        ]
                                    )
                                change_writer.writerow(values)
            after = os.fstat(raw.fileno())
            current = fastq_path.stat()
            signature = lambda stat: (
                stat.st_dev,
                stat.st_ino,
                stat.st_size,
                stat.st_mtime_ns,
            )
            if signature(before) != signature(after) or signature(after) != signature(
                current
            ):
                raise ValueError(
                    "FASTQ changed during analysis; no completed report was produced"
                )
        if changes_handle is not None:
            changes_handle.close()
        if read_count == 0:
            raise ValueError("FASTQ contains no records")
        for mode in MODES:
            if (
                sum(outcomes[mode].values()) != read_count
                or sum(counts[mode]) != outcomes[mode]["unique"]
            ):
                raise RuntimeError("read conservation check failed")
            with (stage / f"{mode}.counts.tsv").open(
                "w", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["sgRNA", "Gene", sample_label])
                writer.writerows(
                    (row.target_id, row.gene, counts[mode][i])
                    for i, row in enumerate(library)
                )
        changed_guides = []
        handle, writer = _writer(
            stage / "guide_deltas.tsv",
            ["target_id", "gene", *MODES, "radius_minus_exact", "best_minus_exact"],
        )
        with handle:
            for i, row in enumerate(library):
                values = {mode: counts[mode][i] for mode in MODES}
                writer.writerow(
                    [
                        row.target_id,
                        row.gene,
                        *(values[mode] for mode in MODES),
                        values["radius_k1"] - values["exact"],
                        values["best_k1"] - values["exact"],
                    ]
                )
                if len(set(values.values())) > 1:
                    changed_guides.append({"target_id": row.target_id, **values})
        changed_guides.sort(
            key=lambda row: (
                -(max(row[mode] for mode in MODES) - min(row[mode] for mode in MODES)),
                row["target_id"],
            )
        )
        handle, writer = _writer(
            stage / "transitions.tsv",
            ["from_policy", "to_policy", "from_status", "to_status", "reads"],
        )
        with handle:
            for pair in PAIRS:
                for left in STATES:
                    for right in STATES:
                        writer.writerow(
                            [*pair, left, right, transitions[pair][(left, right)]]
                        )
        handle, writer = _writer(
            stage / "sample_qc.tsv", ["sample", "policy", "reads", *STATES]
        )
        with handle:
            for mode in MODES:
                writer.writerow(
                    [
                        sample_label,
                        mode,
                        read_count,
                        *(outcomes[mode][state] for state in STATES),
                    ]
                )
        from . import __version__
        from .core import _LIB

        summary = {
            "schema_version": "dotmatch.sensitivity.v1",
            "completion": "complete",
            "software_version": __version__,
            "implementation_sha256": _sha256(Path(__file__)),
            "native_library_sha256": _sha256(Path(_LIB._name)),
            "sample_label": sample_label,
            "read_count": read_count,
            "target_count": len(library),
            "parameters": {
                "target_start": target_start,
                "target_length": target_length,
                "metric": "hamming",
                "orientation": "as_supplied",
                "batch_size": batch_size,
                "write_read_changes": write_read_changes,
            },
            "inputs": {
                "targets": {"name": library_path.name, "sha256": library_digest},
                "reads": {"name": fastq_path.name, "sha256": read_digest.hexdigest()},
            },
            "outcomes": {mode: dict(outcomes[mode]) for mode in MODES},
            "changed_reads": changed_reads,
            "count_comparisons": [
                {
                    "left": left,
                    "right": right,
                    "counts_identical": counts[left] == counts[right],
                    "differing_guides": sum(
                        a != b for a, b in zip(counts[left], counts[right])
                    ),
                    "total_count_delta": sum(counts[right]) - sum(counts[left]),
                }
                for left, right in PAIRS
            ],
            "execution": {
                "fastq_passes": 1,
                "distinct_windows_queried_across_batches": distinct_windows_queried,
                "seconds": round(time.perf_counter() - started, 6),
            },
            "interpretation": "Assignment-policy sensitivity only; no automatic policy selection, biological truth labels or false-assignment estimate.",
            "limitations": [
                "Fixed windows, supplied orientation, substitutions only; no offset inference or indels.",
                "N/IUPAC symbols compare literally. Quality characters are validated but not used to score assignments.",
                "Count identity is distinct from read-assignment identity. No read sequence or quality is copied to this bundle.",
            ],
        }
        (stage / "report.html").write_text(
            _report(summary, changed_guides), encoding="utf-8"
        )
        summary["artifacts"] = {
            path.name: {"sha256": _sha256(path), "bytes": path.stat().st_size}
            for path in sorted(stage.iterdir())
        }
        (stage / "summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        # Hard-link publication is no-clobber, even if another process creates a file.
        # Publish the completion manifest last. All links share this filesystem.
        for path in sorted(
            stage.iterdir(), key=lambda path: (path.name == "summary.json", path.name)
        ):
            destination = output / path.name
            os.link(path, destination)
            published.append(destination)
        shutil.rmtree(stage)
        return summary
    except BaseException:
        if changes_handle is not None:
            changes_handle.close()
        for path in published:
            path.unlink(missing_ok=True)
        shutil.rmtree(stage, ignore_errors=True)
        try:
            output.rmdir()
        except OSError:
            pass
        raise


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="dotmatch sensitivity", description=__doc__)
    parser.add_argument(
        "--targets", required=True, help="Target TSV/CSV, optionally gzipped"
    )
    parser.add_argument("--reads", required=True, help="One FASTQ or FASTQ.gz sample")
    parser.add_argument(
        "--target-start",
        required=True,
        type=int,
        help="Zero-based fixed read-window start",
    )
    parser.add_argument("--target-length", required=True, type=int)
    parser.add_argument(
        "--out-dir",
        required=True,
        help="New output directory; existing directories are refused",
    )
    parser.add_argument("--sample-label", default="sample")
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument(
        "--write-read-changes",
        action="store_true",
        help="Write changed record IDs and policy calls, not raw sequences",
    )
    args = parser.parse_args(argv)
    try:
        summary = run_sensitivity(**vars(args))
    except (OSError, EOFError, ValueError, csv.Error, RuntimeError) as exc:
        print(f"dotmatch sensitivity: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "completion": "complete",
                "report": str(Path(args.out_dir) / "report.html"),
                "summary": str(Path(args.out_dir) / "summary.json"),
                "read_count": summary["read_count"],
                "changed_reads": summary["changed_reads"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
