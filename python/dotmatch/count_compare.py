"""Compare existing raw-count tables without changing either counting workflow.

No FASTQs, normalization, network requests, tracking, or biological accuracy
claims. Run with ``python -m dotmatch.count_compare --help``.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import html
import json
import sys
import tempfile
from pathlib import Path
from typing import Iterator, Sequence

from .count_io import CountTable, read_count_table

SCHEMA_VERSION = 1
TOP_CHANGES = 20
IDENTITY_FIELDS = frozenset({"gene", "gene_id", "gene_symbol", "target_seq", "sequence", "seq"})
BOUNDARY = (
    "Agreement is not accuracy. Differences can come from reference libraries, "
    "extraction windows, filtering, assignment rules, or other preprocessing. "
    "Count tables alone cannot identify reassigned reads, independent molecules, "
    "or which result is biologically correct."
)


def _read_snapshot(path: str | Path) -> tuple[CountTable, dict]:
    """Hash the same private byte snapshot that is parsed, not a later reread."""
    source = Path(path)
    if not source.is_file():
        raise ValueError(f"count input must be a regular file: {source}")
    digest, size = hashlib.sha256(), 0
    with tempfile.TemporaryDirectory(prefix="dotmatch-count-compare-") as temporary:
        snapshot = Path(temporary) / ("counts.tsv.gz" if source.name.lower().endswith(".gz") else "counts.tsv")
        with source.open("rb") as incoming, snapshot.open("xb") as outgoing:
            for chunk in iter(lambda: incoming.read(1024 * 1024), b""):
                digest.update(chunk)
                size += len(chunk)
                outgoing.write(chunk)
        table = read_count_table(snapshot)
    # Do not put filesystem paths or filenames in a potentially shared report.
    return table, {"sha256": digest.hexdigest(), "bytes": size}


def _alignment(baseline: CountTable, candidate: CountTable, shared_only: bool) -> dict:
    left_ids, right_ids = set(baseline.target_ids), set(candidate.target_ids)
    left_samples, right_samples = set(baseline.sample_names), set(candidate.sample_names)
    excluded = {
        "baseline_only_guides": sorted(left_ids - right_ids),
        "candidate_only_guides": sorted(right_ids - left_ids),
        "baseline_only_samples": sorted(left_samples - right_samples),
        "candidate_only_samples": sorted(right_samples - left_samples),
    }
    if any(excluded.values()) and not shared_only:
        counts = ", ".join(f"{key}={len(value)}" for key, value in excluded.items() if value)
        raise ValueError(
            f"guide/sample sets differ ({counts}); check reference and biological sample names. "
            "Use --shared-only only for an explicitly partial comparison; missing counts are never zero-filled."
        )
    guides = [name for name in baseline.target_ids if name in right_ids]
    samples = [name for name in baseline.sample_names if name in right_samples]
    if not guides or not samples:
        raise ValueError("comparison needs at least one shared guide and one shared sample")
    return {
        "guides": guides, "samples": samples, "excluded": excluded,
        "complete": not any(excluded.values()),
        "baseline_input_guides": len(left_ids), "candidate_input_guides": len(right_ids),
        "baseline_input_samples": len(left_samples), "candidate_input_samples": len(right_samples),
    }


def _check_identity(baseline: CountTable, candidate: CountTable, guides: Sequence[str]) -> dict:
    left = {name.casefold(): name for name in baseline.metadata if name.casefold() in IDENTITY_FIELDS}
    right = {name.casefold(): name for name in candidate.metadata if name.casefold() in IDENTITY_FIELDS}
    shared = sorted(left.keys() & right.keys())
    li = {name: i for i, name in enumerate(baseline.target_ids)}
    ri = {name: i for i, name in enumerate(candidate.target_ids)}
    for field in shared:
        for guide in guides:
            if baseline.metadata[left[field]][li[guide]] != candidate.metadata[right[field]][ri[guide]]:
                raise ValueError(
                    f"conflicting {field} annotation for guide {guide!r}; "
                    "resolve reference identity before comparing counts"
                )
    return {
        "checked_fields": shared,
        "baseline_only_fields": sorted(left.keys() - right.keys()),
        "candidate_only_fields": sorted(right.keys() - left.keys()),
        "note": "Only identically named identity fields are checked; matching IDs do not prove matching biological references.",
    }


def _cells(baseline: CountTable, candidate: CountTable, alignment: dict) -> Iterator[tuple]:
    li = {name: i for i, name in enumerate(baseline.target_ids)}
    ri = {name: i for i, name in enumerate(candidate.target_ids)}
    ls = {name: i for i, name in enumerate(baseline.sample_names)}
    rs = {name: i for i, name in enumerate(candidate.sample_names)}
    for sample in alignment["samples"]:
        for guide in alignment["guides"]:
            before = baseline.counts[li[guide]][ls[sample]]
            after = candidate.counts[ri[guide]][rs[sample]]
            yield guide, sample, before, after, after - before


def compare_tables(baseline: CountTable, candidate: CountTable, *, shared_only: bool = False) -> dict:
    """Compare validated CountTables in linear time, with bounded HTML previews.

    IDs and samples, not row/column positions, establish correspondence. The
    complete changed-cell table is streamed separately; JSON contains at most
    TOP_CHANGES preview rows per sample. All count arithmetic uses Python ints.
    """
    alignment = _alignment(baseline, candidate, shared_only)
    identity = _check_identity(baseline, candidate, alignment["guides"])
    summaries, previews = {}, {}
    baseline_totals = dict(zip(baseline.sample_names, map(sum, zip(*baseline.counts))))
    candidate_totals = dict(zip(candidate.sample_names, map(sum, zip(*candidate.counts))))
    for sample in alignment["samples"]:
        summaries[sample] = {
            "sample": sample, "compared_guides": len(alignment["guides"]),
            "baseline_input_total": baseline_totals[sample],
            "candidate_input_total": candidate_totals[sample],
            "baseline_compared_total": 0, "candidate_compared_total": 0,
            "changed_guides": 0, "increased_guides": 0, "decreased_guides": 0,
            "newly_nonzero_guides": 0, "newly_zero_guides": 0,
            "absolute_count_delta": 0,
        }
        previews[sample] = []
    for guide, sample, before, after, delta in _cells(baseline, candidate, alignment):
        summary = summaries[sample]
        summary["baseline_compared_total"] += before
        summary["candidate_compared_total"] += after
        if not delta:
            continue
        summary["changed_guides"] += 1
        summary["increased_guides"] += int(delta > 0)
        summary["decreased_guides"] += int(delta < 0)
        summary["newly_nonzero_guides"] += int(before == 0 and after > 0)
        summary["newly_zero_guides"] += int(before > 0 and after == 0)
        summary["absolute_count_delta"] += abs(delta)
        item = (abs(delta), guide, before, after, delta)
        if len(previews[sample]) < TOP_CHANGES:
            heapq.heappush(previews[sample], item)
        else:
            heapq.heappushpop(previews[sample], item)
    for sample, summary in summaries.items():
        before, after = summary["baseline_compared_total"], summary["candidate_compared_total"]
        summary["total_delta"] = after - before
        summary["baseline_excluded_total"] = summary["baseline_input_total"] - before
        summary["candidate_excluded_total"] = summary["candidate_input_total"] - after
        summary["same_total_different_counts"] = before == after and summary["changed_guides"] > 0
        summary["top_changes"] = [
            {"guide": guide, "baseline": before, "candidate": after, "delta": delta}
            for _, guide, before, after, delta in sorted(previews[sample], reverse=True)
        ]
    return {
        "schema_version": SCHEMA_VERSION,
        "method": "raw_integer_counts_no_normalization",
        "scope": "complete" if alignment["complete"] else "shared_only_partial",
        "alignment": alignment, "identity": identity,
        "any_count_difference": any(row["changed_guides"] for row in summaries.values()),
        "samples": list(summaries.values()), "interpretation": BOUNDARY,
    }


def _report_html(report: dict) -> str:
    escape = lambda value: html.escape(str(value), quote=True)
    changed = sum(row["changed_guides"] for row in report["samples"])
    masked = sum(row["same_total_different_counts"] for row in report["samples"])
    complete = report["scope"] == "complete"
    title = "Totals agree. Guide counts do not." if masked else (
        "The count tables differ." if changed else "The compared counts agree."
    )
    coverage = "Complete guide and sample sets" if complete else "Partial comparison — exclusions must be reviewed"
    sections = []
    for row in report["samples"]:
        body = "".join(
            f'<tr><th scope="row">{escape(item["guide"])}</th>'
            f'<td>{item["baseline"]:,}</td><td>{item["candidate"]:,}</td>'
            f'<td>{item["delta"]:+,}</td></tr>' for item in row["top_changes"]
        )
        notice = '<p class="notice">Equal totals hide different per-guide counts in this sample.</p>' if row["same_total_different_counts"] else ""
        detail = (
            '<div class="scroll"><table><caption>Largest absolute count changes '
            f'(up to {TOP_CHANGES}; all changed cells are in changes.tsv)</caption>'
            '<thead><tr><th scope="col">Guide</th><th scope="col">Baseline</th>'
            '<th scope="col">Candidate</th><th scope="col">Delta</th></tr></thead>'
            f'<tbody>{body}</tbody></table></div>' if body else '<p>No changed counts in the compared guide set.</p>'
        )
        sections.append(
            f'<section><h2>{escape(row["sample"])}</h2>{notice}'
            '<dl class="metrics">'
            f'<div><dt>Baseline count total</dt><dd>{row["baseline_compared_total"]:,}</dd></div>'
            f'<div><dt>Candidate count total</dt><dd>{row["candidate_compared_total"]:,}</dd></div>'
            f'<div><dt>Changed guides</dt><dd>{row["changed_guides"]:,} / {row["compared_guides"]:,}</dd></div>'
            '</dl>'
            f'<p>{row["increased_guides"]:,} increased · {row["decreased_guides"]:,} decreased · '
            f'{row["newly_nonzero_guides"]:,} became nonzero · {row["newly_zero_guides"]:,} became zero.</p>'
            f'<p>Sum of absolute count deltas: {row["absolute_count_delta"]:,}. '
            'This is not a count of reassigned reads.</p>'
            f'<p>Count totals outside the compared guide set: baseline {row["baseline_excluded_total"]:,}; '
            f'candidate {row["candidate_excluded_total"]:,}.</p>{detail}</section>'
        )
    exclusions = "".join(
        f'<li>{escape(key.replace("_", " "))}: {len(value):,}</li>'
        for key, value in report["alignment"]["excluded"].items()
    )
    identities = ", ".join(report["identity"]["checked_fields"]) or "none available in both files"
    hashes = "".join(
        f'<dt>{escape(role)} input SHA-256</dt><dd class="hash">{escape(info["sha256"])}</dd>'
        for role, info in report["inputs"].items()
    )
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'unsafe-inline'; base-uri 'none'; form-action 'none'">
<title>DotMatch · Count comparison</title><style>
:root{{color-scheme:light;--ink:#172d39;--muted:#48606c;--line:#d8e0e4;--accent:#126455}}
*{{box-sizing:border-box}}body{{margin:0;background:#f6f8f9;color:var(--ink);font:17px/1.6 system-ui,sans-serif}}
main{{max-width:1100px;margin:auto;padding:48px 24px}}.brand{{font-weight:700;letter-spacing:.1em;color:var(--accent)}}
h1{{font-size:clamp(2rem,5vw,3.8rem);line-height:1.12;max-width:850px;letter-spacing:-.035em}}h2{{font-size:1.5rem;overflow-wrap:anywhere}}
p{{max-width:85ch}}.muted,dt,caption{{color:var(--muted)}}section{{background:white;padding:28px;margin:24px 0;border:1px solid var(--line);border-radius:12px}}
.metrics{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:20px}}.metrics dd{{font-size:1.8rem;font-weight:650;line-height:1.3}}
dd{{margin:4px 0 20px;overflow-wrap:anywhere}}.notice{{padding:12px 16px;border-left:4px solid var(--accent);background:#edf6f2}}
.scroll{{overflow-x:auto}}table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}}caption{{text-align:left;margin:12px 0}}
th,td{{padding:10px 12px;border-bottom:1px solid var(--line);text-align:right;overflow-wrap:anywhere}}th:first-child{{text-align:left}}
.hash{{font:13px/1.7 ui-monospace,monospace}}a{{color:var(--accent)}}a:focus-visible{{outline:3px solid var(--ink);outline-offset:4px}}
@media(max-width:650px){{main{{padding:24px 16px}}section{{padding:20px 16px}}.metrics{{grid-template-columns:1fr}}.metrics dd{{margin-bottom:4px}}}}
@media print{{body{{background:white}}main{{padding:0}}section{{break-inside:avoid;border-radius:0}}}}
</style></head><body><main><div class="brand">DOTMATCH / COUNT COMPARISON</div>
<h1>{escape(title)}</h1><p class="notice">{escape(coverage)}. {len(report["alignment"]["guides"]):,} guides · {len(report["samples"]):,} samples.</p>
<p>{escape(BOUNDARY)}</p><p class="muted">Baseline means your existing result, not ground truth. Counts are raw integers; nothing is normalized or zero-filled.</p>
{''.join(sections)}
<section><h2>What was compared</h2><ul>{exclusions}</ul>
<p>Checked identity fields: {escape(identities)}. Missing or differently named annotations are not verified.</p>
<p>Full exclusion lists, count summaries and bounded previews: <a href="report.json">report.json</a>. Every changed guide/sample cell: <a href="changes.tsv">changes.tsv</a>. File integrity: <a href="manifest.json">manifest.json</a>.</p>
<dl>{hashes}</dl><p>Private by default: this tool does not send data anywhere. The report still contains guide and sample identifiers and counts. Review it before sharing. Import TSV identifiers as text in spreadsheet software.</p>
</section></main></body></html>'''


def write_comparison(baseline_path: str | Path, candidate_path: str | Path,
                     out_dir: str | Path, *, shared_only: bool = False) -> dict:
    """Validate before writing. Never reuse an output directory or overwrite files.

    manifest.json is written last and marks completed output. An I/O failure can
    leave a partial new directory, but never removes/replaces pre-existing data.
    """
    baseline, left_source = _read_snapshot(baseline_path)
    candidate, right_source = _read_snapshot(candidate_path)
    report = compare_tables(baseline, candidate, shared_only=shared_only)
    report["inputs"] = {"baseline": left_source, "candidate": right_source}
    report["implementation_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    output = Path(out_dir)
    output.mkdir(parents=True, exist_ok=False)
    with (output / "changes.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["guide", "sample", "baseline", "candidate", "delta"])
        for cell in _cells(baseline, candidate, report["alignment"]):
            if cell[-1]:
                writer.writerow(cell)
    with (output / "report.json").open("x", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=True, allow_nan=False)
        handle.write("\n")
    with (output / "report.html").open("x", encoding="utf-8") as handle:
        handle.write(_report_html(report))
    files = {}
    for name in ("changes.tsv", "report.json", "report.html"):
        digest, size = hashlib.sha256(), 0
        with (output / name).open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
                size += len(chunk)
        files[name] = {"sha256": digest.hexdigest(), "bytes": size}
    with (output / "manifest.json").open("x", encoding="utf-8") as handle:
        json.dump({"schema_version": 1, "status": "complete", "files": files}, handle, indent=2)
        handle.write("\n")
    return report


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compare raw MAGeCK/DotMatch count tables locally before changing workflows. Agreement is not accuracy.")
    parser.add_argument("--baseline", required=True, help="existing raw-count TSV or .tsv.gz; not ground truth")
    parser.add_argument("--candidate", required=True, help="candidate raw-count TSV or .tsv.gz")
    parser.add_argument("--out-dir", required=True, help="new directory for HTML, JSON, TSV and integrity manifest")
    parser.add_argument("--shared-only", action="store_true", help="explicitly compare only shared guides/samples and report all exclusions")
    parser.add_argument("--fail-on-difference", action="store_true", help="exit 1 on changed counts OR excluded guides/samples, after writing the report")
    args = parser.parse_args(argv)
    try:
        report = write_comparison(args.baseline, args.candidate, args.out_dir, shared_only=args.shared_only)
    except (ValueError, OSError, csv.Error, UnicodeError, EOFError) as exc:
        parser.exit(2, f"count comparison: {exc}\n")
    changed = sum(row["changed_guides"] for row in report["samples"])
    print(f"Compared {len(report['alignment']['guides'])} guides across {len(report['samples'])} samples: {changed} changed guide/sample cells.")
    if report["scope"] != "complete":
        print("PARTIAL comparison: inspect excluded guides/samples in report.json.")
    print(f"Report: {Path(args.out_dir) / 'report.html'}")
    print("Agreement is not biological accuracy. No data was uploaded.")
    return int(args.fail_on_difference and (report["any_count_difference"] or report["scope"] != "complete"))


if __name__ == "__main__":
    raise SystemExit(main())
