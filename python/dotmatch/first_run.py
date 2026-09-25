"""Installed first-run experience over the existing engines, without telemetry."""
from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import html
import json
import os
import shlex
import sys
import tempfile
from importlib import resources
from pathlib import Path
from typing import Sequence

from .assayspec import AssaySpecError, FASTQ_SUFFIXES, command_assay, scaffold_assay_project
from .count_compare import write_comparison
from .count_io import read_count_table
from .sensitivity import MODES, run_sensitivity


def resolve_fastqs(patterns: Sequence[str]) -> list[Path]:
    """Resolve every supplied argument before writing or sampling any inputs."""
    sources: list[Path] = []
    seen: set[Path] = set()
    basenames: set[str] = set()
    for pattern in patterns:
        expanded = os.path.expanduser(pattern)
        literal = Path(expanded)
        # A real path with brackets in its name is not a glob expression.
        matches = [literal] if os.path.lexists(literal) else [Path(p) for p in sorted(glob.glob(expanded))]
        if not matches:
            raise ValueError(f"--fastq matched no files: {pattern!r}; no project was created")
        for match in matches:
            source = match.resolve()
            if not source.is_file():
                raise ValueError(f"FASTQ is not a regular file: {match}")
            if not any(source.name.lower().endswith(suffix) for suffix in FASTQ_SUFFIXES):
                raise ValueError(f"input does not look like a FASTQ file: {source}")
            if source in seen:
                raise ValueError(f"duplicate FASTQ input: {source}; remove repeated or overlapping patterns")
            if source.name.casefold() in basenames:
                raise ValueError(f"duplicate FASTQ basename: {source.name}; use distinct filenames")
            seen.add(source)
            basenames.add(source.name.casefold())
            sources.append(source)
    if not sources:
        raise ValueError("at least one FASTQ input is required")
    return sources


def _positive(text: str) -> int:
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return value


def _nonnegative(text: str) -> int:
    value = int(text)
    if value < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return value


def quickstart_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="dotmatch crispr quickstart", description=(
        "Prepare a reviewable CRISPR project. Every FASTQ argument must resolve. "
        "Inputs are copied by default; --link-reads keeps large FASTQs in place."))
    parser.add_argument("--library", required=True, help="guide library CSV/TSV")
    parser.add_argument("--fastq", action="append", required=True, help="FASTQ path or quoted glob; repeat for multiple inputs")
    parser.add_argument("--out", required=True, help="new project directory; existing paths are refused")
    parser.add_argument("--link-reads", action="store_true", help="symlink original FASTQs instead of copying; originals must stay available and unchanged")
    parser.add_argument("--threads", type=_positive, default=1)
    parser.add_argument("--max-reads", type=_positive, default=50000)
    parser.add_argument("--max-start", type=_nonnegative, default=32)
    parser.add_argument("--no-run", action="store_true", help="leave the project in draft, even with --accept-inference")
    parser.add_argument("--accept-inference", action="store_true", help="run if inference is ready; uncertain inference still requires manual review")
    args = parser.parse_args(argv)
    try:
        library = Path(args.library).expanduser().resolve()
        if not library.is_file():
            raise ValueError(f"library is not a regular file: {library}")
        sources = resolve_fastqs(args.fastq)
        project = Path(args.out).expanduser().absolute()
        if os.path.lexists(project):
            raise FileExistsError(f"refusing to overwrite existing project: {project}")
        project.parent.mkdir(parents=True, exist_ok=True)
        # Reserve exclusively, before the scaffold checks the empty directory.
        project.mkdir(exist_ok=False)
        with tempfile.TemporaryDirectory(prefix=f".{project.name}.dotmatch-inputs-", dir=project.parent) as temporary:
            staging = Path(temporary)
            for source in sources:
                (staging / source.name).symlink_to(source)
            result = scaffold_assay_project(
                template="crispr", project_dir=project, reads_dir=staging,
                targets=library, link_reads=args.link_reads, threads=args.threads,
                max_reads=args.max_reads, max_start=args.max_start,
                force_draft=args.no_run or not args.accept_inference,
            )
        report_path = project / "inference_report.json"
        report = json.loads(report_path.read_text(encoding="utf-8"))
        original = {source.name: str(source) for source in sources}
        # Keep source identity useful after disposable staging has disappeared.
        for sample in report["samples"]:
            sample["source_fastq"] = original[Path(sample["fastq"]).name]
        report["input_storage"] = "linked" if args.link_reads else "copied"
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"Created reviewable CRISPR project: {result['project']}")
        print(f"Inputs: {len(sources)} FASTQ files; {'linked to originals' if args.link_reads else 'copied into project'}.")
        if args.link_reads:
            print("Keep the original FASTQs available and unchanged. This project is not self-contained.")
        print(f"Review: {report_path}")
        print("Check sample names against your run sheet: filenames are not biological replicate declarations.")
        if args.no_run or not args.accept_inference:
            print('After confirming the window, orientation and samples, change status = "draft" to "ready" '
                  f'in {project / "assay.toml"}.')
            print("Run: " + shlex.join(["dotmatch", "assay", "start", str(project / "assay.toml")]))
            return 0
        if report["status"] != "ready":
            print("Inference still needs review; no run was started. Inspect inference_report.json and assay.toml.", file=sys.stderr)
            return 2
        return command_assay(["start", str(project / "assay.toml")])
    except (AssaySpecError, ValueError, OSError) as exc:
        print(f"dotmatch crispr: {exc}", file=sys.stderr)
        return 2


def _hash_file(path: Path) -> dict:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            size += len(chunk)
            digest.update(chunk)
    return {"sha256": digest.hexdigest(), "bytes": size}


def _demo_html() -> str:
    from . import __version__
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'unsafe-inline'; base-uri 'none'; form-action 'none'">
<title>DotMatch · First run</title><style>
*{{box-sizing:border-box}}body{{margin:0;color:#182c36;background:#f8faf9;font:17px/1.65 system-ui,sans-serif}}
main{{max-width:960px;margin:auto;padding:56px 24px}}h1{{font-size:clamp(2.2rem,6vw,3.8rem);line-height:1.12;letter-spacing:-.04em;max-width:760px}}
h2{{font-size:1.35rem;margin:0 0 12px}}.kicker{{font-size:.85rem;letter-spacing:.12em;color:#3c6257}}p{{max-width:78ch}}
section{{background:white;border:1px solid #d7e2df;border-radius:12px;padding:28px;margin:24px 0}}
a{{color:#145e4b;text-underline-offset:4px}}a:focus-visible{{outline:3px solid #145e4b;outline-offset:4px}}
pre{{white-space:pre-wrap;overflow-wrap:anywhere;background:#f0f5f3;padding:16px;font:14px/1.7 ui-monospace,monospace}}
.note{{border-left:3px solid #466a5c;padding-left:20px}}.detail{{color:#4e636a}}strong{{font-weight:650}}
@media(max-width:600px){{main{{padding:28px 16px}}section{{padding:20px}}}}
</style></head><body><main><p class="kicker">DOTMATCH / SYNTHETIC FIRST RUN</p>
<h1>Your installation produced the expected counts.</h1>
<p>Nine synthetic reads, five target entries, three explicit assignment rules.
The native engine ran locally and its count tables and read outcomes matched the checked fixture.</p>
<p class="note">This checks software behaviour, not biological accuracy, clinical suitability or independent adoption. No research data was downloaded or uploaded.</p>
<section><h2>1. Inspect the assignments</h2><p>Exact matching and radius-one each assign three reads here, but to different guides. Best-distance assigns five. A larger assigned total does not tell you which policy is appropriate.</p>
<p><a href="sensitivity/report.html">Open the interactive assignment review</a></p>
<p class="detail">Inspect ambiguous reads, the short read, duplicate targets, and the literal N. N is not treated as a wildcard.</p></section>
<section><h2>2. Compare the count tables</h2><p>The exact and radius-one totals agree, while three guide counts differ. Reordered rows are matched by guide identifier.</p>
<p><a href="comparison/report.html">Open the count comparison</a></p>
<p class="detail">The complete tables remain beside the reports. <a href="verification.json">Verification record</a> · <a href="manifest.json">File checksums</a></p></section>
<section><h2>3. Prepare your own reads</h2><p>Use only the FASTQ read containing the guide window. Paired-end files and technical lanes are not automatically identified as biological samples.</p>
<pre>dotmatch crispr quickstart \\
  --library guides.csv \\
  --fastq 'fastqs/*_R1*.fastq.gz' \\
  --link-reads --out crispr-screen/</pre>
<p>Replace the paths with your files. Linked reads stay in place and must remain available and unchanged. Review the proposed window, orientation, samples and library before changing the project from draft to ready.</p>
<pre>dotmatch assay start crispr-screen/assay.toml
# After reviewing a completed run:
dotmatch assay handoff crispr-screen/assay.toml</pre>
<p>Keep your existing pipeline while evaluating. Compare its raw counts with the completed DotMatch output:</p>
<pre>dotmatch compare-counts \\
  --baseline existing-counts.tsv \\
  --candidate crispr-screen/assay_out/counts.mageck.tsv \\
  --out-dir comparison/</pre>
<p>Agreement is not accuracy. Differences require inspection, not automatic selection of the result with more counts.</p></section>
<p class="detail">DotMatch {html.escape(__version__)}. This report is fully local. The bundled inputs in <code>inputs/</code> are synthetic and may be inspected and shared. Your own reports may contain private sample identifiers and counts.</p>
</main></body></html>'''


def run_demo(out_dir: str | Path) -> dict:
    """Run the committed nine-read fixture through real engines, not mock outputs."""
    fixture_bytes = resources.files("dotmatch").joinpath("data", "first-run.json").read_bytes()
    fixture = json.loads(fixture_bytes)
    if fixture.get("schema_version") != 1 or fixture.get("synthetic") is not True:
        raise ValueError("unsupported first-run fixture")
    output = Path(out_dir).expanduser().absolute()
    output.mkdir(parents=True, exist_ok=False)
    inputs = output / "inputs"
    inputs.mkdir()
    library, reads = inputs / "targets.tsv", inputs / "reads.fastq"
    library.write_text(fixture["targets_tsv"], encoding="utf-8")
    reads.write_text(fixture["reads_fastq"], encoding="utf-8")
    summary = run_sensitivity(targets=library, reads=reads, target_start=0, target_length=20,
                              sample_label="synthetic", out_dir=output / "sensitivity", write_read_changes=True)
    expected = fixture["expected"]
    if summary["read_count"] != expected["read_count"] or summary["outcomes"] != expected["outcomes"]:
        raise ValueError("first-run read outcomes did not match the checked fixture; no completed demo manifest was written")
    for policy in MODES:
        table = read_count_table(output / "sensitivity" / f"{policy}.counts.tsv")
        actual = {name: row[0] for name, row in zip(table.target_ids, table.counts)}
        if table.sample_names != ("synthetic",) or actual != expected["counts"][policy]:
            raise ValueError(f"{policy} counts did not match the checked fixture; no completed demo manifest was written")
    comparison = write_comparison(output / "sensitivity/exact.counts.tsv",
                                  output / "sensitivity/radius_k1.counts.tsv", output / "comparison")
    if comparison["scope"] != "complete" or not comparison["samples"][0]["same_total_different_counts"]:
        raise ValueError("first-run comparison did not expose the expected equal-total discrepancy")
    verification = {
        "schema_version": 1, "synthetic": True, "status": "expected_results_verified",
        "fixture_sha256": hashlib.sha256(fixture_bytes).hexdigest(),
        "software_version": summary["software_version"],
        "implementation_sha256": _hash_file(Path(__file__))["sha256"],
        "native_library_sha256": summary["native_library_sha256"],
        "expected": expected,
        "boundary": "Software mechanics only, not biological accuracy or evidence of independent adoption.",
    }
    (output / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "index.html").write_text(_demo_html(), encoding="utf-8")
    files = {p.relative_to(output).as_posix(): _hash_file(p) for p in sorted(output.rglob("*")) if p.is_file()}
    # Written last. Failed runs retain their diagnostic outputs but never this marker.
    with (output / "manifest.json").open("x", encoding="utf-8") as handle:
        json.dump({"schema_version": 1, "status": "complete", "synthetic": True, "files": files}, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return verification


def demo_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="dotmatch demo", description="Run the bundled synthetic first-run example locally; no clone or download needed.")
    parser.add_argument("--out-dir", required=True, help="new directory for inputs, reports, verification and checksums")
    args = parser.parse_args(argv)
    try:
        run_demo(args.out_dir)
    except (ValueError, OSError, csv.Error, RuntimeError) as exc:
        print(f"dotmatch demo: {exc}", file=sys.stderr)
        return 2
    print("Synthetic first run verified: 9 reads, 5 target entries, 3 assignment policies.")
    print(f"Open: {Path(args.out_dir) / 'index.html'}")
    print("No data was uploaded. This verifies software behaviour, not biological accuracy.")
    return 0
