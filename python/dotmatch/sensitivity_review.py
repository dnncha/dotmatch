"""Offline, read-only review of DotMatch sensitivity artifacts (standard library).

Rendering is not a new analysis or biological validation. Counts and transitions
are reconciled before presentation; recorded hashes are not authenticated keys.
"""
from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import html
import io
import json
import re
import sys
from pathlib import Path

if __package__:
    from .sensitivity_review_assets import CSS, JS, SHELL
else:  # Source-only report rebuilding does not load DotMatch's native engine.
    from sensitivity_review_assets import CSS, JS, SHELL

MODES = ("exact", "radius_k1", "best_k1")
STATES = ("unique", "ambiguous", "none", "invalid")
PAIRS = (("exact", "radius_k1"), ("exact", "best_k1"), ("radius_k1", "best_k1"))
MAX_INTEGER = 2**53 - 1
MAX_GUIDES = 250_000
MAX_TABLE_BYTES = 64 * 1024 * 1024
MAX_MANIFEST_BYTES = 1024 * 1024
GUIDE_FIELDS = ["target_id", "gene", *MODES, "radius_minus_exact", "best_minus_exact"]
TRANSITION_FIELDS = ["from_policy", "to_policy", "from_status", "to_status", "reads"]


class ReviewCapacityError(ValueError):
    """Viewer capacity exceeded; preserve the producer's full static/TSV outputs."""


def _integer(value, name: str, minimum: int = 0) -> int:
    if type(value) is not int or not minimum <= value <= MAX_INTEGER:
        raise ValueError(f"{name}: expected an integer in {minimum}..{MAX_INTEGER}")
    return value


def _text(value, name: str, *, empty: bool = False) -> str:
    if not isinstance(value, str) or (not empty and not value):
        raise ValueError(f"{name}: missing text")
    if len(value) > 16_384:
        raise ReviewCapacityError(f"{name}: exceeds the portable viewer text limit")
    if any(ord(c) < 32 and c not in "\t\n\r" for c in value):
        raise ValueError(f"{name}: unsupported control character")
    return value


def _number(value: str, name: str, *, signed: bool = False) -> int:
    if not re.fullmatch(r"-?(?:0|[1-9][0-9]*)" if signed else r"(?:0|[1-9][0-9]*)", value):
        raise ValueError(f"{name}: invalid decimal count")
    n = int(value)
    if abs(n) > MAX_INTEGER:
        raise ValueError(f"{name}: exceeds exact browser integer range")
    return n


def _hash(value, name: str) -> str:
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ValueError(f"{name}: invalid SHA-256")
    return value


def _no_duplicate_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _json(text: str):
    def invalid(value):
        raise ValueError(f"Non-finite JSON value: {value}")
    return json.loads(text, object_pairs_hook=_no_duplicate_keys, parse_constant=invalid)


def _read(path: Path, limit: int) -> bytes:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"{path.name}: expected a regular, non-symlink file")
    with path.open("rb") as handle:
        data = handle.read(limit + 1)
    if len(data) > limit:
        raise ReviewCapacityError(f"{path.name}: exceeds this viewer's {limit:,}-byte limit; TSVs remain available")
    return data


def _table(data: bytes, fields: list[str], name: str, limit: int):
    # csv's process-wide field limit is intentionally left unchanged.
    reader = csv.reader(io.StringIO(data.decode("utf-8-sig"), newline=""), delimiter="\t", strict=True)
    if next(reader, None) != fields:
        raise ValueError(f"{name}: unsupported or duplicate column headers")
    for i, row in enumerate(reader, start=1):
        if i > limit or len(row) != len(fields):
            raise ValueError(f"{name}: too many rows or malformed row {i}")
        yield dict(zip(fields, row))


def _source(directory: Path, name: str, manifest: dict, limit: int):
    data = _read(directory / name, limit)
    info = {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}
    if "artifacts" in manifest:
        recorded = manifest["artifacts"].get(name)
        if not isinstance(recorded, dict):
            raise ValueError(f"{name}: missing artifact identity")
        _hash(recorded.get("sha256"), name)
        _integer(recorded.get("bytes"), name + " bytes")
        if recorded != info:
            raise ValueError(f"{name}: bytes/hash do not match the supplied completion manifest")
    return data, info


def build_review_data(summary: dict, directory: str | Path, *, staged: bool = False) -> dict:
    """Read the exact v1 artifacts, reject contradictions, copy no raw reads.

    Called inside the producer's staging directory before publication, or on a
    completed bundle. In the latter case source hashes must match summary.json.
    Size limits bound the portable viewer, not the scientific count exports.
    """
    directory = Path(directory)
    if not isinstance(summary, dict) or summary.get("schema_version") != "dotmatch.sensitivity.v1":
        raise ValueError("Unsupported sensitivity schema; expected dotmatch.sensitivity.v1")
    if summary.get("completion") != "complete":
        raise ValueError("Incomplete sensitivity run: no completed review can be rendered")
    if "artifacts" in summary and not isinstance(summary["artifacts"], dict):
        raise ValueError("Invalid artifact manifest")
    reads = _integer(summary.get("read_count"), "read_count", 1)
    targets = _integer(summary.get("target_count"), "target_count", 1)
    if targets > MAX_GUIDES:
        raise ReviewCapacityError(f"Portable review limit is {MAX_GUIDES:,} guides; count/TSV analysis is not truncated")
    changed = _integer(summary.get("changed_reads"), "changed_reads")
    if changed > reads:
        raise ValueError("changed_reads exceeds read_count")
    _text(summary.get("sample_label"), "sample_label")
    params = summary.get("parameters", {})
    if not isinstance(params, dict) or params.get("metric") != "hamming" or params.get("orientation") != "as_supplied":
        raise ValueError("Unsupported metric or orientation; this review requires fixed-window Hamming")
    start = _integer(params.get("target_start"), "target_start")
    length = _integer(params.get("target_length"), "target_length", 1)
    _integer(start + length, "window end")
    outcomes = summary.get("outcomes", {})
    if not isinstance(outcomes, dict) or set(outcomes) != set(MODES):
        raise ValueError("Missing or unknown outcome policy")
    for mode in MODES:
        if not isinstance(outcomes[mode], dict) or set(outcomes[mode]) != set(STATES):
            raise ValueError(f"{mode}: missing or unknown outcome state")
        if sum(_integer(outcomes[mode][state], f"{mode}.{state}") for state in STATES) != reads:
            raise ValueError(f"{mode}: outcomes do not conserve reads")
    for name in ("implementation_sha256", "native_library_sha256"):
        _hash(summary.get(name), name)
    inputs = summary.get("inputs", {})
    if not isinstance(inputs, dict):
        raise ValueError("Missing declared input identities")
    for name in ("targets", "reads"):
        entry = inputs.get(name)
        if not isinstance(entry, dict):
            raise ValueError(f"Missing {name} input identity")
        _text(entry.get("name"), name)
        _hash(entry.get("sha256"), name)

    data, guide_info = _source(directory, "guide_deltas.tsv", summary, MAX_TABLE_BYTES)
    guides, seen, sums = [], set(), [0, 0, 0]
    for row in _table(data, GUIDE_FIELDS, "guide_deltas.tsv", MAX_GUIDES):
        key = _text(row["target_id"], "target_id")
        if key in seen:
            raise ValueError(f"Duplicate target ID: {key}")
        seen.add(key)
        values = [_number(row[mode], mode) for mode in MODES]
        if (_number(row["radius_minus_exact"], "radius delta", signed=True) != values[1] - values[0]
                or _number(row["best_minus_exact"], "best delta", signed=True) != values[2] - values[0]):
            raise ValueError("Guide deltas contradict recorded counts")
        guides.append([key, _text(row["gene"], "gene", empty=True), *values])
        sums = [a + b for a, b in zip(sums, values)]
    if len(guides) != targets or any(sums[i] != outcomes[mode]["unique"] for i, mode in enumerate(MODES)):
        raise ValueError("Guide counts do not reconcile with target_count / unique outcomes")
    comparisons = summary.get("count_comparisons")
    if not isinstance(comparisons, list) or len(comparisons) != 3:
        raise ValueError("Expected all three count comparisons")
    indexed = {}
    for item in comparisons:
        if not isinstance(item, dict):
            raise ValueError("Malformed count comparison")
        pair = (item.get("left"), item.get("right"))
        if pair not in PAIRS or pair in indexed:
            raise ValueError("Unknown or duplicate policy comparison")
        indexed[pair] = item
    for left, right in PAIRS:
        a, b = MODES.index(left) + 2, MODES.index(right) + 2
        differing = sum(g[a] != g[b] for g in guides)
        item = indexed[(left, right)]
        if (type(item.get("counts_identical")) is not bool or item["counts_identical"] != (differing == 0)
                or _integer(item.get("differing_guides"), "differing_guides") != differing
                or type(item.get("total_count_delta")) is not int
                or item["total_count_delta"] != sums[b - 2] - sums[a - 2]):
            raise ValueError("Recorded count comparison contradicts guide counts")

    data, transition_info = _source(directory, "transitions.tsv", summary, MAX_TABLE_BYTES)
    transitions, matrix = [], {}
    for row in _table(data, TRANSITION_FIELDS, "transitions.tsv", 48):
        pair = (row["from_policy"], row["to_policy"])
        a, b = row["from_status"], row["to_status"]
        key = (*pair, a, b)
        if pair not in PAIRS or a not in STATES or b not in STATES or key in matrix:
            raise ValueError("Unknown or duplicate read-state transition")
        n = _number(row["reads"], "transition reads")
        matrix[key] = n
        transitions.append([*key, n])
    if len(matrix) != 48:
        raise ValueError("Expected all 48 transitions, including zero cells")
    for left, right in PAIRS:
        for state in STATES:
            if (sum(matrix[left, right, state, b] for b in STATES) != outcomes[left][state]
                    or sum(matrix[left, right, a, state] for a in STATES) != outcomes[right][state]):
                raise ValueError("Transition margins contradict the policy outcomes")
        off_diagonal = sum(matrix[left, right, a, b] for a in STATES for b in STATES if a != b)
        if off_diagonal > changed:
            raise ValueError("Read-state transitions exceed changed_reads")
    sources = {"guide_deltas.tsv": guide_info, "transitions.tsv": transition_info}
    # The optional read-ID artifact is hashed but NEVER embedded in the HTML.
    changes_path = directory / "read_changes.tsv"
    if "artifacts" in summary and "read_changes.tsv" in summary["artifacts"]:
        info = summary["artifacts"]["read_changes.tsv"]
        if not isinstance(info, dict):
            raise ValueError("Invalid read_changes artifact declaration")
        sources["read_changes.tsv"] = {"sha256": _hash(info.get("sha256"), "read_changes.tsv"),
                                       "bytes": _integer(info.get("bytes"), "read_changes bytes")}
    elif changes_path.exists():
        if changes_path.is_symlink() or not changes_path.is_file():
            raise ValueError("read_changes.tsv must be a regular file")
        digest, size = hashlib.sha256(), 0
        with changes_path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
        sources["read_changes.tsv"] = {"sha256": digest.hexdigest(), "bytes": size}
    # A small explicit projection avoids copying arbitrary extension fields.
    snapshot = {key: summary[key] for key in ("schema_version", "completion", "sample_label", "read_count", "target_count",
                "changed_reads", "outcomes", "count_comparisons", "implementation_sha256", "native_library_sha256")}
    snapshot["parameters"] = {key: params[key] for key in ("target_start", "target_length", "metric", "orientation")}
    snapshot["inputs"] = {name: {"name": inputs[name]["name"], "sha256": inputs[name]["sha256"]} for name in ("targets", "reads")}
    snapshot["software_version"] = _text(summary.get("software_version"), "software_version")
    assets_path = Path(__file__).with_name("sensitivity_review_assets.py")
    renderer_hash = hashlib.sha256(Path(__file__).read_bytes() + b"\0" + assets_path.read_bytes()).hexdigest()
    return {"schema_version": "dotmatch.review.v1", "summary": snapshot, "guides": guides,
            "transitions": transitions, "sources": sources, "renderer_sha256": renderer_hash,
            "source_manifest_checked": "artifacts" in summary and not staged}


def render_sensitivity_report(summary: dict, directory: str | Path, *, staged: bool = False) -> str:
    """Return one network-free HTML document with a useful no-JavaScript fallback."""
    data = build_review_data(summary, directory, staged=staged)
    payload = json.dumps(data, ensure_ascii=True, separators=(",", ":"), allow_nan=False).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    script_hash = base64.b64encode(hashlib.sha256(JS.encode()).digest()).decode()
    csp = f"default-src 'none'; script-src 'sha256-{script_hash}'; style-src 'unsafe-inline'; connect-src 'none'; img-src 'none'; font-src 'none'; base-uri 'none'; form-action 'none'; object-src 'none'"
    fallback = '<h2>Recorded outcomes</h2><table><caption>All reads under each policy</caption><thead><tr><th>Policy</th>' + ''.join(f'<th>{s}</th>' for s in STATES) + '</tr></thead><tbody>'
    for mode in MODES:
        fallback += f'<tr><th>{mode}</th>' + ''.join(f'<td>{summary["outcomes"][mode][s]:,}</td>' for s in STATES) + '</tr>'
    fallback += '</tbody></table><h2>Guide counts (first 50 of ' + str(len(data['guides'])) + ')</h2><table><thead><tr><th>Target</th><th>Gene</th><th>Exact</th><th>Radius k=1</th><th>Best k=1</th></tr></thead><tbody>'
    for row in data['guides'][:50]:
        fallback += '<tr>' + ''.join('<td>' + html.escape(str(v)) + '</td>' for v in row) + '</tr>'
    fallback += '</tbody></table><p>JavaScript enables full-library search, policy comparison and read-state review. The TSV artifacts retain every target.</p>'
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<meta http-equiv="Content-Security-Policy" content="{html.escape(csp, quote=True)}">'
            '<meta name="color-scheme" content="light"><title>DotMatch · Assignment review</title>'
            f'<style>{CSS}</style></head><body>{SHELL}'
            f'<section id="fallback" class="fallback"><h1>DotMatch · {html.escape(summary["sample_label"])}</h1>{fallback}</section>'
            f'<script type="application/json" id="review-data">{payload}</script>'
            f'<script type="text/javascript">{JS}</script></body></html>')


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True, help="Completed sensitivity output directory")
    parser.add_argument("--out", type=Path, required=True, help="New HTML file; never overwritten")
    args = parser.parse_args(argv)
    try:
        summary = _json(_read(args.bundle / "summary.json", MAX_MANIFEST_BYTES).decode("utf-8"))
        if not isinstance(summary, dict) or not isinstance(summary.get("artifacts"), dict):
            raise ValueError("A completed artifact manifest is required to rebuild a review")
        report = render_sensitivity_report(summary, args.bundle)
        with args.out.open("x", encoding="utf-8", newline="\n") as handle:
            handle.write(report)
    except (OSError, ValueError, TypeError, KeyError, csv.Error) as exc:
        print(f"dotmatch review: {exc}", file=sys.stderr)
        return 2
    print(str(args.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
