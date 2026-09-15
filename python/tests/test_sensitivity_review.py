"""Portable renderer contract tests; independent of native engine availability.

Fixtures enumerate the documented synthetic reads, not real biological data.
The existing test_sensitivity.py remains the native producer integration gate.
"""
import copy
import csv
import hashlib
import importlib.util
import io
import json
import re
import sys
from collections import Counter
from pathlib import Path

import pytest

SOURCE = Path(__file__).resolve().parents[1] / "dotmatch"
# Loading by path keeps these report-only tests independent of native libraries.
for name in ("sensitivity_review_assets", "sensitivity_review"):
    spec = importlib.util.spec_from_file_location(name, SOURCE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
review = sys.modules["sensitivity_review"]
MODES, STATES, PAIRS = review.MODES, review.STATES, review.PAIRS


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def fixture(directory, extra_zero_guides=0, identifiers=None, read_ids=None):
    """Reference-enumerated synthetic fixture, NOT a native-engine execution."""
    directory.mkdir(parents=True, exist_ok=True)
    sequences = ["ACGTACGTACGTACGTACGT", "CCGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCA",
                 "GGTAGGTAGGTAGGTAGGTA", "GGTAGGTAGGTAGGTAGGTA"]
    ids = identifiers or ["guide_A", "guide_B", "guide_C", "guide_D", "guide_D_duplicate"]
    reads = [sequences[0], sequences[1], "GCGTACGTACGTACGTACGT", sequences[2], "CGCATGCATGCATGCATGCA",
             "GGGGGGGGGGGGGGGGGGGG", "AC", sequences[3], "NGCATGCATGCATGCATGCA"]
    read_ids = read_ids or [f"synthetic_record_{i + 1}" for i in range(9)]
    counts = {m: Counter() for m in MODES}
    outcomes = {m: Counter({s: 0 for s in STATES}) for m in MODES}
    transitions = {p: Counter() for p in PAIRS}
    changes = []
    for ordinal, sequence in enumerate(reads, 1):
        calls = {}
        for mode in MODES:
            if len(sequence) != 20:
                call = ("invalid", "")
            else:
                candidates = [(sum(a != b for a, b in zip(sequence, target)), i) for i, target in enumerate(sequences)]
                candidates = [(n, i) for n, i in candidates if n <= (0 if mode == "exact" else 1)]
                if mode == "best_k1" and candidates:
                    minimum = min(n for n, _ in candidates)
                    candidates = [(n, i) for n, i in candidates if n == minimum]
                call = ("none", "") if not candidates else ("ambiguous", "") if len(candidates) > 1 else ("unique", ids[candidates[0][1]])
            calls[mode] = call
            outcomes[mode][call[0]] += 1
            if call[0] == "unique":
                counts[mode][call[1]] += 1
        for pair in PAIRS:
            transitions[pair][calls[pair[0]][0], calls[pair[1]][0]] += 1
        if len(set(calls.values())) > 1:
            changes.append([ordinal, read_ids[ordinal-1], *(value for m in MODES for value in calls[m])])
    ids = ids + [f"zero_{i:06d}" for i in range(extra_zero_guides)]
    guide_rows = [[key, "GENE_" + key, *(counts[m][key] for m in MODES), counts["radius_k1"][key]-counts["exact"][key], counts["best_k1"][key]-counts["exact"][key]] for key in ids]
    write_tsv(directory / "guide_deltas.tsv", review.GUIDE_FIELDS, guide_rows)
    write_tsv(directory / "transitions.tsv", review.TRANSITION_FIELDS, [[*p, a, b, transitions[p][a,b]] for p in PAIRS for a in STATES for b in STATES])
    write_tsv(directory / "read_changes.tsv", ["record_index", "read_id", *(f"{m}_{k}" for m in MODES for k in ("status", "target_id"))], changes)
    summary = {"schema_version": "dotmatch.sensitivity.v1", "completion": "complete", "software_version": "synthetic-reference-fixture",
               "sample_label": "SYNTHETIC nine-read example (reference enumeration, not a native run)", "read_count": 9,
               "target_count": len(ids), "changed_reads": len(changes), "outcomes": outcomes,
               "parameters": {"target_start": 0, "target_length": 20, "metric": "hamming", "orientation": "as_supplied"},
               "inputs": {"targets": {"name": "synthetic targets", "sha256": sha("\n".join(sequences).encode())}, "reads": {"name": "synthetic reads", "sha256": sha("\n".join(reads).encode())}},
               "implementation_sha256": sha(Path(__file__).read_bytes()),
               "native_library_sha256": sha(b"TEST FIXTURE: NO NATIVE LIBRARY WAS EXECUTED"),
               "count_comparisons": [{"left": a, "right": b, "counts_identical": all(counts[a][key] == counts[b][key] for key in ids),
                                      "differing_guides": sum(counts[a][key] != counts[b][key] for key in ids),
                                      "total_count_delta": sum(counts[b].values())-sum(counts[a].values())} for a,b in PAIRS]}
    return summary


def seal(directory, summary):
    summary["artifacts"] = {p.name: {"bytes": p.stat().st_size, "sha256": sha(p.read_bytes())} for p in directory.iterdir() if p.is_file() and p.name != "summary.json"}
    (directory / "summary.json").write_text(json.dumps(summary))
    return summary


def test_reconciles_documented_fixture(tmp_path):
    s = fixture(tmp_path)
    data = review.build_review_data(s, tmp_path)
    assert s["outcomes"]["exact"] == dict(zip(STATES, (3,1,4,1)))
    assert s["outcomes"]["radius_k1"] == dict(zip(STATES, (3,4,1,1)))
    assert s["outcomes"]["best_k1"] == dict(zip(STATES, (5,2,1,1)))
    assert s["changed_reads"] == 5
    assert len(data["guides"]) == 5 and len(data["transitions"]) == 48
    assert data["summary"]["count_comparisons"][0]["total_count_delta"] == 0
    assert data["summary"]["count_comparisons"][0]["differing_guides"] == 3


def test_all_guides_including_zero_beyond_first_page(tmp_path):
    s = fixture(tmp_path, extra_zero_guides=110)
    document = review.render_sensitivity_report(s, tmp_path)
    match = re.search(r'<script type="application/json" id="review-data">(.*?)</script>', document, re.S)
    assert len(json.loads(match[1])["guides"]) == 115
    assert "zero_000109" in document


def test_html_escapes_and_embeds_no_read_ids(tmp_path):
    attack = '</script><script>alert("x")</script>&\u2028'
    s = fixture(tmp_path, identifiers=[attack,"b","c","d","e"], read_ids=["PRIVATE_READ_ID_"+str(i) for i in range(9)])
    document = review.render_sensitivity_report(s, tmp_path)
    assert attack not in document and "&lt;script&gt;" in document
    assert '<script>' not in document
    assert "PRIVATE_READ_ID" not in document
    assert "ACGTACGTACGTACGTACGT" not in document
    assert 'connect-src &#x27;none&#x27;' in document
    assert "script-src &#x27;sha256-" in document
    assert "fetch(" not in document and "XMLHttpRequest" not in document


def test_csp_hash_matches_actual_script(tmp_path):
    import base64
    document = review.render_sensitivity_report(fixture(tmp_path), tmp_path)
    javascript = re.search(r'<script type="text/javascript">(.*?)</script>', document, re.S)[1]
    assert base64.b64encode(hashlib.sha256(javascript.encode()).digest()).decode() in document


@pytest.mark.parametrize("change", [
    lambda s:s.update(completion="pending"), lambda s:s.update(schema_version="future.v2"),
    lambda s:s.update(read_count=True), lambda s:s.update(read_count=2**53), lambda s:s.update(changed_reads=10),
    lambda s:s["outcomes"]["exact"].update(unique=4), lambda s:s["parameters"].update(metric="levenshtein"),
    lambda s:s["parameters"].update(orientation="reverse_complement"), lambda s:s["parameters"].update(target_start=-1),
    lambda s:s["parameters"].update(target_length=0), lambda s:s["parameters"].update(target_start=2**53-1),
    lambda s:s.update(native_library_sha256="not-a-hash"), lambda s:s["count_comparisons"][0].update(counts_identical=True),
    lambda s:s["count_comparisons"][0].update(counts_identical=0), lambda s:s["count_comparisons"][0].update(differing_guides=99),
    lambda s:s["count_comparisons"][0].update(total_count_delta=True), lambda s:s.update(count_comparisons=[]),
    lambda s:s["inputs"]["reads"].update(sha256="x"), lambda s:s["outcomes"].update(other={}),
])
def test_rejects_invalid_summaries(tmp_path, change):
    s = fixture(tmp_path)
    change(s)
    with pytest.raises(ValueError):
        review.build_review_data(s, tmp_path)


@pytest.mark.parametrize("mutation", ["delta", "count", "duplicate", "missing", "header", "extra", "fraction", "negative"])
def test_rejects_contradictory_guide_tables(tmp_path, mutation):
    s = fixture(tmp_path)
    path = tmp_path / "guide_deltas.tsv"
    rows = list(csv.reader(io.StringIO(path.read_text()), delimiter="\t"))
    if mutation == "delta": rows[1][-2] = "99"
    if mutation == "count": rows[1][2] = "99"
    if mutation == "duplicate": rows[2][0] = rows[1][0]
    if mutation == "missing": rows.pop()
    if mutation == "header": rows[0][1] = rows[0][0]
    if mutation == "extra": rows[1].append("unexpected")
    if mutation == "fraction": rows[1][2] = "1.5"
    if mutation == "negative": rows[1][2] = "-1"
    write_tsv(path, rows[0], rows[1:])
    with pytest.raises(ValueError): review.build_review_data(s, tmp_path)


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "count", "policy", "state", "fraction"])
def test_rejects_invalid_transitions(tmp_path, mutation):
    s = fixture(tmp_path)
    path = tmp_path / "transitions.tsv"
    rows = list(csv.reader(io.StringIO(path.read_text()), delimiter="\t"))
    if mutation == "missing": rows.pop()
    if mutation == "duplicate": rows[2] = rows[1]
    if mutation == "count": rows[1][-1] = "88"
    if mutation == "policy": rows[1][0] = "other"
    if mutation == "state": rows[1][2] = "other"
    if mutation == "fraction": rows[1][-1] = "NaN"
    write_tsv(path, rows[0], rows[1:])
    with pytest.raises(ValueError): review.build_review_data(s, tmp_path)


def test_completed_bundle_hash_checks(tmp_path):
    s = seal(tmp_path, fixture(tmp_path))
    assert review.build_review_data(s, tmp_path)["source_manifest_checked"] is True
    path = tmp_path / "guide_deltas.tsv"
    path.write_bytes(path.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="hash"): review.build_review_data(s, tmp_path)


def test_missing_source_hash_is_not_a_verification(tmp_path):
    s = seal(tmp_path, fixture(tmp_path))
    del s["artifacts"]["transitions.tsv"]
    with pytest.raises(ValueError): review.build_review_data(s, tmp_path)


def test_read_attachment_need_not_be_copied_with_portable_rebuild(tmp_path):
    s = seal(tmp_path, fixture(tmp_path))
    (tmp_path / "read_changes.tsv").unlink()
    assert "read_changes.tsv" in review.build_review_data(s, tmp_path)["sources"]


def test_no_read_attachment_when_not_recorded(tmp_path):
    s = fixture(tmp_path)
    (tmp_path / "read_changes.tsv").unlink()
    data = review.build_review_data(s, tmp_path)
    assert "read_changes.tsv" not in data["sources"]


def test_rejects_symlink(tmp_path):
    s = fixture(tmp_path)
    path = tmp_path / "guide_deltas.tsv"
    path.rename(tmp_path / "other.tsv")
    path.symlink_to(tmp_path / "other.tsv")
    with pytest.raises(ValueError, match="non-symlink"): review.build_review_data(s, tmp_path)


def test_capacity_error_is_distinct_from_corruption(tmp_path, monkeypatch):
    s = fixture(tmp_path)
    monkeypatch.setattr(review, "MAX_GUIDES", 4)
    with pytest.raises(review.ReviewCapacityError): review.build_review_data(s, tmp_path)


def test_source_size_is_bounded(tmp_path, monkeypatch):
    s = fixture(tmp_path)
    monkeypatch.setattr(review, "MAX_TABLE_BYTES", 20)
    with pytest.raises(review.ReviewCapacityError): review.build_review_data(s, tmp_path)


def test_cli_no_overwrite_or_input_mutation(tmp_path):
    bundle = tmp_path / "bundle"
    seal(bundle, fixture(bundle))
    before = {p.name: p.read_bytes() for p in bundle.iterdir()}
    destination = tmp_path / "new-review.html"
    assert review.main(["--bundle", str(bundle), "--out", str(destination)]) == 0
    original = destination.read_bytes()
    assert review.main(["--bundle", str(bundle), "--out", str(destination)]) == 2
    assert destination.read_bytes() == original
    assert {p.name: p.read_bytes() for p in bundle.iterdir()} == before


@pytest.mark.parametrize("text", ['{"schema_version":"a","schema_version":"b"}', '{"a":NaN}', '{"a":Infinity}'])
def test_strict_json(text):
    with pytest.raises(ValueError): review._json(text)


def test_renderer_hash_binds_both_modules(tmp_path):
    data = review.build_review_data(fixture(tmp_path), tmp_path)
    assert data["renderer_sha256"] == sha((SOURCE / "sensitivity_review.py").read_bytes()+b"\0"+(SOURCE / "sensitivity_review_assets.py").read_bytes())


def test_untrusted_extra_metadata_is_not_embedded(tmp_path):
    s = fixture(tmp_path)
    s["private_notes"] = "DO_NOT_EXPORT_THIS"
    s["inputs"]["reads"]["private_path"] = "/secret/unpublished"
    document = review.render_sensitivity_report(s, tmp_path)
    assert "DO_NOT_EXPORT_THIS" not in document and "/secret/unpublished" not in document


def native_inputs(directory):
    directory.mkdir(parents=True, exist_ok=True)
    targets = directory / "targets.tsv"
    targets.write_text("target_id\tsequence\ng1\tAAAAAAAA\ng2\tCAAAAAAA\ng3\tTTTTTTTT\n")
    reads = directory / "reads.fastq"
    reads.write_text("".join(f"@private_read_{i}\n{seq}\n+\n{'I'*len(seq)}\n" for i,seq in enumerate(["AAAAAAAA","CAAAAAAA","TTTTTTTT","ATTTTTTT"])))
    return dict(targets=targets, reads=reads, target_start=0, target_length=8, write_read_changes=True)


def test_native_producer_embeds_complete_review_and_preserves_hashes(tmp_path):
    engine = pytest.importorskip("dotmatch.sensitivity", reason="Native engine integration requires a built DotMatch checkout")
    args = native_inputs(tmp_path / "input")
    result = tmp_path / "result"
    summary = engine.run_sensitivity(**args, out_dir=result)
    document = (result / "report.html").read_text()
    assert 'dotmatch.review.v1' in document
    assert 'private_read_' not in document
    assert summary["execution"]["fastq_passes"] == 1
    for name, info in summary["artifacts"].items():
        assert sha((result / name).read_bytes()) == info["sha256"]
    assert json.loads((result / "summary.json").read_text()) == summary


def test_native_capacity_fallback_does_not_discard_counts(tmp_path, monkeypatch):
    engine = pytest.importorskip("dotmatch.sensitivity", reason="Native engine integration requires a built DotMatch checkout")
    import dotmatch.sensitivity_review as renderer
    monkeypatch.setattr(renderer, "MAX_GUIDES", 2)
    summary = engine.run_sensitivity(**native_inputs(tmp_path / "input"), out_dir=tmp_path / "result")
    assert summary["interactive_review"]["available"] is False
    assert "Interactive review unavailable" in (tmp_path / "result/report.html").read_text()
    assert len(list(csv.DictReader((tmp_path / "result/guide_deltas.tsv").open(), delimiter="\t"))) == 3


def test_native_report_corruption_fails_closed(tmp_path, monkeypatch):
    engine = pytest.importorskip("dotmatch.sensitivity", reason="Native engine integration requires a built DotMatch checkout")
    import dotmatch.sensitivity_review as renderer
    def fail(*args, **kwargs):
        raise ValueError("contradictory evidence")
    monkeypatch.setattr(renderer, "build_review_data", fail)
    result = tmp_path / "result"
    with pytest.raises(ValueError, match="contradictory"):
        engine.run_sensitivity(**native_inputs(tmp_path / "input"), out_dir=result)
    assert not result.exists()


def test_staged_hashes_are_not_labelled_published_manifest(tmp_path):
    s = seal(tmp_path, fixture(tmp_path))
    assert review.build_review_data(s, tmp_path, staged=True)["source_manifest_checked"] is False


def test_overlong_identifier_is_viewer_capacity_not_lost_analysis(tmp_path):
    s = fixture(tmp_path, identifiers=["x"*16385,"b","c","d","e"])
    with pytest.raises(review.ReviewCapacityError): review.build_review_data(s, tmp_path)


def test_artifact_size_must_be_integer(tmp_path):
    s = seal(tmp_path, fixture(tmp_path))
    s["artifacts"]["guide_deltas.tsv"]["bytes"] = float(s["artifacts"]["guide_deltas.tsv"]["bytes"])
    with pytest.raises(ValueError): review.build_review_data(s, tmp_path)
