"""Independent arithmetic and I/O contracts for the side-by-side evaluator."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import random
import subprocess
import sys
from pathlib import Path

import pytest

from dotmatch.count_compare import TOP_CHANGES, compare_tables, main, write_comparison
from dotmatch.count_io import read_count_table


def table(tmp_path, name, content):
    path = tmp_path / name
    path.write_text(content, encoding="utf-8")
    return path


def pair(tmp_path, before="10\n", after="11\n"):
    return (table(tmp_path, "before.tsv", "sgRNA\tGene\tsample\ng1\tG1\t" + before),
            table(tmp_path, "after.tsv", "sgRNA\tGene\tsample\ng1\tG1\t" + after))


def test_equal_totals_do_not_hide_redistribution(tmp_path):
    a = table(tmp_path, "a.tsv", "sgRNA\tGene\ts\ng1\tG1\t10\ng2\tG2\t0\n")
    b = table(tmp_path, "b.tsv", "sgRNA\tGene\ts\ng2\tG2\t10\ng1\tG1\t0\n")
    report = write_comparison(a, b, tmp_path / "out")
    s = report["samples"][0]
    assert s["baseline_compared_total"] == s["candidate_compared_total"] == 10
    assert s["same_total_different_counts"] is True
    assert s["changed_guides"] == 2
    assert s["absolute_count_delta"] == 20
    assert s["newly_nonzero_guides"] == s["newly_zero_guides"] == 1
    assert "Totals agree. Guide counts do not." in (tmp_path / "out/report.html").read_text()


def test_row_and_sample_reordering_is_not_a_difference(tmp_path):
    a = table(tmp_path, "a.tsv", "sgRNA\tGene\ta\tb\nNA\tG1\t0\t12\n001\tG2\t13\t0\n")
    b = table(tmp_path, "b.tsv", "sgRNA\tGene\tb\ta\n001\tG2\t0\t13\nNA\tG1\t12\t0\n")
    result = compare_tables(read_count_table(a), read_count_table(b))
    assert result["any_count_difference"] is False
    assert result["alignment"]["guides"] == ["NA", "001"]
    assert [row["baseline_compared_total"] for row in result["samples"]] == [13, 12]


def test_detailed_count_components_are_not_extra_samples(tmp_path):
    a = table(tmp_path, "a.tsv", "sgRNA\tGene\ts\ng1\tG1\t10\n")
    b = table(tmp_path, "b.tsv", "target_id\ttarget_seq\tgene\ts_count_exact\ts_count_corrected\ts_count_total\ng1\tAAAA\tG1\t8\t2\t10\n")
    report = compare_tables(read_count_table(a), read_count_table(b))
    assert len(report["samples"]) == 1
    assert report["any_count_difference"] is False
    assert report["identity"]["candidate_only_fields"] == ["target_seq"]


def test_explicit_sample_map_preserves_auditable_source_names(tmp_path):
    a = table(tmp_path, "a.tsv", "sgRNA\tGene\tpatient_1\tpatient_2\ng1\tG1\t10\t0\n")
    b = table(tmp_path, "b.tsv", "sgRNA\tGene\tL002\tL001\ng1\tG1\t0\t11\n")
    mapping = table(tmp_path, "samples.tsv", "baseline\tcandidate\npatient_1\tL001\npatient_2\tL002\n")
    assert main(["--baseline", str(a), "--candidate", str(b), "--out-dir", str(tmp_path / "out"),
                 "--sample-map", str(mapping), "--fail-on-difference"]) == 1
    report = json.loads((tmp_path / "out/report.json").read_text())
    assert [s["sample"] for s in report["samples"]] == ["patient_1", "patient_2"]
    assert [s["total_delta"] for s in report["samples"]] == [1, 0]
    assert report["sample_map"]["sha256"] == hashlib.sha256(mapping.read_bytes()).hexdigest()
    assert report["sample_map"]["pairs"] == [
        {"baseline": "patient_1", "candidate": "L001"},
        {"baseline": "patient_2", "candidate": "L002"},
    ]
    assert "asserted biological sample identity" in (tmp_path / "out/report.html").read_text()


def test_guide_map_aligns_reordered_ids_and_checks_counts(tmp_path):
    a = table(tmp_path, "a.tsv", "sgRNA\tGene\ts\ng1\tG1\t10\ng2\tG2\t0\n")
    b = table(tmp_path, "b.tsv", "sgRNA\tGene\ts\nnew2\tG2\t1\nnew1\tG1\t9\n")
    mapping = table(tmp_path, "guides.tsv", "baseline\tcandidate\ng1\tnew1\ng2\tnew2\n")
    assert main(["--baseline", str(a), "--candidate", str(b), "--out-dir", str(tmp_path / "out"),
                 "--guide-map", str(mapping), "--fail-on-difference"]) == 1
    report = json.loads((tmp_path / "out/report.json").read_text())
    assert report["alignment"]["guides"] == ["g1", "g2"]
    assert report["samples"][0]["same_total_different_counts"]
    assert report["guide_map"]["sha256"] == hashlib.sha256(mapping.read_bytes()).hexdigest()
    assert report["guide_map"]["pairs"] == [{"baseline": "g1", "candidate": "new1"},
                                             {"baseline": "g2", "candidate": "new2"}]
    with (tmp_path / "out/changes.tsv").open(newline="") as handle:
        assert [row["guide"] for row in csv.DictReader(handle, delimiter="\t")] == ["g1", "g2"]


@pytest.mark.parametrize("mapping", [
    "baseline\tcandidate\ng1\tnew1\ng1\tnew2\n",
    "baseline\tcandidate\ng1\tnew1\ng2\tnew1\n",
    "baseline\tcandidate\ng1\tmissing\n",
    "baseline\tcandidate\ng1\tnew1\textra\n",
    "baseline\tcandidate\n",
])
def test_invalid_guide_map_refused_before_output(tmp_path, mapping):
    a = table(tmp_path, "a.tsv", "sgRNA\tGene\ts\ng1\tG1\t1\ng2\tG2\t2\n")
    b = table(tmp_path, "b.tsv", "sgRNA\tGene\ts\nnew1\tG1\t1\nnew2\tG2\t2\n")
    m = table(tmp_path, "guides.tsv", mapping)
    with pytest.raises(ValueError):
        write_comparison(a, b, tmp_path / "out", guide_map=m)
    assert not (tmp_path / "out").exists()


def test_guide_map_collision_and_annotation_conflict_refused(tmp_path):
    a = table(tmp_path, "a.tsv", "sgRNA\tGene\ts\ng1\tG1\t1\ng2\tG2\t2\n")
    b = table(tmp_path, "b.tsv", "sgRNA\tGene\ts\nnew1\tG1\t1\ng2\tG2\t2\n")
    m = table(tmp_path, "guides.tsv", "baseline\tcandidate\ng2\tnew1\n")
    with pytest.raises(ValueError, match="collides"):
        write_comparison(a, b, tmp_path / "out", guide_map=m)
    m.write_text("baseline\tcandidate\ng1\tnew1\n")
    b.write_text("sgRNA\tGene\ts\nnew1\tWRONG\t1\ng2\tG2\t2\n")
    with pytest.raises(ValueError, match="conflicting gene"):
        write_comparison(a, b, tmp_path / "out", guide_map=m)
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("mapping", [
    "baseline\tcandidate\npatient_1\tL001\npatient_1\tL002\n",
    "baseline\tcandidate\npatient_1\tL001\npatient_2\tL001\n",
    "baseline\tcandidate\npatient_1\tmissing\n",
    "baseline\tcandidate\npatient_1\tL001\textra\n",
    "baseline\tcandidate\n",
])
def test_invalid_sample_map_refused_before_output(tmp_path, mapping):
    a = table(tmp_path, "a.tsv", "sgRNA\tGene\tpatient_1\tpatient_2\ng1\tG1\t10\t0\n")
    b = table(tmp_path, "b.tsv", "sgRNA\tGene\tL001\tL002\ng1\tG1\t10\t0\n")
    m = table(tmp_path, "samples.tsv", mapping)
    with pytest.raises(ValueError):
        write_comparison(a, b, tmp_path / "out", sample_map=m)
    assert not (tmp_path / "out").exists()


def test_partial_map_collision_refused(tmp_path):
    a = table(tmp_path, "a.tsv", "sgRNA\tGene\ta\tb\ng1\tG1\t1\t2\n")
    b = table(tmp_path, "b.tsv", "sgRNA\tGene\ta\tx\ng1\tG1\t1\t2\n")
    m = table(tmp_path, "samples.tsv", "baseline\tcandidate\na\tx\n")
    with pytest.raises(ValueError, match="collides"):
        write_comparison(a, b, tmp_path / "out", sample_map=m)
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("bad", ["", "NaN", "inf", "-1", "0.5", "1e1000", "bad"])
def test_invalid_counts_refused_before_any_output(tmp_path, bad):
    a, b = pair(tmp_path, after=bad + "\n")
    with pytest.raises(ValueError):
        write_comparison(a, b, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_large_integer_difference_never_rounded_away(tmp_path):
    huge = 10 ** 999 + 2 ** 53
    a, b = pair(tmp_path, str(huge) + "\n", str(huge + 1) + "\n")
    report = write_comparison(a, b, tmp_path / "out")
    assert report["samples"][0]["total_delta"] == 1
    parsed = json.loads((tmp_path / "out/report.json").read_text())
    assert parsed["samples"][0]["candidate_compared_total"] == huge + 1


def test_integral_decimal_notation_is_accepted(tmp_path):
    a, b = pair(tmp_path, "10.0\n", "1e1\n")
    assert not compare_tables(read_count_table(a), read_count_table(b))["any_count_difference"]


def test_zero_tables_have_no_nan_or_accuracy_claim(tmp_path):
    a, b = pair(tmp_path, "0\n", "0\n")
    result = write_comparison(a, b, tmp_path / "out")
    assert not result["any_count_difference"]
    assert not result["samples"][0]["same_total_different_counts"]
    assert "Agreement is not accuracy" in result["interpretation"]


def test_mismatched_guide_sets_refused_by_default(tmp_path):
    a, b = pair(tmp_path)
    b.write_text("sgRNA\tGene\tsample\ng2\tG2\t1\n")
    with pytest.raises(ValueError, match="guide/sample sets differ"):
        write_comparison(a, b, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_shared_only_reports_excluded_axes_and_count_mass(tmp_path):
    a = table(tmp_path, "a.tsv", "sgRNA\tGene\ts\textra_a\ng1\tG1\t10\t1\ng2\tG2\t7\t9\n")
    b = table(tmp_path, "b.tsv", "sgRNA\tGene\ts\textra_b\ng1\tG1\t10\t2\ng3\tG3\t11\t9\n")
    report = write_comparison(a, b, tmp_path / "out", shared_only=True)
    assert report["scope"] == "shared_only_partial"
    assert report["alignment"]["excluded"] == {
        "baseline_only_guides": ["g2"], "candidate_only_guides": ["g3"],
        "baseline_only_samples": ["extra_a"], "candidate_only_samples": ["extra_b"],
    }
    s = report["samples"][0]
    assert s["baseline_compared_total"] == s["candidate_compared_total"] == 10
    assert s["baseline_input_total"] == 17
    assert s["candidate_input_total"] == 21
    assert s["baseline_excluded_total"] == 7
    assert s["candidate_excluded_total"] == 11
    assert "Partial comparison" in (tmp_path / "out/report.html").read_text()
    assert main(["--baseline", str(a), "--candidate", str(b), "--out-dir", str(tmp_path / "ci"),
                 "--shared-only", "--fail-on-difference"]) == 1


@pytest.mark.parametrize("content", ["sgRNA\tGene\tsample\ng2\tG2\t2\n", "sgRNA\tGene\tother\ng1\tG1\t2\n"])
def test_shared_only_refuses_no_overlap(tmp_path, content):
    a, b = pair(tmp_path)
    b.write_text(content)
    with pytest.raises(ValueError, match="at least one shared"):
        write_comparison(a, b, tmp_path / "out", shared_only=True)


@pytest.mark.parametrize("field", ["Gene", "gene_id", "target_seq"])
def test_identity_conflicts_are_not_silently_compared(tmp_path, field):
    a = table(tmp_path, "a.tsv", f"sgRNA\t{field}\ts\ng1\tAAAA\t1\n")
    b = table(tmp_path, "b.tsv", f"sgRNA\t{field.lower()}\ts\ng1\tCCCC\t1\n")
    with pytest.raises(ValueError, match="conflicting"):
        write_comparison(a, b, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_no_metadata_is_explicitly_unverified(tmp_path):
    a = table(tmp_path, "a.tsv", "guide_id\ts\ng1\t1\n")
    report = write_comparison(a, a, tmp_path / "out")
    assert report["identity"]["checked_fields"] == []
    assert "none available in both files" in (tmp_path / "out/report.html").read_text()


@pytest.mark.parametrize("bad", ["sgRNA\tGene\ts\ng1\tG1\t1\ng1\tG1\t2\n",
                                  "sgRNA\tGene\ts\ts\ng1\tG1\t1\t2\n",
                                  "sgRNA\tGene\ts\n", "sgRNA\tGene\ts\ng1\tG1\n"])
def test_duplicate_or_malformed_axes_rejected(tmp_path, bad):
    a, b = pair(tmp_path)
    b.write_text(bad)
    with pytest.raises(ValueError):
        write_comparison(a, b, tmp_path / "out")


def test_gzip_bom_and_hashes_are_from_exact_original_bytes(tmp_path):
    a, b = pair(tmp_path)
    packed = tmp_path / "a.tsv.gz"
    raw = gzip.compress(b"\xef\xbb\xbf" + a.read_bytes(), mtime=0)
    packed.write_bytes(raw)
    report = write_comparison(packed, b, tmp_path / "out")
    assert report["inputs"]["baseline"] == {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}
    assert str(tmp_path) not in (tmp_path / "out/report.json").read_text()
    manifest = json.loads((tmp_path / "out/manifest.json").read_text())
    assert manifest["status"] == "complete"
    for name, info in manifest["files"].items():
        data = (tmp_path / "out" / name).read_bytes()
        assert info == {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}


def test_truncated_gzip_is_readable_cli_error(tmp_path, capsys):
    a, b = pair(tmp_path)
    packed = tmp_path / "broken.gz"
    packed.write_bytes(gzip.compress(a.read_bytes())[:-5])
    with pytest.raises(SystemExit) as exc:
        main(["--baseline", str(packed), "--candidate", str(b), "--out-dir", str(tmp_path / "out")])
    assert exc.value.code == 2
    assert "Traceback" not in capsys.readouterr().err


def test_existing_output_directory_and_inputs_preserved(tmp_path):
    a, b = pair(tmp_path)
    output = tmp_path / "out"
    output.mkdir()
    sentinel = output / "keep.txt"
    sentinel.write_text("do not touch")
    before = (a.read_bytes(), b.read_bytes())
    with pytest.raises(FileExistsError):
        write_comparison(a, b, output)
    assert sentinel.read_text() == "do not touch"
    assert sorted(p.name for p in output.iterdir()) == ["keep.txt"]
    assert (a.read_bytes(), b.read_bytes()) == before
    with pytest.raises(FileExistsError):
        write_comparison(a, b, a)
    assert a.read_bytes() == before[0]


def test_failed_write_never_publishes_completion_manifest(tmp_path, monkeypatch):
    a, b = pair(tmp_path)
    original = Path.open
    def failing_open(self, *args, **kwargs):
        if self.name == "report.html":
            raise OSError("simulated disk failure")
        return original(self, *args, **kwargs)
    monkeypatch.setattr(Path, "open", failing_open)
    with pytest.raises(OSError):
        write_comparison(a, b, tmp_path / "out")
    assert not (tmp_path / "out/manifest.json").exists()


def test_user_identifiers_escaped_in_html_and_preserved_in_tsv(tmp_path):
    label = '<script>alert("x")</script>'
    a = table(tmp_path, "a.tsv", f"sgRNA\tGene\t{label}\n{label}\tG1\t1\n")
    b = table(tmp_path, "b.tsv", f"sgRNA\tGene\t{label}\n{label}\tG1\t2\n")
    write_comparison(a, b, tmp_path / "out")
    markup = (tmp_path / "out/report.html").read_text()
    assert "<script>" not in markup
    assert "&lt;script&gt;" in markup
    assert "default-src 'none'" in markup
    with (tmp_path / "out/changes.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["guide"] == rows[0]["sample"] == label


def test_html_previews_bounded_but_tsv_complete(tmp_path):
    n = TOP_CHANGES + 13
    a = table(tmp_path, "a.tsv", "sgRNA\tGene\ts\n" + "".join(f"g{i}\tG{i}\t0\n" for i in range(n)))
    b = table(tmp_path, "b.tsv", "sgRNA\tGene\ts\n" + "".join(f"g{i}\tG{i}\t{i+1}\n" for i in range(n)))
    report = write_comparison(a, b, tmp_path / "out")
    assert len(report["samples"][0]["top_changes"]) == TOP_CHANGES
    assert report["samples"][0]["top_changes"][0]["delta"] == n
    with (tmp_path / "out/changes.tsv").open(newline="") as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == n


def test_repeated_outputs_are_byte_identical(tmp_path):
    a, b = pair(tmp_path)
    write_comparison(a, b, tmp_path / "first")
    write_comparison(a, b, tmp_path / "second")
    for first in (tmp_path / "first").iterdir():
        assert first.read_bytes() == (tmp_path / "second" / first.name).read_bytes()


@pytest.mark.parametrize("after,expected", [("10\n", 0), ("11\n", 1)])
def test_cli_exit_codes_and_reports(tmp_path, after, expected):
    a, b = pair(tmp_path, after=after)
    assert main(["--baseline", str(a), "--candidate", str(b), "--out-dir", str(tmp_path / "out"),
                 "--fail-on-difference"]) == expected
    assert (tmp_path / "out/manifest.json").exists()


def test_module_entrypoint_help():
    result = subprocess.run([sys.executable, "-m", "dotmatch.count_compare", "--help"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "--shared-only" in result.stdout


def test_seeded_independent_arithmetic_oracle(tmp_path):
    rng = random.Random(25092026)
    for iteration in range(100):
        ids, samples = [f"g{i}" for i in range(17)], ["a", "b", "c"]
        left = {(g, s): rng.randrange(0, 1000) for g in ids for s in samples}
        right = {(g, s): (left[g, s] if rng.random() < .3 else rng.randrange(0, 1000)) for g in ids for s in samples}
        def write(name, values, guides, names):
            return table(tmp_path, name, "sgRNA\tGene\t" + "\t".join(names) + "\n" + "".join(
                g + "\tG" + g + "\t" + "\t".join(str(values[g, s]) for s in names) + "\n" for g in guides))
        shuffled = ids[:]
        rng.shuffle(shuffled)
        a = read_count_table(write("a.tsv", left, ids, samples))
        b = read_count_table(write("b.tsv", right, shuffled, samples[::-1]))
        result = compare_tables(a, b)
        for summary in result["samples"]:
            s = summary["sample"]
            deltas = [right[g, s] - left[g, s] for g in ids]
            assert summary["total_delta"] == sum(deltas), iteration
            assert summary["absolute_count_delta"] == sum(abs(d) for d in deltas)
            assert summary["changed_guides"] == sum(d != 0 for d in deltas)
            assert summary["same_total_different_counts"] == (sum(deltas) == 0 and any(deltas))
