import csv
import gzip
import hashlib
import itertools
import json
import random
from pathlib import Path

import pytest

from dotmatch.core import (
    Matcher,
    MATCH_UNIQUE,
    MATCH_NONE,
    MATCH_AMBIGUOUS,
    MATCH_INVALID,
    load_targets,
)
from dotmatch.sensitivity import (
    compare_hamming_policies,
    run_sensitivity,
    MODES,
    STATES,
    main,
)
from dotmatch.target_io import read_target_table


def oracle(window, targets, k, policy):
    if window is None:
        return MATCH_INVALID, -1
    distances = [
        (sum(a != b for a, b in zip(window, seq)), i)
        for i, seq in enumerate(targets)
        if len(seq) == len(window)
    ]
    candidates = [(distance, i) for distance, i in distances if distance <= k]
    if not candidates:
        return MATCH_NONE, -1
    if policy == "best":
        best = min(distance for distance, _ in candidates)
        candidates = [(d, i) for d, i in candidates if d == best]
    return (
        (MATCH_UNIQUE, candidates[0][1])
        if len(candidates) == 1
        else (MATCH_AMBIGUOUS, -1)
    )


def write_fastq(path, sequences):
    text = "".join(
        f'@record_{i}\n{s}\n+\n{"I" * len(s)}\n' for i, s in enumerate(sequences)
    )
    if str(path).endswith(".gz"):
        path.write_bytes(gzip.compress(text.encode(), mtime=0))
    else:
        path.write_text(text)


def inputs(tmp_path):
    library = tmp_path / "library.csv.gz"
    library.write_bytes(
        gzip.compress(
            b"gene,sequence,id\nG1,AAAAAAAA,a\nG2,CAAAAAAA,b\nG3,TTTTTTTT,c\nG4,CCCCCCCC,d\nG5,CCCCCCCC,e\n",
            mtime=0,
        )
    )
    reads = tmp_path / "sample.fastq.gz"
    write_fastq(
        reads,
        [
            "AAAAAAAA",
            "CAAAAAAA",
            "GAAAAAAA",
            "TTTTTTTT",
            "ATTTTTTT",
            "GGGGGGGG",
            "AA",
            "CCCCCCCC",
            "NTTTTTTT",
        ],
    )
    return library, reads


def test_fused_policies_exhaustive_and_literal_symbols():
    targets = ["AAAA", "CAAA", "ACAA", "CCCC", "CCCC", "TTTT", "NNNN"]
    windows = ["".join(chars) for chars in itertools.product("ACGTN", repeat=4)] + [
        None
    ]
    with Matcher(targets) as matcher:
        results = compare_hamming_policies(matcher, windows)
    for window, result in zip(windows, results):
        for mode, k, policy in [
            ("exact", 0, "radius"),
            ("radius_k1", 1, "radius"),
            ("best_k1", 1, "best"),
        ]:
            call = getattr(result, mode)
            assert (call.status, call.target_index) == oracle(
                window, targets, k, policy
            )


def test_seeded_longer_libraries_and_duplicates():
    rng = random.Random(60109)
    for length in (8, 20, 32, 64):
        targets = ["".join(rng.choices("ACGT", k=length)) for _ in range(40)]
        targets += [targets[0], "A" + targets[1][1:]]
        windows = targets + ["N" + s[1:] for s in targets] + [None]
        with Matcher(targets) as matcher:
            results = compare_hamming_policies(matcher, windows)
        for window, result in zip(windows, results):
            for mode, k, policy in [
                ("exact", 0, "radius"),
                ("radius_k1", 1, "radius"),
                ("best_k1", 1, "best"),
            ]:
                call = getattr(result, mode)
                assert (call.status, call.target_index) == oracle(
                    window, targets, k, policy
                )


def test_one_native_query_for_duplicate_windows(monkeypatch):
    with Matcher(["AAAAAAAA", "CAAAAAAA"]) as matcher:
        original = matcher.assign_hamming
        seen = []

        def spy(windows, **kwargs):
            seen.append(list(windows))
            return original(windows, **kwargs)

        monkeypatch.setattr(matcher, "assign_hamming", spy)
        results = compare_hamming_policies(
            matcher, ["AAAAAAAA"] * 100 + [None, "GGGGGGGG"]
        )
        assert seen == [["AAAAAAAA", "GGGGGGGG"]]
        assert len(results) == 102
        assert results[0].exact.status == MATCH_UNIQUE
        assert results[0].radius_k1.status == MATCH_AMBIGUOUS
        assert results[0].best_k1.status == MATCH_UNIQUE


def test_complete_bundle_conserves_reads_checksums_and_changes(tmp_path):
    library, reads = inputs(tmp_path)
    output = tmp_path / "result"
    summary = run_sensitivity(
        targets=library,
        reads=reads,
        target_start=0,
        target_length=8,
        out_dir=output,
        batch_size=2,
        write_read_changes=True,
    )
    assert summary["read_count"] == 9
    assert summary["changed_reads"] == 5
    assert (
        summary["inputs"]["reads"]["sha256"]
        == hashlib.sha256(reads.read_bytes()).hexdigest()
    )
    assert (
        summary["inputs"]["targets"]["sha256"]
        == hashlib.sha256(library.read_bytes()).hexdigest()
    )
    assert summary["execution"]["fastq_passes"] == 1
    assert summary["outcomes"]["exact"] == {
        "unique": 3,
        "ambiguous": 1,
        "none": 4,
        "invalid": 1,
    }
    assert summary["outcomes"]["radius_k1"] == {
        "unique": 3,
        "ambiguous": 4,
        "none": 1,
        "invalid": 1,
    }
    assert summary["outcomes"]["best_k1"] == {
        "unique": 5,
        "ambiguous": 2,
        "none": 1,
        "invalid": 1,
    }
    assert summary["count_comparisons"][0]["total_count_delta"] == 0
    assert not summary["count_comparisons"][0]["counts_identical"]
    for mode in MODES:
        assert sum(summary["outcomes"][mode].values()) == 9
        rows = list(csv.reader((output / f"{mode}.counts.tsv").open(), delimiter="\t"))
        assert len(rows) == 6
        assert (
            sum(int(row[2]) for row in rows[1:]) == summary["outcomes"][mode]["unique"]
        )
    for name, metadata in summary["artifacts"].items():
        assert (
            metadata["sha256"]
            == hashlib.sha256((output / name).read_bytes()).hexdigest()
        )
    changes = list(csv.DictReader((output / "read_changes.tsv").open(), delimiter="\t"))
    assert len(changes) == 5
    assert all("sequence" not in key for key in changes[0])
    assert json.loads((output / "summary.json").read_text()) == summary
    assert not any(p.name.startswith(".pending") for p in output.iterdir())


def test_plain_and_gzip_results_are_equivalent_across_batch_sizes(tmp_path):
    library, reads = inputs(tmp_path)
    plain = tmp_path / "sample.fastq"
    plain.write_bytes(gzip.decompress(reads.read_bytes()))
    args = dict(targets=library, target_start=0, target_length=8)
    a = run_sensitivity(**args, reads=reads, out_dir=tmp_path / "a", batch_size=1)
    b = run_sensitivity(**args, reads=plain, out_dir=tmp_path / "b", batch_size=100)
    for key in ["outcomes", "count_comparisons", "changed_reads"]:
        assert a[key] == b[key]
    for name in ["guide_deltas.tsv", "transitions.tsv", "exact.counts.tsv"]:
        assert (tmp_path / "a" / name).read_bytes() == (
            tmp_path / "b" / name
        ).read_bytes()


@pytest.mark.parametrize(
    "content",
    [
        "",
        "@a\nAAAA\n+\n",
        "@\nAAAA\n+\nIIII\n",
        "@a\nAAAA\n+\nIII\n",
        "@a\nAAAA\n+\nI I!\n",
    ],
)
def test_invalid_fastq_leaves_no_completed_output(tmp_path, content):
    library, _ = inputs(tmp_path)
    reads = tmp_path / "bad.fastq"
    reads.write_text(content)
    with pytest.raises(ValueError):
        run_sensitivity(
            targets=library,
            reads=reads,
            target_start=0,
            target_length=8,
            out_dir=tmp_path / "result",
            write_read_changes=True,
        )
    assert not (tmp_path / "result").exists()


def test_corrupt_gzip_fails_closed(tmp_path):
    library, reads = inputs(tmp_path)
    reads.write_bytes(reads.read_bytes()[:-6])
    with pytest.raises((EOFError, OSError)):
        run_sensitivity(
            targets=library,
            reads=reads,
            target_start=0,
            target_length=8,
            out_dir=tmp_path / "result",
        )
    assert not (tmp_path / "result").exists()


def test_existing_directory_is_not_touched(tmp_path):
    library, reads = inputs(tmp_path)
    output = tmp_path / "result"
    output.mkdir()
    (output / "keep").write_text("important")
    with pytest.raises(FileExistsError):
        run_sensitivity(
            targets=library,
            reads=reads,
            target_start=0,
            target_length=8,
            out_dir=output,
        )
    assert (output / "keep").read_text() == "important"


@pytest.mark.parametrize(
    "override",
    [
        {"target_start": -1},
        {"target_length": 7},
        {"sample_label": "Gene"},
        {"sample_label": "bad\tname"},
        {"batch_size": 0},
        {"batch_size": 65537},
    ],
)
def test_invalid_configuration_fails_before_output(tmp_path, override):
    library, reads = inputs(tmp_path)
    args = dict(
        targets=library,
        reads=reads,
        target_start=0,
        target_length=8,
        out_dir=tmp_path / "result",
    )
    args.update(override)
    with pytest.raises(ValueError):
        run_sensitivity(**args)
    assert not (tmp_path / "result").exists()


def test_report_escapes_identifiers(tmp_path):
    library = tmp_path / "targets.tsv"
    library.write_text(
        "target_id\tsequence\n<script>alert(1)</script>\tAAAAAAAA\nb\tCAAAAAAA\n"
    )
    reads = tmp_path / "reads.fastq"
    write_fastq(reads, ["AAAAAAAA"])
    run_sensitivity(
        targets=library,
        reads=reads,
        target_start=0,
        target_length=8,
        out_dir=tmp_path / "report",
    )
    report = (tmp_path / "report" / "report.html").read_text()
    assert "<script>" not in report and "&lt;script&gt;" in report


def test_cli_routes_without_native_delegation(tmp_path, capsys):
    from dotmatch.cli import main as cli_main

    library, reads = inputs(tmp_path)
    assert (
        cli_main(
            [
                "sensitivity",
                "--targets",
                str(library),
                "--reads",
                str(reads),
                "--target-start",
                "0",
                "--target-length",
                "8",
                "--out-dir",
                str(tmp_path / "out"),
            ]
        )
        == 0
    )
    assert json.loads(capsys.readouterr().out)["completion"] == "complete"
