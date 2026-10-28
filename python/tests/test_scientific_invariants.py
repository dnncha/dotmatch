"""Scientific invariants: oracle equivalence, count conservation, determinism."""

from __future__ import annotations

import csv
import gzip
import tempfile
from pathlib import Path

import dotmatch as dm
import pytest
from dotmatch.count_io import read_count_table
from dotmatch.feature_matrix import build_feature_matrix
from dotmatch.panel import (
    AssignmentConfig,
    DesignOptions,
    PanelConstraints,
    PanelRecord,
    check_panel,
    design_panel,
    hamming_distance,
)
from dotmatch.sensitivity import run_sensitivity
from hypothesis import given, settings
from hypothesis import strategies as st

pytestmark = pytest.mark.scientific

ALPHABET = "ACGT"
DNA = st.text(alphabet=ALPHABET, min_size=1, max_size=8)
TARGETS = st.lists(DNA, min_size=1, max_size=6)
READS = st.lists(DNA, min_size=1, max_size=6)


def _levenshtein(left: str, right: str) -> int:
    previous = list(range(len(right) + 1))
    for i, a in enumerate(left, start=1):
        current = [i]
        for j, b in enumerate(right, start=1):
            current.append(
                min(current[-1] + 1, previous[j] + 1, previous[j - 1] + (a != b))
            )
        previous = current
    return previous[-1]


def _distance(read: str, target: str, metric: str) -> int | None:
    if metric == "exact":
        return 0 if read == target else None
    if metric == "hamming":
        return sum(a != b for a, b in zip(read, target)) if len(read) == len(target) else None
    return _levenshtein(read, target)


def _reference(
    read: str, targets: list[str], metric: str, k: int, policy: str
) -> tuple[dm.MatchStatus, int, int, int, int]:
    candidates = [
        (index, distance)
        for index, target in enumerate(targets)
        if (distance := _distance(read, target, metric)) is not None and distance <= k
    ]
    if not candidates:
        return dm.MATCH_NONE, -1, -1, -1, 0
    best = min(distance for _index, distance in candidates)
    best_matches = [index for index, distance in candidates if distance == best]
    second = min((distance for _index, distance in candidates if distance > best), default=-1)
    if policy == "best":
        status = dm.MATCH_UNIQUE if len(best_matches) == 1 else dm.MATCH_AMBIGUOUS
    else:
        status = dm.MATCH_AMBIGUOUS if len(candidates) > 1 else dm.MATCH_UNIQUE
    return status, best_matches[0], best, second, len(candidates)


@settings(max_examples=250, deadline=None)
@given(TARGETS, READS, st.sampled_from(("levenshtein", "hamming", "exact")), st.integers(0, 2), st.sampled_from(("radius", "best")))
def test_native_assignment_matches_scalar_oracle(
    targets: list[str], reads: list[str], metric: dm.Metric, k: int, policy: dm.AssignmentPolicy
) -> None:
    if metric == "exact":
        k = 0
        results = dm.assign_exact(reads, targets, policy=policy)
    elif metric == "hamming":
        results = dm.assign_hamming(reads, targets, k=k, policy=policy)
    else:
        results = dm.assign(reads, targets, k=k, policy=policy)

    assert [
        (
            result.status,
            result.target_index,
            result.best_distance,
            result.second_best_distance,
            result.match_count,
        )
        for result in results
    ] == [_reference(read, targets, metric, k, policy) for read in reads]


def _fixed_length_case() -> st.SearchStrategy[tuple[list[str], list[str]]]:
    return st.integers(1, 8).flatmap(
        lambda length: st.tuples(
            st.lists(st.text(alphabet=ALPHABET, min_size=length, max_size=length), min_size=1, max_size=6),
            st.lists(st.text(alphabet=ALPHABET, min_size=length, max_size=length), min_size=1, max_size=6),
        )
    )


@settings(max_examples=60, deadline=None)
@given(_fixed_length_case(), st.sampled_from(("levenshtein", "hamming", "exact")))
def test_stream_summary_conserves_every_read(
    case: tuple[list[str], list[str]], metric: dm.Metric
) -> None:
    targets, reads = case
    k = 0 if metric == "exact" else 1
    results = (
        dm.assign_exact(reads, targets)
        if metric == "exact"
        else dm.assign_hamming(reads, targets, k=k)
        if metric == "hamming"
        else dm.assign(reads, targets, k=k)
    )
    import tempfile

    with tempfile.TemporaryDirectory(prefix="dotmatch-science-") as temporary:
        fastq = Path(temporary) / "reads.fastq"
        fastq.write_text(
            "".join(f"@r{i}\n{read}\n+\n{'I' * len(read)}\n" for i, read in enumerate(reads)),
            encoding="ascii",
        )
        streamed = list(dm.stream_assign(fastq, targets, metric=metric, k=k))
    summary = dm.assignment_summary(streamed)

    assert summary["total_reads"] == len(reads)
    assert sum(summary[key] for key in ("assigned_unique", "ambiguous", "unmatched", "invalid")) == len(reads)
    assert summary["assigned_unique"] == sum(1 for result in results if result.status == dm.MATCH_UNIQUE)
    assert summary["ambiguous"] == sum(1 for result in results if result.status == dm.MATCH_AMBIGUOUS)
    assert summary["unmatched"] == sum(1 for result in results if result.status == dm.MATCH_NONE)
    assert summary["invalid"] == sum(1 for result in results if result.status == dm.MATCH_INVALID)
    for key in ("assignment_rate", "ambiguous_rate", "no_match_rate", "invalid_rate"):
        assert 0.0 <= summary[key] <= 1.0


@settings(max_examples=50, deadline=None)
@given(
    st.lists(st.integers(min_value=0, max_value=10**9), min_size=1, max_size=6),
    st.lists(st.integers(min_value=0, max_value=10**9), min_size=1, max_size=6),
)
def test_count_table_round_trip_is_lossless(
    first: list[int], second: list[int]
) -> None:
    rows = [f"g{i}\tGENE{i}\t{a}\t{b}" for i, (a, b) in enumerate(zip(first, second))]
    text = "sgRNA\tGene\tS0\tS1\n" + "\n".join(rows[: min(len(first), len(second))]) + "\n"
    import tempfile

    with tempfile.TemporaryDirectory(prefix="dotmatch-counts-") as temporary:
        path = Path(temporary) / "counts.tsv"
        path.write_text(text, encoding="utf-8")
        table = read_count_table(path)
    expected = list(zip(first, second))[: len(table.target_ids)]
    assert table.counts == tuple((a, b) for a, b in expected)
    assert table.sample_names == ("S0", "S1")
    assert table.target_ids == tuple(f"g{i}" for i in range(len(expected)))


def test_count_table_preserves_integers_above_float_precision(tmp_path: Path) -> None:
    path = tmp_path / "counts.tsv"
    value = 2**53 + 1
    path.write_text(f"sgRNA\tGene\tS0\ng\tG\t{value}\n", encoding="utf-8")

    table = read_count_table(path)

    assert table.counts == ((value,),)


def test_gzip_and_plain_fastq_stream_identically(tmp_path: Path) -> None:
    text = "@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nACGTACGA\n+\nIIIIIIII\n"
    plain = tmp_path / "reads.fastq"
    compressed = tmp_path / "reads.fastq.gz"
    plain.write_text(text, encoding="ascii")
    with gzip.open(compressed, "wt", encoding="ascii") as handle:
        handle.write(text)
    targets = ["ACGTACGT", "TTTTTTTT"]

    assert list(dm.stream_assign(plain, targets, metric="hamming", k=1)) == list(
        dm.stream_assign(compressed, targets, metric="hamming", k=1)
    )


def test_panel_design_is_deterministic_for_a_seed() -> None:
    options = DesignOptions(
        count=4,
        length=8,
        seed=17,
        engine="greedy",
        restarts=3,
        candidate_pool_size=1200,
    )
    constraints = PanelConstraints(
        length=8,
        min_hamming_distance=4,
        min_levenshtein_distance=3,
        min_sequence_levenshtein_distance=3,
    )
    assignment = AssignmentConfig(metric="hamming", k=1)

    first_records, first_trace = design_panel(options, assignment, constraints)
    second_records, second_trace = design_panel(options, assignment, constraints)

    assert first_records == second_records
    assert first_trace == second_trace
    assert len(first_records) == 4


@settings(max_examples=80, deadline=None)
@given(TARGETS, READS, st.sampled_from(("levenshtein", "hamming", "exact")), st.integers(0, 2), st.sampled_from(("radius", "best")))
def test_indexed_and_scan_paths_match(
    targets: list[str], reads: list[str], metric: dm.Metric, k: int, policy: dm.AssignmentPolicy
) -> None:
    if metric == "exact":
        k = 0
        scan = dm.assign_exact(reads, targets, policy=policy)
        with dm.Matcher(targets) as matcher:
            indexed = matcher.assign_exact(reads, policy=policy)
    elif metric == "hamming":
        scan = dm.assign_hamming(reads, targets, k=k, policy=policy)
        with dm.Matcher(targets) as matcher:
            indexed = matcher.assign_hamming(reads, k=k, policy=policy)
    else:
        scan = dm.assign(reads, targets, k=k, policy=policy)
        with dm.Matcher(targets) as matcher:
            indexed = matcher.assign(reads, k=k, policy=policy)

    assert scan == indexed


def _panel_records() -> st.SearchStrategy[list[str]]:
    return st.integers(min_value=2, max_value=4).flatmap(
        lambda length: st.lists(
            st.text(alphabet=ALPHABET, min_size=length, max_size=length),
            min_size=2,
            max_size=4,
        )
    )


@settings(max_examples=40, deadline=None)
@given(_panel_records())
def test_panel_k1_safety_matches_pairwise_distance(sequences: list[str]) -> None:
    records = [PanelRecord(f"BC{i}", sequence) for i, sequence in enumerate(sequences)]
    constraints = PanelConstraints(
        length=len(sequences[0]),
        alphabet=ALPHABET,
        min_hamming_distance=0,
        min_levenshtein_distance=0,
        min_sequence_levenshtein_distance=0,
        gc_min=0.0,
        gc_max=1.0,
        max_homopolymer=len(sequences[0]),
        max_dinucleotide_repeat=len(sequences[0]),
        avoid_self_complementarity=False,
        avoid_ambiguous_bases=False,
        avoid_reverse_complements=False,
        avoid_prefix_collisions=False,
        avoid_suffix_collisions=False,
        forbidden_motifs=[],
        cycle_balance_enabled=False,
    )
    with tempfile.TemporaryDirectory(prefix="dotmatch-panel-science-") as temporary:
        certificate = check_panel(
            records,
            AssignmentConfig(metric="hamming", k=1),
            constraints,
            Path(temporary),
            source_panel=Path("barcodes.tsv"),
        )

    minimum = min(
        hamming_distance(left, right)
        for i, left in enumerate(sequences)
        for right in sequences[i + 1 :]
    )
    assert certificate["minimum_hamming_distance"] == minimum
    assert certificate["safe_for_k1_hamming"] is (minimum >= 3)


@settings(max_examples=25, deadline=None)
@given(_fixed_length_case())
def test_feature_matrix_counts_only_unique_assignments(
    case: tuple[list[str], list[str]]
) -> None:
    targets, reads = case
    target_rows = [f"t{i}\t{sequence}" for i, sequence in enumerate(targets)]
    observation_rows = [f"o{i}\tcell_{i}\t{read}" for i, read in enumerate(reads)]
    with tempfile.TemporaryDirectory(prefix="dotmatch-feature-science-") as temporary:
        root = Path(temporary)
        target_path = root / "targets.tsv"
        observation_path = root / "observations.tsv"
        target_path.write_text("target_id\ttarget_seq\n" + "\n".join(target_rows) + "\n", encoding="utf-8")
        observation_path.write_text(
            "observation_id\tcell_barcode\tfeature_seq\n"
            + "\n".join(observation_rows)
            + "\n",
            encoding="utf-8",
        )
        result = build_feature_matrix(
            observation_path,
            target_path,
            root / "matrix",
            cell_column="cell_barcode",
            sequence_column="feature_seq",
            id_column="observation_id",
            k=1,
            metric="hamming",
        )
        rows = list(csv.DictReader((result.output_dir / "cell_feature_counts.tsv").open(encoding="utf-8"), delimiter="\t"))

    results = dm.assign_hamming(reads, targets, k=1)
    expected: dict[tuple[str, str], int] = {}
    for index, assignment in enumerate(results):
        if assignment.status == dm.MATCH_UNIQUE:
            key = (f"cell_{index}", f"t{assignment.target_index}")
            expected[key] = expected.get(key, 0) + 1
    observed = {(row["cell_barcode"], row["target_id"]): int(row["count"]) for row in rows}
    assert observed == expected


def test_sensitivity_artifacts_are_deterministic_except_wall_clock(tmp_path: Path) -> None:
    targets = tmp_path / "targets.tsv"
    targets.write_text("target_id\ttarget_seq\nA\tACGTACGT\nB\tTTTTTTTT\n", encoding="utf-8")
    reads = tmp_path / "reads.fastq"
    reads.write_text(
        "@r1\nACGTACGT\n+\nIIIIIIII\n"
        "@r2\nACGTACGA\n+\nIIIIIIII\n"
        "@r3\nTTTTTTTT\n+\nIIIIIIII\n",
        encoding="ascii",
    )
    first = run_sensitivity(
        targets=targets,
        reads=reads,
        target_start=0,
        target_length=8,
        out_dir=tmp_path / "first",
        batch_size=2,
    )
    second = run_sensitivity(
        targets=targets,
        reads=reads,
        target_start=0,
        target_length=8,
        out_dir=tmp_path / "second",
        batch_size=2,
    )

    first.pop("execution", {}).pop("seconds", None)
    second.pop("execution", {}).pop("seconds", None)
    assert first == second
    for name in sorted(path.name for path in (tmp_path / "first").iterdir()):
        if name == "summary.json":
            continue
        assert (tmp_path / "first" / name).read_bytes() == (tmp_path / "second" / name).read_bytes(), name
