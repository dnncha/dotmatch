"""Property-based invariants for the assignment core and streaming adapter."""

from __future__ import annotations

import tempfile
from pathlib import Path

import dotmatch as dm
from dotmatch.panel import hamming_distance
from hypothesis import given, settings
from hypothesis import strategies as st

ALPHABET = "ACGT"
DNA = st.text(alphabet=ALPHABET, min_size=1, max_size=12)
TARGET_LISTS = st.lists(
    st.text(alphabet=ALPHABET, min_size=1, max_size=12),
    min_size=1,
    max_size=8,
)
READ_LISTS = st.lists(
    st.text(alphabet=ALPHABET, min_size=1, max_size=12),
    min_size=1,
    max_size=8,
)


def _equal_length_pairs() -> st.SearchStrategy[tuple[str, str]]:
    return st.integers(min_value=1, max_value=12).flatmap(
        lambda length: st.tuples(
            st.text(alphabet=ALPHABET, min_size=length, max_size=length),
            st.text(alphabet=ALPHABET, min_size=length, max_size=length),
        )
    )


@settings(max_examples=100, deadline=None)
@given(DNA, DNA)
def test_edit_distance_is_symmetric_and_threshold_matches(a: str, b: str) -> None:
    distance = dm.distance(a, b)

    assert distance == dm.distance(b, a)
    for k in range(0, 5):
        assert dm.distance_leq(a, b, k) is (distance <= k)


@settings(max_examples=100, deadline=None)
@given(_equal_length_pairs())
def test_hamming_equals_levenshtein_for_equal_length_substitutions(
    pair: tuple[str, str],
) -> None:
    left, right = pair

    assert dm.distance(left, right) <= hamming_distance(left, right)


@settings(max_examples=60, deadline=None)
@given(TARGET_LISTS, READ_LISTS)
def test_exact_and_hamming_k0_agree(targets: list[str], reads: list[str]) -> None:
    exact = dm.assign_exact(reads, targets)
    hamming = dm.assign_hamming(reads, targets, k=0)

    assert [
        (result.status, result.target_index, result.best_distance, result.match_count)
        for result in exact
    ] == [
        (result.status, result.target_index, result.best_distance, result.match_count)
        for result in hamming
    ]


@settings(max_examples=80, deadline=None)
@given(TARGET_LISTS, DNA, st.integers(min_value=0, max_value=2))
def test_radius_is_never_more_permissive_than_best(
    targets: list[str], read: str, k: int
) -> None:
    best = dm.assign([read], targets, k=k, policy="best")[0]
    radius = dm.assign([read], targets, k=k, policy="radius")[0]

    if best.status in {dm.MATCH_INVALID, dm.MATCH_NONE}:
        assert radius.status == best.status
    elif best.status == dm.MATCH_UNIQUE:
        assert radius.status in {dm.MATCH_UNIQUE, dm.MATCH_AMBIGUOUS}
        if radius.status == dm.MATCH_UNIQUE:
            assert radius.target_index == best.target_index


def _fixed_length_case() -> st.SearchStrategy[tuple[list[str], list[str]]]:
    return st.integers(min_value=1, max_value=12).flatmap(
        lambda length: st.tuples(
            st.lists(
                st.text(alphabet=ALPHABET, min_size=length, max_size=length),
                min_size=1,
                max_size=8,
            ),
            st.lists(
                st.text(alphabet=ALPHABET, min_size=length, max_size=length),
                min_size=1,
                max_size=8,
            ),
        )
    )


@settings(max_examples=30, deadline=None)
@given(_fixed_length_case())
def test_stream_batch_size_does_not_change_assignments(
    case: tuple[list[str], list[str]]
) -> None:
    targets, reads = case
    with tempfile.TemporaryDirectory(prefix="dotmatch-invariants-") as temporary:
        fastq = Path(temporary) / "reads.fastq"
        fastq.write_text(
            "".join(
                f"@r{index}\n{read}\n+\n{'I' * len(read)}\n"
                for index, read in enumerate(reads)
            ),
            encoding="ascii",
        )
        one = list(dm.stream_assign(fastq, targets, metric="hamming", k=1, batch_size=1))
        many = list(dm.stream_assign(fastq, targets, metric="hamming", k=1, batch_size=17))

    assert one == many
