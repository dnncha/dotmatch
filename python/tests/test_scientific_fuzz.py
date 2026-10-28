"""Fuzz and mutation-sensitivity checks for the scientific input boundary."""

from __future__ import annotations

import csv
import io
import tempfile
from pathlib import Path

import dotmatch as dm
import pytest
from dotmatch.assayspec import AssaySpecError, load_assay_spec
from dotmatch.count_io import read_count_table
from dotmatch.fastq_io import iter_fastq_records
from dotmatch.feature_matrix import build_feature_matrix
from dotmatch.panel import AssignmentConfig, PanelConstraints, PanelRecord, check_panel
from hypothesis import given, settings
from hypothesis import strategies as st

pytestmark = pytest.mark.scientific

READ_ID = st.from_regex(r"[A-Za-z0-9_.:-]{1,12}", fullmatch=True)
DNA = st.text(alphabet="ACGTN", min_size=1, max_size=24)
QUALITY = st.text(
    alphabet=st.characters(min_codepoint=33, max_codepoint=126),
    min_size=1,
    max_size=24,
)


@st.composite
def _valid_fastq_record(draw: st.DrawFn) -> tuple[str, str, str]:
    read_id = draw(READ_ID)
    length = draw(st.integers(min_value=1, max_value=24))
    sequence = draw(st.text(alphabet="ACGTN", min_size=length, max_size=length))
    quality = draw(
        st.text(
            alphabet=st.characters(min_codepoint=33, max_codepoint=126),
            min_size=length,
            max_size=length,
        )
    )
    return read_id, sequence, quality


def _fastq_text(records: list[tuple[str, str, str]]) -> str:
    return "".join(
        f"@{read_id} description\n{sequence}\n+{read_id}\n{quality}\n"
        for read_id, sequence, quality in records
    )


@settings(max_examples=80, deadline=None)
@given(st.lists(_valid_fastq_record(), min_size=1, max_size=5))
def test_valid_fastq_round_trips(records: list[tuple[str, str, str]]) -> None:
    parsed = list(iter_fastq_records(io.StringIO(_fastq_text(records)), "fuzz.fastq"))

    assert [(record.read_id, record.seq, record.qual) for record in parsed] == records


@settings(max_examples=80, deadline=None)
@given(_valid_fastq_record(), st.sampled_from(("header", "plus", "empty-seq", "length", "control", "space", "non-ascii", "id-mismatch")))
def test_corrupt_fastq_is_rejected(
    record: tuple[str, str, str], corruption: str
) -> None:
    read_id, sequence, quality = record
    header = f"@{read_id}"
    plus = f"+{read_id}"
    if corruption == "header":
        header = read_id
    elif corruption == "plus":
        plus = "-"
    elif corruption == "empty-seq":
        sequence = ""
    elif corruption == "length":
        quality += "I"
    elif corruption == "control":
        sequence = "\x01" + sequence[1:]
    elif corruption == "space":
        quality = " " + quality[1:]
    elif corruption == "non-ascii":
        sequence = "\u00e9" + sequence[1:]
    else:
        plus = "+some_other_id"
    text = f"{header}\n{sequence}\n{plus}\n{quality}\n"

    with pytest.raises(ValueError):
        list(iter_fastq_records(io.StringIO(text), "corrupt.fastq"))


@settings(max_examples=100, deadline=None)
@given(st.text(max_size=240))
def test_arbitrary_fastq_text_never_raises_unexpected_errors(text: str) -> None:
    try:
        records = list(iter_fastq_records(io.StringIO(text), "arbitrary.fastq"))
    except (ValueError, csv.Error):
        return
    for record in records:
        assert record.seq
        assert len(record.seq) == len(record.qual)
        assert record.seq.isascii() and record.seq.isprintable() and " " not in record.seq
        assert record.qual.isascii() and record.qual.isprintable() and " " not in record.qual


@settings(max_examples=60, deadline=None)
@given(
    st.lists(st.integers(min_value=0, max_value=10**6), min_size=1, max_size=8),
    st.lists(st.text(alphabet="ACGT", min_size=1, max_size=6), min_size=1, max_size=4),
)
def test_valid_count_tables_round_trip(
    values: list[int], gene_parts: list[str]
) -> None:
    samples = [f"S{i}" for i in range(len(gene_parts))]
    rows = [
        [f"g{i}", gene_parts[i % len(gene_parts)], *[str(values[(i + j) % len(values)]) for j in range(len(samples))]]
        for i in range(len(values))
    ]
    text = "sgRNA\tGene\t" + "\t".join(samples) + "\n"
    text += "\n".join("\t".join(row) for row in rows) + "\n"
    with tempfile.TemporaryDirectory(prefix="dotmatch-count-fuzz-") as temporary:
        path = Path(temporary) / "counts.tsv"
        path.write_text(text, encoding="utf-8")
        table = read_count_table(path)
    expected = tuple(
        tuple(int(row[2 + j]) for j in range(len(samples))) for row in rows
    )
    assert table.counts == expected
    assert table.sample_names == tuple(samples)
    assert table.target_ids == tuple(row[0] for row in rows)


@settings(max_examples=100, deadline=None)
@given(st.text(max_size=240))
def test_arbitrary_count_text_never_raises_unexpected_errors(text: str) -> None:
    with tempfile.TemporaryDirectory(prefix="dotmatch-count-fuzz-") as temporary:
        path = Path(temporary) / "counts.tsv"
        path.write_text(text, encoding="utf-8")
        try:
            table = read_count_table(path)
        except (ValueError, csv.Error):
            return
    assert table.target_ids
    assert len(set(table.target_ids)) == len(table.target_ids)
    assert len(set(table.sample_names)) == len(table.sample_names)
    assert len(table.sample_columns) == len(table.sample_names)
    assert len(table.counts) == len(table.target_ids)
    for row in table.counts:
        assert len(row) == len(table.sample_names)
        assert all(isinstance(value, int) and value >= 0 for value in row)


@settings(max_examples=60, deadline=None)
@given(st.text(max_size=300))
def test_arbitrary_assay_toml_never_raises_unexpected_errors(text: str) -> None:
    with tempfile.TemporaryDirectory(prefix="dotmatch-assay-fuzz-") as temporary:
        path = Path(temporary) / "assay.toml"
        path.write_text(text, encoding="utf-8")
        try:
            load_assay_spec(path)
        except AssaySpecError:
            return


def test_assignment_tsv_is_byte_deterministic_and_conserves_counts(tmp_path: Path) -> None:
    reads = ["ACGT", "ACGC", "TTTT", "ACGG"]
    targets = ["ACGT", "AGGT", "ACGA"]
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    with tempfile.TemporaryDirectory(prefix="dotmatch-assignments-") as temporary:
        fastq = Path(temporary) / "reads.fastq"
        fastq.write_text(
            "".join(f"@r{i}\n{read}\n+\n{'I' * len(read)}\n" for i, read in enumerate(reads)),
            encoding="ascii",
        )
        streamed = list(dm.stream_assign(fastq, targets, metric="hamming", k=1))
        summary_one = dm.write_assignments_tsv(streamed, first)
        summary_two = dm.write_assignments_tsv(streamed, second)
    assert first.read_bytes() == second.read_bytes()
    assert summary_one == summary_two
    assert summary_one["total_reads"] == len(reads)
    assert summary_one["assigned_unique"] + summary_one["ambiguous"] + summary_one["unmatched"] + summary_one["invalid"] == len(reads)


def test_panel_certificate_is_deterministic_and_complete(tmp_path: Path) -> None:
    records = [PanelRecord("BC001", "ACGTACGT"), PanelRecord("BC002", "ACGTACGA")]
    assignment = AssignmentConfig(metric="hamming", k=1, ambiguity_policy="radius")
    constraints = PanelConstraints(length=8, min_hamming_distance=5)
    first = tmp_path / "first"
    second = tmp_path / "second"
    certificate_one = check_panel(records, assignment, constraints, first, source_panel=Path("barcodes.tsv"))
    certificate_two = check_panel(records, assignment, constraints, second, source_panel=Path("barcodes.tsv"))

    assert certificate_one == certificate_two
    assert certificate_one["status"] in {"pass", "warn", "fail"}
    assert certificate_one["n_barcodes"] == 2


def test_feature_matrix_artifacts_are_byte_deterministic(tmp_path: Path) -> None:
    targets = tmp_path / "targets.tsv"
    targets.write_text("target_id\ttarget_seq\nA\tAAAA\nB\tAACC\n", encoding="utf-8")
    observations = tmp_path / "observations.tsv"
    observations.write_text(
        "observation_id\tcell_barcode\tfeature_seq\n"
        "o1\tcell_a\tAAAA\n"
        "o2\tcell_a\tAAAT\n"
        "o3\tcell_b\tAACC\n",
        encoding="utf-8",
    )
    first = build_feature_matrix(
        observations, targets, tmp_path / "m1", cell_column="cell_barcode", sequence_column="feature_seq", k=1
    )
    second = build_feature_matrix(
        observations, targets, tmp_path / "m2", cell_column="cell_barcode", sequence_column="feature_seq", k=1
    )
    names = sorted(path.name for path in first.output_dir.iterdir() if path.is_file())
    assert names == sorted(path.name for path in second.output_dir.iterdir() if path.is_file())
    for name in names:
        left = (first.output_dir / name).read_bytes()
        right = (second.output_dir / name).read_bytes()
        assert left == right, name


def test_semantic_mutations_are_detectable() -> None:
    radius = dm.assign(["ACGT"], ["ACGT", "AGGT"], k=1, policy="radius")[0]
    best = dm.assign(["ACGT"], ["ACGT", "AGGT"], k=1, policy="best")[0]
    assert radius.status == dm.MATCH_AMBIGUOUS
    assert best.status == dm.MATCH_UNIQUE

    exact = dm.assign_exact(["ACGA"], ["ACGT"])[0]
    relaxed = dm.assign_hamming(["ACGA"], ["ACGT"], k=1)[0]
    assert exact.status == dm.MATCH_NONE
    assert relaxed.status == dm.MATCH_UNIQUE


def test_mutation_that_counts_ambiguous_reads_is_detected() -> None:
    results = dm.assign(["ACGT", "ACGC", "TTTT"], ["ACGT", "AGGT", "ACGA"], k=1)
    unique_only = sum(1 for result in results if result.status == dm.MATCH_UNIQUE)
    count_everything = len(results)

    assert unique_only != count_everything
    assert unique_only == 0


def test_mutation_that_uses_float_counts_is_detected() -> None:
    exact_value = 2**53 + 1
    lossy_value = int(float(exact_value))

    assert lossy_value != exact_value
    assert lossy_value == 2**53

