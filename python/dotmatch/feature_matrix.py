"""Build deterministic cell-by-feature matrices from pre-extracted observations.

This module deliberately operates on a tabular observation stream rather than
FASTQ pairs.  Each input row must already contain an explicit cell identifier
and a feature sequence window.  It assigns the window to a known feature
library and counts only unique assignments.  Cell calling, barcode correction,
UMI deduplication, and read-pair extraction are outside this command's scope.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence, TextIO

from .core import MATCH_AMBIGUOUS, MATCH_INVALID, MATCH_NONE, MATCH_UNIQUE, Matcher, MatchResult, load_targets, status_name


@dataclass(frozen=True)
class FeatureMatrixResult:
    """Locations and summary data produced by :func:`build_feature_matrix`."""

    output_dir: Path
    summary: dict[str, Any]


@dataclass(frozen=True)
class _Observation:
    observation_id: str
    cell_barcode: str
    sequence: str


def _open_text(path: str | Path, mode: str = "rt") -> TextIO:
    source = Path(path)
    if str(source).endswith(".gz"):
        import gzip

        return gzip.open(source, mode, encoding="utf-8", newline="")
    return source.open(mode, encoding="utf-8", newline="")


def _delimiter(path: str | Path) -> str:
    name = Path(path).name.lower()
    return "," if name.endswith(".csv") or name.endswith(".csv.gz") else "\t"


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require_safe_field(value: str, *, field: str, row_number: int) -> str:
    normalized = value.strip()
    if not normalized:
        raise ValueError(f"row {row_number} has an empty {field}")
    if any(character in normalized for character in "\t\n\r"):
        raise ValueError(f"row {row_number} has a tab or newline in {field}")
    return normalized


def _read_observations(
    path: str | Path,
    *,
    cell_column: str,
    sequence_column: str,
    id_column: str | None,
) -> Iterator[_Observation]:
    source = Path(path)
    with _open_text(source) as fh:
        reader = csv.DictReader(fh, delimiter=_delimiter(source))
        if not reader.fieldnames:
            raise ValueError(f"observation table has no header: {source}")
        fieldnames = {
            field.strip(): field
            for field in reader.fieldnames
            if field is not None and field.strip()
        }
        for field in (cell_column, sequence_column):
            if field not in fieldnames:
                raise ValueError(f"observation table is missing required column '{field}': {source}")
        if id_column is not None and id_column not in fieldnames:
            raise ValueError(f"observation table is missing requested id column '{id_column}': {source}")
        cell_key = fieldnames[cell_column]
        sequence_key = fieldnames[sequence_column]
        id_key = fieldnames[id_column] if id_column is not None else None

        for row_number, row in enumerate(reader, start=2):
            cell = _require_safe_field(row.get(cell_key, "") or "", field=cell_column, row_number=row_number)
            sequence = (row.get(sequence_key, "") or "").strip().upper()
            if any(character in sequence for character in "\t\n\r"):
                raise ValueError(f"row {row_number} has a tab or newline in {sequence_column}")
            if id_column is None:
                observation_id = f"row_{row_number - 1}"
            else:
                observation_id = _require_safe_field(
                    row.get(id_key or id_column, "") or "",
                    field=id_column,
                    row_number=row_number,
                )
            yield _Observation(observation_id=observation_id, cell_barcode=cell, sequence=sequence)


def _chunks(items: Iterable[_Observation], size: int) -> Iterator[list[_Observation]]:
    if size <= 0:
        raise ValueError("batch_size must be positive")
    batch: list[_Observation] = []
    for item in items:
        batch.append(item)
        if len(batch) == size:
            yield batch
            batch = []
    if batch:
        yield batch


def _assign(
    matcher: Matcher,
    sequences: Sequence[str],
    *,
    metric: str,
    k: int,
    ambiguity_policy: str,
) -> list[MatchResult]:
    if metric == "hamming":
        return matcher.assign_hamming(sequences, k=k, policy=ambiguity_policy)
    if metric == "exact":
        return matcher.assign_exact(sequences, policy=ambiguity_policy)
    return matcher.assign(sequences, k=k, policy=ambiguity_policy)


def _write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _write_matrix_market(path: Path, *, counts: Counter[tuple[str, str]], cells: Sequence[str], features: Sequence[str]) -> None:
    cell_index = {cell: index + 1 for index, cell in enumerate(cells)}
    feature_index = {feature: index + 1 for index, feature in enumerate(features)}
    ordered_counts = sorted(
        ((cell_index[cell], feature_index[feature], count) for (cell, feature), count in counts.items()),
        key=lambda row: (row[0], row[1]),
    )
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n")
        fh.write("% DotMatch cell-by-feature unique-assignment counts\n")
        fh.write(f"{len(cells)} {len(features)} {len(ordered_counts)}\n")
        for row, column, count in ordered_counts:
            fh.write(f"{row} {column} {count}\n")


def _new_qc() -> dict[str, int]:
    return {
        "total_observations": 0,
        "assigned_unique": 0,
        "ambiguous": 0,
        "unmatched": 0,
        "invalid": 0,
    }


def _increment_status(qc: dict[str, int], status: int) -> None:
    qc["total_observations"] += 1
    if status == MATCH_UNIQUE:
        qc["assigned_unique"] += 1
    elif status == MATCH_AMBIGUOUS:
        qc["ambiguous"] += 1
    elif status == MATCH_NONE:
        qc["unmatched"] += 1
    else:
        qc["invalid"] += 1


def _validate_targets(targets: Sequence[tuple[str, str]]) -> list[tuple[str, str]]:
    seen: set[str] = set()
    normalized: list[tuple[str, str]] = []
    for target_id, sequence in targets:
        target_id = _require_safe_field(target_id, field="target_id", row_number=len(normalized) + 1)
        if target_id in seen:
            raise ValueError(f"target library contains duplicate target_id '{target_id}'")
        seen.add(target_id)
        normalized.append((target_id, sequence.upper()))
    return normalized


def build_feature_matrix(
    observations: str | Path,
    targets: str | Path,
    output_dir: str | Path,
    *,
    cell_column: str,
    sequence_column: str,
    id_column: str | None = None,
    k: int = 1,
    metric: str = "hamming",
    ambiguity_policy: str = "radius",
    batch_size: int = 4096,
) -> FeatureMatrixResult:
    """Assign a pre-extracted observation table and write a sparse cell matrix.

    ``observations`` must be a headered TSV (or CSV) with an explicit cell
    column and sequence column.  The resulting Matrix Market matrix has cells
    on rows and features on columns.  Only ``unique`` assignments add a count;
    every outcome is retained in ``assignments.tsv`` and ``cell_qc.tsv``.

    The output directory must not already exist.  This avoids mixing artifacts
    from different assignments in one run directory.
    """
    if metric not in {"hamming", "levenshtein", "exact"}:
        raise ValueError("metric must be 'hamming', 'levenshtein', or 'exact'")
    if k < 0:
        raise ValueError("k must be non-negative")
    if metric == "exact" and k != 0:
        raise ValueError("metric='exact' requires k=0")
    if metric == "hamming" and k > 3:
        raise ValueError("hamming assignment supports k between 0 and 3")
    if ambiguity_policy not in {"radius", "best"}:
        raise ValueError("ambiguity_policy must be 'radius' or 'best'")
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")

    observation_path = Path(observations)
    target_path = Path(targets)
    final_dir = Path(output_dir)
    if not observation_path.is_file():
        raise ValueError(f"observation table does not exist: {observation_path}")
    if not target_path.is_file():
        raise ValueError(f"target library does not exist: {target_path}")
    if final_dir.exists():
        raise ValueError(f"output directory already exists: {final_dir}")

    normalized_targets = _validate_targets(load_targets(target_path))
    target_ids = [target_id for target_id, _sequence in normalized_targets]
    target_sequences = [sequence for _target_id, sequence in normalized_targets]
    if metric == "hamming" and len({len(sequence) for sequence in target_sequences}) != 1:
        raise ValueError("hamming feature assignment requires target sequences with one shared length")
    target_by_index = dict(enumerate(normalized_targets))

    final_dir.parent.mkdir(parents=True, exist_ok=True)
    staging_dir = Path(tempfile.mkdtemp(prefix=f".{final_dir.name}.tmp-", dir=final_dir.parent))
    assignments_path = staging_dir / "assignments.tsv"
    counts: Counter[tuple[str, str]] = Counter()
    qc_by_cell: defaultdict[str, dict[str, int]] = defaultdict(_new_qc)
    summary = {
        "schema_version": 1,
        "workflow": "feature_matrix",
        "matrix_orientation": "cells_by_features",
        "metric": metric,
        "k": k,
        "ambiguity_policy": ambiguity_policy,
        "cell_column": cell_column,
        "sequence_column": sequence_column,
        "id_column": id_column,
        "inputs": {
            "observations": str(observation_path),
            "observations_sha256": _sha256(observation_path),
            "targets": str(target_path),
            "targets_sha256": _sha256(target_path),
        },
        "scope": {
            "input": "pre-extracted observations with explicit cell identifiers",
            "cell_calling": "not_performed",
            "barcode_correction": "not_performed",
            "umi_deduplication": "not_performed",
            "paired_read_extraction": "not_performed",
        },
        "total_observations": 0,
        "valid_observations": 0,
        "assigned_unique": 0,
        "assigned_exact": 0,
        "assigned_corrected": 0,
        "ambiguous": 0,
        "unmatched": 0,
        "invalid": 0,
        "features": len(target_ids),
    }

    try:
        with assignments_path.open("w", encoding="utf-8", newline="") as assignment_fh:
            writer = csv.writer(assignment_fh, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "observation_id",
                    "cell_barcode",
                    "observed_seq",
                    "target_id",
                    "target_seq",
                    "distance",
                    "status",
                    "match_count",
                    "second_best_distance",
                ]
            )
            with Matcher(target_sequences) as matcher:
                for batch in _chunks(
                    _read_observations(
                        observation_path,
                        cell_column=cell_column,
                        sequence_column=sequence_column,
                        id_column=id_column,
                    ),
                    batch_size,
                ):
                    valid_positions = [index for index, observation in enumerate(batch) if observation.sequence]
                    valid_sequences = [batch[index].sequence for index in valid_positions]
                    results_by_position: dict[int, MatchResult] = {
                        index: MatchResult(-1, -1, -1, 0, MATCH_INVALID)
                        for index, observation in enumerate(batch)
                        if not observation.sequence
                    }
                    if valid_sequences:
                        results_by_position.update(
                            zip(
                                valid_positions,
                                _assign(
                                    matcher,
                                    valid_sequences,
                                    metric=metric,
                                    k=k,
                                    ambiguity_policy=ambiguity_policy,
                                ),
                            )
                        )

                    for index, observation in enumerate(batch):
                        result = results_by_position[index]
                        status = result.status
                        target_id = ""
                        target_sequence = ""
                        if status == MATCH_UNIQUE and 0 <= result.target_index < len(normalized_targets):
                            target_id, target_sequence = target_by_index[result.target_index]
                        summary["total_observations"] += 1
                        _increment_status(qc_by_cell[observation.cell_barcode], status)
                        if status == MATCH_UNIQUE and target_id:
                            counts[(observation.cell_barcode, target_id)] += 1
                            summary["assigned_unique"] += 1
                            if result.best_distance == 0:
                                summary["assigned_exact"] += 1
                            else:
                                summary["assigned_corrected"] += 1
                        elif status == MATCH_AMBIGUOUS:
                            summary["ambiguous"] += 1
                        elif status == MATCH_NONE:
                            summary["unmatched"] += 1
                        else:
                            summary["invalid"] += 1
                        writer.writerow(
                            [
                                observation.observation_id,
                                observation.cell_barcode,
                                observation.sequence,
                                target_id,
                                target_sequence,
                                result.best_distance,
                                status_name(status),
                                result.match_count,
                                result.second_best_distance,
                            ]
                        )

        cells = sorted(qc_by_cell)
        features = sorted(target_ids)
        target_sequence_by_id = dict(normalized_targets)
        cell_feature_sets: defaultdict[str, set[str]] = defaultdict(set)
        for cell, feature in counts:
            cell_feature_sets[cell].add(feature)

        summary["valid_observations"] = summary["total_observations"] - summary["invalid"]
        summary["cells"] = len(cells)
        summary["nonzero_entries"] = len(counts)
        summary["assignment_rate"] = (
            summary["assigned_unique"] / summary["valid_observations"]
            if summary["valid_observations"]
            else 0.0
        )
        summary["artifacts"] = [
            "assignments.tsv",
            "barcodes.tsv",
            "cell_feature_counts.tsv",
            "cell_qc.tsv",
            "features.tsv",
            "matrix.mtx",
            "summary.json",
        ]

        _write_matrix_market(staging_dir / "matrix.mtx", counts=counts, cells=cells, features=features)
        _write_tsv(staging_dir / "barcodes.tsv", ["cell_barcode"], ((cell,) for cell in cells))
        _write_tsv(
            staging_dir / "features.tsv",
            ["target_id", "target_seq"],
            ((feature, target_sequence_by_id[feature]) for feature in features),
        )
        _write_tsv(
            staging_dir / "cell_feature_counts.tsv",
            ["cell_barcode", "target_id", "count"],
            (
                (cell, feature, count)
                for (cell, feature), count in sorted(counts.items(), key=lambda item: (item[0][0], item[0][1]))
            ),
        )
        _write_tsv(
            staging_dir / "cell_qc.tsv",
            [
                "cell_barcode",
                "total_observations",
                "valid_observations",
                "assigned_unique",
                "ambiguous",
                "unmatched",
                "invalid",
                "unique_features",
                "assignment_rate",
            ],
            (
                (
                    cell,
                    qc_by_cell[cell]["total_observations"],
                    qc_by_cell[cell]["total_observations"] - qc_by_cell[cell]["invalid"],
                    qc_by_cell[cell]["assigned_unique"],
                    qc_by_cell[cell]["ambiguous"],
                    qc_by_cell[cell]["unmatched"],
                    qc_by_cell[cell]["invalid"],
                    len(cell_feature_sets[cell]),
                    (
                        qc_by_cell[cell]["assigned_unique"]
                        / (qc_by_cell[cell]["total_observations"] - qc_by_cell[cell]["invalid"])
                        if qc_by_cell[cell]["total_observations"] - qc_by_cell[cell]["invalid"]
                        else 0.0
                    ),
                )
                for cell in cells
            ),
        )
        with (staging_dir / "summary.json").open("w", encoding="utf-8", newline="") as summary_fh:
            json.dump(summary, summary_fh, indent=2, sort_keys=True)
            summary_fh.write("\n")
        os.replace(staging_dir, final_dir)
    except Exception:
        shutil.rmtree(staging_dir, ignore_errors=True)
        raise

    return FeatureMatrixResult(output_dir=final_dir, summary=summary)
