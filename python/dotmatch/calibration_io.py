"""File interfaces for experimental calibrated decoding."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Mapping

from .calibration import ErrorModel, decode, fit_error_model


def fit_model_tsv(path: str | Path, *, prior_strength: float = 100.0) -> ErrorModel:
    rows = _dict_rows(path)
    observations: list[tuple[str, str, str]] = []
    for line_number, row in enumerate(rows, start=2):
        try:
            observed = row["observed"].strip().upper()
            expected = row["expected"].strip().upper()
            quality = row["quality"].strip()
        except KeyError as exc:
            raise ValueError("training TSV requires observed, expected, and quality columns") from exc
        if not observed or not expected or not quality:
            raise ValueError(f"training TSV line {line_number} contains an empty required value")
        observations.append((observed, expected, quality))
    if not observations:
        raise ValueError("training TSV contains no observations")
    return fit_error_model(observations, prior_strength=prior_strength)


def write_model(model: ErrorModel, path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "status": "experimental",
        **model.to_dict(),
    }
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return output


def read_model(path: str | Path) -> ErrorModel:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("schema_version") != 1:
        raise ValueError("unsupported calibration model schema")
    substitutions = data.get("substitution_counts")
    if not isinstance(substitutions, dict):
        raise ValueError("calibration model substitution_counts must be an object")
    normalized: dict[str, dict[str, int]] = {}
    for truth, row in substitutions.items():
        if not isinstance(truth, str) or not isinstance(row, dict):
            raise ValueError("invalid calibration substitution row")
        normalized[truth] = {str(called): int(count) for called, count in row.items()}
    return ErrorModel(
        cycle_totals=tuple(int(value) for value in data["cycle_totals"]),
        cycle_errors=tuple(int(value) for value in data["cycle_errors"]),
        substitution_counts=normalized,
        prior_strength=float(data["prior_strength"]),
    )


def decode_tsv(
    reads_path: str | Path,
    targets_path: str | Path,
    model: ErrorModel,
    output_path: str | Path,
    *,
    posterior_min: float = 0.99,
    likelihood_ratio_min: float = 10.0,
) -> dict[str, int]:
    targets = _read_targets(targets_path)
    rows = _dict_rows(reads_path)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    summary = {"total": 0, "unique": 0, "ambiguous": 0, "none": 0}
    with output.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            lineterminator="\n",
            fieldnames=[
                "read_id",
                "observed",
                "target",
                "status",
                "posterior",
                "second_posterior",
                "likelihood_ratio",
            ],
        )
        writer.writeheader()
        for line_number, row in enumerate(rows, start=2):
            observed = str(row.get("observed", "")).strip().upper()
            quality = str(row.get("quality", "")).strip()
            read_id = str(row.get("read_id", f"read_{line_number - 2}")).strip()
            if not observed or not quality:
                raise ValueError(f"reads TSV line {line_number} requires observed and quality")
            call = decode(
                observed,
                quality,
                targets,
                model,
                posterior_min=posterior_min,
                likelihood_ratio_min=likelihood_ratio_min,
            )
            summary["total"] += 1
            summary[call.status] += 1
            writer.writerow(
                {
                    "read_id": read_id,
                    "observed": observed,
                    "target": call.target or "",
                    "status": call.status,
                    "posterior": f"{call.posterior:.12g}",
                    "second_posterior": f"{call.second_posterior:.12g}",
                    "likelihood_ratio": "inf"
                    if call.likelihood_ratio == float("inf")
                    else f"{call.likelihood_ratio:.12g}",
                }
            )
    return summary


def _dict_rows(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing TSV header")
        return [dict(row) for row in reader]


def _read_targets(path: str | Path) -> list[str]:
    rows = _dict_rows(path)
    if not rows:
        raise ValueError("target table contains no targets")
    candidates = ["target_seq", "guide_seq", "barcode_seq", "sequence", "seq"]
    fieldnames = list(rows[0])
    sequence_column = next((column for column in candidates if column in fieldnames), None)
    if sequence_column is None:
        if len(fieldnames) < 2:
            raise ValueError("target table requires a sequence column")
        sequence_column = fieldnames[1]
    targets = [str(row.get(sequence_column, "")).strip().upper() for row in rows]
    if any(not target for target in targets):
        raise ValueError("target table contains an empty sequence")
    return targets
