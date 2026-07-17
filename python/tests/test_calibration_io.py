from __future__ import annotations

import csv
import json
from pathlib import Path

from assaycode import cli
from dotmatch.calibration_io import decode_tsv, fit_model_tsv, read_model, write_model


def _training(path: Path) -> Path:
    path.write_text(
        "observed\texpected\tquality\n"
        + "".join("ACGT\tACGT\tIIII\n" for _ in range(20))
        + "AGGT\tACGT\tIIII\n",
        encoding="utf-8",
    )
    return path


def _targets(path: Path) -> Path:
    path.write_text(
        "target_id\ttarget_seq\nexact\tACGT\nneighbor\tAGGT\n",
        encoding="utf-8",
    )
    return path


def _reads(path: Path) -> Path:
    path.write_text(
        "read_id\tobserved\tquality\nr1\tACGT\tIIII\nr2\tANGT\t!!!!\n",
        encoding="utf-8",
    )
    return path


def test_calibration_model_round_trip(tmp_path: Path) -> None:
    model = fit_model_tsv(_training(tmp_path / "training.tsv"))
    output = write_model(model, tmp_path / "model.json")
    restored = read_model(output)

    assert restored == model
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["status"] == "experimental"
    assert payload["schema_version"] == 1


def test_decode_tsv_writes_selective_calls(tmp_path: Path) -> None:
    model = fit_model_tsv(_training(tmp_path / "training.tsv"))
    output = tmp_path / "calls.tsv"

    summary = decode_tsv(
        _reads(tmp_path / "reads.tsv"),
        _targets(tmp_path / "targets.tsv"),
        model,
        output,
    )

    assert summary == {"total": 2, "unique": 1, "ambiguous": 1, "none": 0}
    with output.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["target"] == "ACGT"
    assert rows[0]["status"] == "unique"
    assert rows[1]["target"] == ""
    assert rows[1]["status"] == "ambiguous"


def test_assaycode_calibrate_and_decode_quality(tmp_path: Path, capsys) -> None:
    model = tmp_path / "model.json"
    assert cli.main(
        [
            "calibrate",
            str(_training(tmp_path / "training.tsv")),
            "--out",
            str(model),
        ]
    ) == 0
    calibration_summary = json.loads(capsys.readouterr().out)
    assert calibration_summary["status"] == "experimental"
    assert calibration_summary["cycles"] == 4

    output = tmp_path / "calls.tsv"
    assert cli.main(
        [
            "decode-quality",
            "--reads",
            str(_reads(tmp_path / "reads.tsv")),
            "--targets",
            str(_targets(tmp_path / "targets.tsv")),
            "--model",
            str(model),
            "--out",
            str(output),
        ]
    ) == 0
    decode_summary = json.loads(capsys.readouterr().out)
    assert decode_summary["unique"] == 1
    assert decode_summary["ambiguous"] == 1


def test_training_tsv_requires_explicit_trusted_columns(tmp_path: Path) -> None:
    bad = tmp_path / "bad.tsv"
    bad.write_text("observed\texpected\nACGT\tACGT\n", encoding="utf-8")

    assert cli.main(["calibrate", str(bad), "--out", str(tmp_path / "model.json")]) == 2
