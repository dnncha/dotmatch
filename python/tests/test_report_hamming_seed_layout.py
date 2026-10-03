import csv
import hashlib
import importlib.util
import json
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("hamming_layout_report", ROOT / "scripts/report_hamming_seed_layout.py")
REPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORT)


def test_report_reproduces_paired_evidence():
    assert REPORT.render() == REPORT.REPORT.read_text()


def mutate_native(tmp_path, monkeypatch, mutate):
    for suffix in ("csv", "json"):
        shutil.copyfile(REPORT.RAW / f"hamming_seed_layout.{suffix}", tmp_path / f"hamming_seed_layout.{suffix}")
    path = tmp_path / "hamming_seed_layout.csv"
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    mutate(rows)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = json.loads(path.with_suffix(".json").read_text())
    metadata["csv_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_suffix(".json").write_text(json.dumps(metadata))
    monkeypatch.setattr(REPORT, "RAW", tmp_path)


def test_gate_rejects_oracle_errors_even_with_matching_file_hash(tmp_path, monkeypatch):
    def mutate(rows):
        next(r for r in rows if r["tool"] == "candidate")["oracle_mismatches"] = "1"
    mutate_native(tmp_path, monkeypatch, mutate)
    with pytest.raises(AssertionError, match="candidate oracle failure"):
        REPORT.load("hamming_seed_layout")


def test_gate_rejects_missing_paired_engine(tmp_path, monkeypatch):
    mutate_native(tmp_path, monkeypatch, lambda rows: rows.pop(0))
    with pytest.raises(AssertionError, match="missing engine/repeat"):
        REPORT.load("hamming_seed_layout")
