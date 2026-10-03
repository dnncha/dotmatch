import csv
import copy
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("guide_counter_report", ROOT / "scripts/report_guide_counter.py")
REPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORT)


def test_report_reproduces_paired_evidence():
    metadata, rows, stats = REPORT.load()
    assert len(rows) == 80
    assert REPORT.render(metadata, rows, stats) == REPORT.REPORT.read_text()


@pytest.mark.parametrize("mutation", [
    lambda record: record.update(figures={}),
    lambda record: record.update(commands=0),
    lambda record: record.update(count_mismatches=1),
    lambda record: record['figures'].pop('guide_counter_memory.svg'),
    lambda record: record['figures'].update({'../unexpected.svg': '0' * 64}),
])
def test_gate_requires_complete_validation_record(mutation):
    _, rows, _ = REPORT.load()
    recorded = json.loads((REPORT.RAW / f'{REPORT.STEM}_validation.json').read_text())
    REPORT.check_validation_record(recorded, rows)
    damaged = copy.deepcopy(recorded)
    mutation(damaged)
    with pytest.raises(AssertionError):
        REPORT.check_validation_record(damaged, rows)


def test_gate_rejects_stale_website_headline(tmp_path, monkeypatch):
    _, _, stats = REPORT.load()
    source = (ROOT / 'app/guide-counter-benchmark.tsx').read_text()
    path = tmp_path / 'app/guide-counter-benchmark.tsx'
    path.parent.mkdir()
    mismatch = [s for s in stats if s['k'] == 1 and s['case']['scope'] == 'controlled_simulation']
    speed = f"{min(s['speedup'] for s in mismatch):.1f}–{max(s['speedup'] for s in mismatch):.1f}×"
    path.write_text(source.replace(speed, '9.9–99.9×'))
    monkeypatch.setattr(REPORT, 'ROOT', tmp_path)
    monkeypatch.setattr(REPORT, 'CALLOUTS', ('app/guide-counter-benchmark.tsx',))
    with pytest.raises(AssertionError, match='headline differs'):
        REPORT.check_callout_headlines(stats)
    path.write_text(source)
    REPORT.check_callout_headlines(stats)


def mutate_evidence(tmp_path, monkeypatch, mutate):
    for suffix in ("csv", "json"):
        name = f"{REPORT.STEM}.{suffix}"
        shutil.copyfile(REPORT.RAW / name, tmp_path / name)
    path = tmp_path / f"{REPORT.STEM}.csv"
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


def test_gate_rejects_count_errors_even_with_matching_file_hash(tmp_path, monkeypatch):
    mutate_evidence(tmp_path, monkeypatch, lambda rows: rows[0].update(count_mismatches="1"))
    with pytest.raises(AssertionError, match="output disagreement"):
        REPORT.load()


def test_gate_rejects_missing_paired_run(tmp_path, monkeypatch):
    mutate_evidence(tmp_path, monkeypatch, lambda rows: rows.pop(0))
    with pytest.raises(AssertionError, match="missing paired run"):
        REPORT.load()


def test_gate_rejects_different_offset_parameters(tmp_path, monkeypatch):
    def mutate(rows):
        command = json.loads(rows[0]["command"])
        command[command.index("--offset-min-fraction") + 1] = ".5"
        rows[0]["command"] = json.dumps(command)
    mutate_evidence(tmp_path, monkeypatch, mutate)
    with pytest.raises(AssertionError):
        REPORT.load()


def test_gate_rejects_different_full_count_digest(tmp_path, monkeypatch):
    mutate_evidence(tmp_path, monkeypatch, lambda rows: rows[0].update(counts_sha256="0" * 64))
    with pytest.raises(AssertionError, match="count digests differ"):
        REPORT.load()


def test_optimized_python_still_rejects_corrupt_evidence(tmp_path, monkeypatch):
    mutate_evidence(tmp_path, monkeypatch, lambda rows: rows[0].update(count_mismatches="1"))
    code = (f"import runpy; from pathlib import Path; "
            f"m=runpy.run_path({str(ROOT / 'scripts/report_guide_counter.py')!r}); "
            f"m['load'].__globals__['RAW']=Path({str(tmp_path)!r}); m['load']()")
    result = subprocess.run([sys.executable, "-O", "-c", code], capture_output=True, text=True)
    assert result.returncode != 0
    assert "output disagreement" in result.stderr
