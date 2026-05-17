import csv
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "check_barcode_science_readiness.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_barcode_science_readiness", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_current_repo_barcode_science_readiness_passes():
    checker = _load_checker()

    result = checker.audit(ROOT)

    assert result.ok, result.failures
    assert any("public fixed-window evidence datasets" in item for item in result.passed)


def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "tool",
        "workflow",
        "status",
        "n_reads",
        "n_targets",
        "target_start",
        "target_length",
        "k",
        "metric",
        "assigned_unique",
        "assigned_exact",
        "corrected_reads",
        "ambiguous_reads",
        "unmatched_reads",
        "validation_mismatches",
        "validation_notes",
        "exit_code",
        "command",
        "metadata",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _minimal_repo(tmp_path: Path, *, dataset_count: int = 5, mismatch: int = 0, include_baseline: bool = True) -> Path:
    root = tmp_path / "repo"
    (root / "docs").mkdir(parents=True)
    (root / "benchmarks" / "raw").mkdir(parents=True)
    (root / "examples" / "data").mkdir(parents=True)
    (root / "Makefile").write_text(
        "\n".join(f"gate-{i}:\n\ttrue" for i in range(dataset_count)) + "\n",
        encoding="utf-8",
    )
    datasets = []
    for i in range(dataset_count):
        metadata = root / "examples" / "data" / f"metadata_{i}.json"
        metadata.write_text(json.dumps({"evidence_ready": True}), encoding="utf-8")
        raw = root / "benchmarks" / "raw" / f"dataset_{i}.csv"
        rows = [
            {
                "tool": "dotmatch_count",
                "workflow": f"public_dataset_{i}",
                "status": "supported",
                "n_reads": 100,
                "n_targets": 2,
                "target_start": 0,
                "target_length": 4,
                "k": 0,
                "metric": "hamming",
                "assigned_unique": 80,
                "assigned_exact": 80,
                "corrected_reads": 0,
                "ambiguous_reads": 0,
                "unmatched_reads": 20,
                "validation_mismatches": mismatch,
                "validation_notes": "",
                "exit_code": 0,
                "command": "dotmatch count ...",
                "metadata": f"examples/data/metadata_{i}.json",
            }
        ]
        if include_baseline:
            rows.append({
                **rows[0],
                "tool": "exact_slice_hash",
                "metric": "exact",
                "validation_mismatches": 0,
                "validation_notes": "transparent exact substring baseline",
            })
        _write_rows(raw, rows)
        datasets.append(
            {
                "id": f"dataset_{i}",
                "label": f"Dataset {i}",
                "raw_artifact": f"benchmarks/raw/dataset_{i}.csv",
                "metadata": f"examples/data/metadata_{i}.json",
                "gate": f"make gate-{i}",
                "comparator_semantics": "DotMatch k=0 against transparent exact substring baseline.",
                "claim_boundary": "Fixed-window known-target assignment only.",
            }
        )
    (root / "docs" / "barcode-science-readiness.json").write_text(
        json.dumps({"schema_version": 1, "datasets": datasets}),
        encoding="utf-8",
    )
    return root


def test_barcode_science_readiness_requires_five_public_datasets(tmp_path):
    checker = _load_checker()
    root = _minimal_repo(tmp_path, dataset_count=4)

    result = checker.audit(root)

    assert not result.ok
    assert any("at least 5" in failure for failure in result.failures)


def test_barcode_science_readiness_requires_comparator_baseline(tmp_path):
    checker = _load_checker()
    root = _minimal_repo(tmp_path, include_baseline=False)

    result = checker.audit(root)

    assert not result.ok
    assert any("transparent exact baseline" in failure for failure in result.failures)


def test_barcode_science_readiness_rejects_validation_mismatches(tmp_path):
    checker = _load_checker()
    root = _minimal_repo(tmp_path, mismatch=1)

    result = checker.audit(root)

    assert not result.ok
    assert any("validation_mismatches" in failure for failure in result.failures)
