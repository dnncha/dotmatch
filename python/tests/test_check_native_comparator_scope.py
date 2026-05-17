import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_native_comparator_scope.py"


VALID_SCOPE = """# Native Comparator Scope

DotMatch currently has one native alignment-library comparator: Edlib exhaustive global edit-distance assignment. The generated native report records Edlib through `EDLIB_MODE_NW`, `EDLIB_TASK_DISTANCE`, fixed threshold `k`, and zero assignment mismatches before speedups are reported.

SeqAn and Parasail are not part of the checked README, website, or release-note
performance comparison set yet. Before either name belongs in those comparisons,
the repository needs all of the following:

- equivalent global edit-distance or documented semi-global scoring semantics for the exact workload being claimed;
- fixed threshold `k` and identical assignment policy for unique, ambiguous, no-match, and invalid reads;
- native dependency name, version, build flags, and platform in the raw artifact or generated report;
- raw CSV rows under `benchmarks/raw/` plus a generated report under `docs/benchmarks/`;
- zero assignment mismatches against DotMatch and the selected comparator;
- a gate script that fails when only scaffold, smoke, or unmatched-scoring rows are present.

Until that evidence exists, the supported native comparison wording is limited to Edlib exhaustive global edit-distance assignment scans plus the exact-hash and BK-tree baselines recorded in `docs/benchmarks/native/README.md`.
"""


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_native_comparator_scope", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_native_comparator_scope_accepts_documented_boundary(tmp_path):
    checker = _load_checker()
    path = tmp_path / "docs" / "native-comparator-scope.md"
    path.parent.mkdir(parents=True)
    path.write_text(VALID_SCOPE, encoding="utf-8")

    assert checker.check(tmp_path) == []


def test_native_comparator_scope_reports_missing_document(tmp_path):
    checker = _load_checker()

    failures = checker.check(tmp_path)

    assert any("missing native comparator scope document" in failure for failure in failures)


def test_native_comparator_scope_rejects_unbounded_seqan_parasail_claim(tmp_path):
    checker = _load_checker()
    path = tmp_path / "docs" / "native-comparator-scope.md"
    path.parent.mkdir(parents=True)
    path.write_text(
        VALID_SCOPE.replace(
            "SeqAn and Parasail are not part of the checked README, website, or release-note",
            "SeqAn and Parasail comparisons are complete",
        ),
        encoding="utf-8",
    )

    failures = checker.check(tmp_path)

    assert any("SeqAn and Parasail" in failure for failure in failures)


def test_native_comparator_scope_requires_zero_mismatch_evidence(tmp_path):
    checker = _load_checker()
    path = tmp_path / "docs" / "native-comparator-scope.md"
    path.parent.mkdir(parents=True)
    path.write_text(VALID_SCOPE.replace("zero assignment mismatches", "assignment validation"), encoding="utf-8")

    failures = checker.check(tmp_path)

    assert any("zero assignment mismatches" in failure for failure in failures)
