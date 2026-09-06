"""Regression tests for the public-release hardening pass."""
import json
import subprocess
import sys

import pytest

from editwitness import analyze
from editwitness.design import expand_deletions
from editwitness.io import InputError
from editwitness.models import Manifest
from editwitness.selftest import self_test


def heterozygous(demo, reverse=False):
    data = demo.model_dump(mode="json")
    expected = next(h for h in data["hypotheses"] if h["id"] == data["expected_hypothesis"])
    expected["alleles"] = ["intended", "reference"][::(-1 if reverse else 1)]
    data["deletion_scan"] = dict(start_min=195, start_max=195, end_min=215, end_max=215)
    return Manifest.model_validate(data)


@pytest.mark.parametrize("reverse", [False, True])
def test_heterozygous_generation_cannot_choose_by_allele_order(demo, reverse):
    with pytest.raises(InputError, match="explicit --fixed-allele"):
        expand_deletions(heterozygous(demo, reverse))


def test_explicit_fixed_allele_generates_same_challenges_after_order_reversal(demo):
    first = expand_deletions(heterozygous(demo), fixed_allele="intended")
    second = expand_deletions(heterozygous(demo, reverse=True), fixed_allele="intended")
    assert first.alleles == second.alleles
    assert first.hypotheses[-1] == second.hypotheses[-1]
    assert first.generation.fixed_allele == second.generation.fixed_allele == "intended"


def test_sequence_identical_expected_aliases_are_not_false_heterozygotes(demo):
    data = heterozygous(demo).model_dump(mode="json")
    intended = next(a for a in data["alleles"] if a["id"] == "intended")
    data["alleles"].append(dict(intended, id="alias"))
    expected = next(h for h in data["hypotheses"] if h["id"] == data["expected_hypothesis"])
    expected["alleles"] = ["intended", "alias"]
    result = expand_deletions(Manifest.model_validate(data))
    assert result.generation.fixed_allele == "alias"


@pytest.mark.parametrize("bad", [[], 5, True])
def test_fixed_allele_types_fail_cleanly(demo, bad):
    with pytest.raises(InputError, match="expected hypothesis"):
        expand_deletions(demo, fixed_allele=bad)


def test_generation_budget_counts_duplicate_work_before_returning(monkeypatch, demo):
    import editwitness.design as design
    monkeypatch.setattr(design, "MAX_GENERATION_BASES", 0)
    with pytest.raises(InputError, match="reconstruction budget"):
        expand_deletions(heterozygous(demo), fixed_allele="intended")


def test_signal_evidence_budget_is_cumulative_and_exact(monkeypatch, demo):
    import editwitness.engine as engine
    expected = analyze(demo)
    total = sum(len(o.signal_ids) for o in expected.hypothesis_observations)
    monkeypatch.setattr(engine, "MAX_HYPOTHESIS_SIGNAL_REFERENCES", total)
    assert analyze(demo) == expected
    monkeypatch.setattr(engine, "MAX_HYPOTHESIS_SIGNAL_REFERENCES", total - 1)
    with pytest.raises(InputError, match="No partial result"):
        analyze(demo)


def test_self_test_exercises_both_readout_models_and_replay():
    report = self_test()
    assert report["passed"] is True
    assert report["network_used"] is False
    assert len(report["checks"]) == 2
    assert all(len(case["checks"]) == 7 for case in report["checks"])
    assert "not assay sensitivity" in report["scope"]


def test_self_test_reports_runtime_failure_without_false_pass(monkeypatch):
    import editwitness.selftest as module
    def broken(_manifest):
        raise RuntimeError("injected engine failure")
    monkeypatch.setattr(module, "analyze", broken)
    report = module.self_test()
    assert report["passed"] is False
    assert all(c["error_type"] == "RuntimeError" for c in report["checks"])


def test_self_test_cli_is_single_json_with_no_progress_noise():
    result = subprocess.run([sys.executable, "-m", "editwitness", "self-test"],
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0
    assert result.stderr == ""
    assert json.loads(result.stdout)["passed"] is True


def test_self_test_failure_has_distinct_exit_code(monkeypatch, capsys):
    from editwitness.cli import main
    import editwitness.selftest as module
    monkeypatch.setattr(module, "self_test", lambda: {"passed": False})
    assert main(["self-test"]) == 6
    assert json.loads(capsys.readouterr().out) == {"passed": False}
