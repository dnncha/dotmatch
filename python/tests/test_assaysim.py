import json

import pytest

from assaycode.cli import command_simulate
from dotmatch.assaysim import simulate_panel


def test_simulation_is_reproducible_and_conserves_reads():
    targets = {"a": "AAAAAAAA", "b": "CCCCCCCC"}
    first = simulate_panel(targets, reads_per_target=100, error_rate=0.05, k=1, seed=42)
    second = simulate_panel(targets, reads_per_target=100, error_rate=0.05, k=1, seed=42)

    assert first == second
    assert first.total_reads == 200
    assert (
        first.correct_unique
        + first.misassigned_unique
        + first.ambiguous
        + first.none
        == first.total_reads
    )
    assert 0.0 <= first.usable_yield <= 1.0
    assert 0.0 <= first.false_discovery_rate <= 1.0


def test_zero_error_recovers_well_separated_targets():
    result = simulate_panel(
        {"a": "AAAA", "b": "CCCC"},
        reads_per_target=10,
        error_rate=0.0,
        k=0,
    )
    assert result.correct_unique == 20
    assert result.misassigned_unique == 0
    assert result.ambiguous == 0
    assert result.none == 0
    assert result.usable_yield == 1.0
    assert result.false_discovery_rate == 0.0


def test_close_targets_expose_ambiguity():
    result = simulate_panel(
        {"a": "AAAA", "b": "AAAT"},
        reads_per_target=5,
        error_rate=0.0,
        k=1,
    )
    assert result.ambiguous == 10
    assert result.usable_yield == 0.0


@pytest.mark.parametrize(
    ("targets", "kwargs"),
    [
        ({}, {}),
        ({"a": "AAAA", "b": "AAAA"}, {}),
        ({"a": "AAAA", "b": "CCC"}, {}),
        ({"a": "AAAN"}, {}),
        ({"a": "AAAA"}, {"reads_per_target": 0}),
        ({"a": "AAAA"}, {"k": -1}),
        ({"a": "AAAA"}, {"error_rate": 1.1}),
        ({"a": "AAAA"}, {"error_rate": [0.1, 0.1]}),
    ],
)
def test_invalid_simulations_are_rejected(targets, kwargs):
    with pytest.raises(ValueError):
        simulate_panel(targets, **kwargs)


def test_simulate_cli_writes_machine_readable_result(tmp_path, capsys):
    targets = tmp_path / "targets.tsv"
    targets.write_text("target_id\ttarget_seq\na\tAAAA\nb\tCCCC\n", encoding="utf-8")
    output = tmp_path / "simulation.json"

    assert command_simulate([
        "--targets", str(targets),
        "--out", str(output),
        "--reads-per-target", "5",
        "--error-rate", "0",
        "-k", "0",
        "--seed", "7",
    ]) == 0

    payload = json.loads(output.read_text(encoding="utf-8"))
    summary = json.loads(capsys.readouterr().out)
    assert payload["correct_unique"] == 10
    assert payload["seed"] == 7
    assert summary["status"] == "experimental"
    assert summary["usable_yield"] == 1.0
