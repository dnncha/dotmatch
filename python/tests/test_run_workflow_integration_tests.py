from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "scripts" / "run_workflow_integration_tests.py"


def _load_runner():
    spec = importlib.util.spec_from_file_location("run_workflow_integration_tests", RUNNER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_planemo_lints_all_wrappers_and_runs_the_scoped_crispr_test(monkeypatch) -> None:
    runner = _load_runner()
    calls: list[tuple[list[str], Path, dict[str, str]]] = []

    monkeypatch.setattr(runner, "_tool", lambda name: name)
    monkeypatch.setattr(
        runner,
        "_run",
        lambda command, *, cwd, env: calls.append((command, cwd, env)),
    )

    environment = {"PATH": ""}
    runner.run_planemo(environment)

    assert calls == [
        (["planemo", "lint", *[str(path) for path in runner.GALAXY_WRAPPERS]], runner.ROOT, environment),
        (
            ["planemo", "test", "--install_galaxy", str(runner.GALAXY_CRISPR_WRAPPER)],
            runner.ROOT,
            environment,
        ),
    ]
