"""Offline installation checks using explicitly synthetic bundled fixtures.

This checks the installed executable, model semantics and deterministic replay.
It is not independent scientific validation and never evaluates a user's sample.
"""
from __future__ import annotations

import json
from importlib.resources import files
from typing import Any

from ._version import EXACT_MODEL_VERSION, __version__
from .engine import analyze
from .io import digest
from .models import Manifest

_CASES = (
    ("demo.json", {"hidden_primer_deletion", "hidden_window_deletion"}),
    ("paired_end.json", {"intended_reference", "hidden_primer_deletion",
                         "hidden_window_deletion", "interior_deletion"}),
)


def self_test() -> dict[str, Any]:
    """Return machine-readable checks; a caught failure is never called a pass."""
    checks: list[dict[str, Any]] = []
    for filename, expected_witnesses in _CASES:
        try:
            data = json.loads(files("editwitness").joinpath("data", filename).read_text(encoding="utf-8"))
            manifest = Manifest.model_validate(data)
            first = analyze(manifest)
            replay = analyze(manifest)
            actual = {w.hypothesis_id for w in first.witnesses}
            assertions = {
                "synthetic_fixture": manifest.reference.synthetic,
                "exact_local_model": first.model_version == EXACT_MODEL_VERSION,
                "known_counterexamples": actual == expected_witnesses,
                "outer_assay_selected": first.plan.selected_assays == ("outer",),
                "whole_window_stays_unresolved": first.plan.unresolved_hypotheses == ("hidden_window_deletion",),
                "checksum_matches": first.result_sha256 == digest(first.model_dump(mode="json", exclude={"result_sha256"})),
                "replay_matches": first.result_sha256 == replay.result_sha256,
            }
            checks.append({"fixture": filename, "passed": all(assertions.values()),
                           "checks": assertions, "analysis_sha256": first.result_sha256})
        except Exception as error:
            # An installation diagnostic must report a broken engine, not mask it.
            # KeyboardInterrupt and SystemExit are not intercepted.
            checks.append({"fixture": filename, "passed": False,
                           "error_type": type(error).__name__, "message": str(error)[:500]})
    return {
        "kind": "editwitness.software_self_test", "schema_version": "1.0",
        "package_version": __version__, "passed": all(c["passed"] for c in checks),
        "network_used": False, "checks": checks,
        "scope": "Bundled synthetic software checks, not assay sensitivity or biological validation.",
    }
