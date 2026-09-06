"""Expose dependence on observation assumptions; neither response model is biological truth."""
from __future__ import annotations

from typing import Any

from ._version import EXACT_MODEL_VERSION, MODEL_VERSION, SCHEMA_VERSION, __version__
from .engine import analyze
from .io import digest
from .models import Manifest, validated_manifest


def compare_models(manifest: Manifest) -> dict[str, Any]:
    manifest = validated_manifest(manifest)
    results = []
    for model in (MODEL_VERSION, EXACT_MODEL_VERSION):
        data = manifest.model_dump(mode="python")
        data["observation_model"] = model
        results.append(analyze(Manifest.model_validate(data)))
    original, exact = results
    a = {w.hypothesis_id for w in original.witnesses}
    b = {w.hypothesis_id for w in exact.witnesses}
    return {
        "kind": "editwitness.model_comparison", "schema_version": SCHEMA_VERSION,
        "package_version": __version__, "input_manifest_sha256": digest(manifest.model_dump(mode="json")),
        "models": [{"model_version": r.model_version, "analysis_sha256": r.result_sha256,
                    "conclusion": r.conclusion, "equivalent_alternatives": [w.hypothesis_id for w in r.witnesses],
                    "selected_assays": list(r.plan.selected_assays)} for r in results],
        "shared_witnesses": sorted(a & b), "original_only": sorted(a-b), "exact_only": sorted(b-a),
        "witnesses_changed": a != b,
        "panel_changed": original.plan.selected_assays != exact.plan.selected_assays,
        "caveat": "Sensitivity to two idealized observation functions, not validation or a probability of correctness. "
                  "Exact rematching is representation-invariant for modeled sequence observations; it still omits "
                  "mismatch tolerance, nonlocal amplification, sampling and real assay failure.",
    }
