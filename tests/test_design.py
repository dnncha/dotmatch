import copy
import json

import pytest

from editwitness import analyze
from editwitness.design import expand_deletions
from editwitness.io import InputError
from editwitness.models import Manifest
from editwitness.sequence import apply_edits


def small_grid(demo, **overrides):
    data = demo.model_dump(mode="json")
    data.update(schema_version="1.1", observation_model="exact-local-sequence-presence-v2")
    data["deletion_scan"] = dict(start_min=195, start_max=201, end_min=215, end_max=221, step=3, **overrides)
    return Manifest.model_validate(data)


def test_expansion_is_explicit_deterministic_and_leaves_input_unchanged(demo):
    original = small_grid(demo)
    before = copy.deepcopy(original.model_dump(mode="json"))
    expanded = expand_deletions(original)
    assert expanded == expand_deletions(original)
    assert original.model_dump(mode="json") == before
    assert expanded.generation.valid_deletions == 9
    assert expanded.generation.added_hypotheses == 9
    assert expanded.generation.added_alleles == 9
    assert len(expanded.hypotheses) == len(original.hypotheses) + 9
    result = analyze(expanded)
    assert result.generation == expanded.generation
    assert len(result.witnesses) > len(analyze(original).witnesses)
    for h in expanded.hypotheses[len(original.hypotheses):]:
        assert h.alleles[0] == "intended"


def test_expansion_fails_without_partial_manifest_when_capped(demo):
    original = small_grid(demo)
    with pytest.raises(InputError, match="No silent sampling"):
        expand_deletions(original, max_new_hypotheses=2)
    assert len(original.alleles) == len(demo.alleles)


def test_expansion_deduplicates_repeat_shifted_local_states(demo):
    data = demo.model_dump(mode="json")
    data["deletion_scan"] = dict(start_min=0, start_max=0, end_min=900, end_max=900)
    original = Manifest.model_validate(data)
    expanded = expand_deletions(original)
    assert expanded.generation.valid_deletions == 1
    assert expanded.generation.deduplicated_states == 1
    assert expanded.generation.added_hypotheses == 0
    assert expanded.alleles == original.alleles


def test_expansion_reuses_existing_deleted_allele_without_duplicate(demo):
    data = demo.model_dump(mode="json")
    data["hypotheses"] = [h for h in data["hypotheses"] if h["id"] != "hidden_window_deletion"]
    data["deletion_scan"] = dict(start_min=0, start_max=0, end_min=900, end_max=900)
    expanded = expand_deletions(Manifest.model_validate(data))
    assert expanded.generation.added_alleles == 0
    assert expanded.generation.added_hypotheses == 1
    assert expanded.hypotheses[-1].alleles == ("intended", "window_deleted")


def test_expansion_avoids_identifier_collision(demo):
    data = small_grid(demo).model_dump(mode="json")
    data["alleles"].append({"id": "ewdel_195_215", "edits": []})
    expanded = expand_deletions(Manifest.model_validate(data))
    assert any(a.id == "ewdel_195_215_1" for a in expanded.alleles)


def test_expansion_requires_grid_and_expected_fixed_allele(demo):
    data = demo.model_dump(mode="json")
    data.pop("deletion_scan")
    with pytest.raises(InputError, match="explicit deletion_scan"):
        expand_deletions(Manifest.model_validate(data))
    with pytest.raises(InputError, match="expected hypothesis"):
        expand_deletions(demo, fixed_allele="not_an_expected_allele")
    for invalid in (0, 1000, True, "100"):
        with pytest.raises(InputError, match="integer"):
            expand_deletions(demo, max_new_hypotheses=invalid)


def test_repeated_expansion_adds_no_duplicate_sequence_states(demo):
    once = expand_deletions(small_grid(demo))
    twice = expand_deletions(once)
    assert twice.hypotheses == once.hypotheses
    assert twice.alleles == once.alleles
    assert twice.generation.added_hypotheses == 0
    assert twice.generation.deduplicated_states == 9


def test_empty_filtered_grid_fails(demo):
    data = demo.model_dump(mode="json")
    data["deletion_scan"] = dict(start_min=200, start_max=201, end_min=100, end_max=101)
    with pytest.raises(InputError, match="no valid deletions"):
        expand_deletions(Manifest.model_validate(data))
