import pytest

from editwitness.io import InputError
from editwitness.models import Allele, DeletionScan, Edit, Manifest
from editwitness.observations import observe_allele
from editwitness.scan import scan_deletions


def test_scan_geometry_matches_reconstruction_for_every_deletion(demo):
    data = demo.model_dump(mode="json")
    data["deletion_scan"] = dict(start_min=0, start_max=890, end_min=1, end_max=900, step=17)
    manifest = Manifest.model_validate(data)
    counts = {a.id: [0, 0, 0] for a in manifest.assays}
    total = blind = 0
    for start in range(0, 891, 17):
        for end in range(1, 901, 17):
            if start >= end:
                continue
            total += 1
            statuses = []
            for a in manifest.assays:
                result = observe_allele(manifest.reference.sequence, Allele(id="d", edits=(Edit(start=start, end=end),)), a)
                index = {"potentially_observable": 0, "original_binding_site_disrupted": 1, "outside_product_bounds": 2}[result.status]
                counts[a.id][index] += 1
                statuses.append(index)
            blind += all(status != 0 for status in statuses)
    result = scan_deletions(manifest)
    assert result.enumerated_deletions == total
    assert result.all_existing_assays_structurally_blind == blind
    for item in result.assays:
        assert [item.potentially_amplifiable, item.binding_site_disrupted, item.outside_product_bounds] == counts[item.assay_id]
    assert len(result.blind_examples) <= 20


def test_scan_bounds_and_empty_grid(demo):
    data = demo.model_dump(mode="json")
    data["deletion_scan"] = None
    with pytest.raises(InputError):
        scan_deletions(Manifest.model_validate(data))
    data["deletion_scan"] = dict(start_min=700, start_max=800, end_min=100, end_max=200)
    with pytest.raises(InputError, match="no valid deletions"):
        scan_deletions(Manifest.model_validate(data))


def test_scan_no_frequency_claim_and_size_bounds(demo):
    data = demo.model_dump(mode="json")
    data["assays"][0]["min_product_bp"] = 300
    data["deletion_scan"].update(min_length=100, max_length=350)
    result = scan_deletions(Manifest.model_validate(data))
    assert result.assays[0].outside_product_bounds > 0
    assert "not event probabilities" in result.caveat
    assert result.enumerated_deletions == sum((result.assays[0].potentially_amplifiable,
                                              result.assays[0].binding_site_disrupted,
                                              result.assays[0].outside_product_bounds))


def test_scan_rejects_huge_work_before_execution():
    with pytest.raises(ValueError, match="500,000"):
        DeletionScan(start_min=0, start_max=1000, end_min=1, end_max=2000)
