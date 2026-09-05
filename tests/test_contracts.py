import copy
import json
from pathlib import Path

import jsonschema
import pytest
from pydantic import ValidationError

from editwitness import analyze
from editwitness.cli import schema_for
from editwitness.models import Assay, DeletionScan, Edit, Interval, Manifest
from editwitness.scan import scan_deletions
from editwitness.sequence import reverse_complement


def mutate(data, path, value):
    for field in path[:-1]:
        data = data[field]
    data[path[-1]] = value


@pytest.mark.parametrize("path,value", [
    (("schema_version",), "2.0"), (("coordinate_system",), "1-based"),
    (("unexpected",), 1), (("reference", "sequence"), "ACNT"),
    (("reference", "sequence"), "acgt"), (("reference", "synthetic"), "true"),
    (("assays", 0, "cost_units"), True), (("assays", 0, "cost_units"), 1.5),
    (("assays", 0, "cost_units"), "3"), (("assays", 0, "cost_units"), 0),
    (("assays", 0, "readout"), "unknown"), (("assays", 0, "read_bases"), 100),
    (("assays", 0, "left_primer", "start"), -1),
    (("assays", 0, "right_primer", "end"), 9999),
    (("assays", 0, "right_primer", "start"), 210),
    (("assays", 0, "max_product_bp"), 0),
    (("alleles", 1, "edits", 0, "start"), 1000),
    (("alleles", 1, "edits", 0, "end"), 1000),
    (("alleles", 1, "id"), "reference"), (("alleles", 1, "id"), "bad/name"),
    (("hypotheses", 0, "alleles"), ["missing", "reference"]),
    (("hypotheses", 0, "alleles"), ["intended"]),
    (("expected_hypothesis",), "missing"),
    (("candidates", 0, "id"), "inner"),
    (("deletion_scan", "end_max"), 1000),
    (("deletion_scan", "start_max"), 900),
    (("assays", 0, "left_oligo"), "TT"),
    (("assays", 0, "right_oligo"), "TT"),
])
def test_strict_invalid_manifest_cases(demo, path, value):
    data = demo.model_dump(mode="json")
    mutate(data, path, value)
    with pytest.raises(ValidationError):
        Manifest.model_validate(data)


def test_sorted_nonoverlapping_edits_and_noop_rejected(demo):
    for edits in (
        [{"start": 400, "end": 450}, {"start": 440, "end": 460}],
        [{"start": 450, "end": 451, "sequence": "A"}, {"start": 400, "end": 410}],
        [{"start": 400, "end": 400, "sequence": "A"}, {"start": 400, "end": 401, "sequence": "T"}],
        [{"start": 450, "end": 451, "sequence": demo.reference.sequence[450]}],
    ):
        data = demo.model_dump(mode="json")
        data["alleles"][1]["edits"] = edits
        with pytest.raises(ValidationError):
            Manifest.model_validate(data)


def test_optional_oligos_use_5_prime_orientation(demo):
    data = demo.model_dump(mode="json")
    a = data["assays"][0]
    seq = demo.reference.sequence
    a["left_oligo"] = seq[200:220]
    a["right_oligo"] = reverse_complement(seq[680:700])
    Manifest.model_validate(data)
    a["right_oligo"] = seq[680:700]
    with pytest.raises(ValidationError, match="reverse complement"):
        Manifest.model_validate(data)


def test_nested_models_are_frozen_and_sequences_are_tuples(demo):
    assert isinstance(demo.alleles, tuple)
    with pytest.raises(ValidationError):
        demo.reference.sequence = "A"
    with pytest.raises(ValidationError):
        demo.assays[0].cost_units = 10


def test_basic_interval_edit_and_assay_invariants():
    with pytest.raises(ValidationError):
        Interval(start=10, end=5)
    with pytest.raises(ValidationError):
        Edit(start=4, end=4)
    with pytest.raises(ValidationError):
        Edit(start=5, end=4)
    with pytest.raises(ValidationError):
        Assay(id="a", left_primer=Interval(start=0, end=2), right_primer=Interval(start=8, end=10), readout="paired_end")
    with pytest.raises(ValidationError):
        Assay(id="a", left_primer=Interval(start=0, end=2), right_primer=Interval(start=8, end=10), min_product_bp=100, max_product_bp=10)
    for kwargs in [dict(start_min=10, start_max=5, end_min=1, end_max=20),
                   dict(start_min=0, start_max=5, end_min=1, end_max=20, min_length=9, max_length=2)]:
        with pytest.raises(ValidationError):
            DeletionScan(**kwargs)


def test_large_insertions_count_toward_work_budget(demo):
    data = demo.model_dump(mode="json")
    data["alleles"][1]["edits"] = [{"start": i, "end": i, "sequence": "A" * 20000} for i in range(50)]
    data["candidates"] = [dict(data["candidates"][0], id=f"c{i}") for i in range(24)]
    with pytest.raises(ValidationError, match="work budget"):
        Manifest.model_validate(data)


def test_json_schemas_validate_real_payloads(demo):
    for kind, value in [("manifest", demo.model_dump(mode="json")),
                        ("analysis", analyze(demo).model_dump(mode="json")),
                        ("scan", scan_deletions(demo).model_dump(mode="json"))]:
        schema = schema_for(kind)
        jsonschema.Draft202012Validator.check_schema(schema)
        jsonschema.validate(value, schema)
        invalid = copy.deepcopy(value)
        invalid["unknown_contract_field"] = 1
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(invalid, schema)


def test_checked_in_schemas_match_runtime():
    root = Path(__file__).resolve().parents[1]
    for kind in ("manifest", "analysis", "scan"):
        assert json.loads((root / f"src/editwitness/schemas/{kind}.schema.json").read_text()) == schema_for(kind)
