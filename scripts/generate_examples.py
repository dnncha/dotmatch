"""Regenerate synthetic examples deterministically; never fetch biological data."""
import json
import random
from pathlib import Path

from editwitness.models import Manifest

root = Path(__file__).resolve().parents[1]
rng = random.Random(20260905)
sequence = "".join(rng.choice("ACGT") for _ in range(900))
alt = next(base for base in "ACGT" if base != sequence[450])
data = {
    "schema_version": "1.1", "observation_model": "exact-local-sequence-presence-v2", "coordinate_system": "0-based-half-open",
    "reference": {"name": "Synthetic 900-bp teaching locus (not biological data)", "sequence": sequence, "synthetic": True},
    "alleles": [
        {"id": "reference", "edits": []},
        {"id": "intended", "edits": [{"start": 450, "end": 451, "sequence": alt}]},
        {"id": "primer_site_deletion", "description": "Hypothetical deletion removes the inner left-primer site",
         "edits": [{"start": 180, "end": 480, "sequence": ""}]},
        {"id": "interior_deletion", "description": "Hypothetical deletion within the paired-end read gap",
         "edits": [{"start": 400, "end": 500, "sequence": ""}]},
        {"id": "window_deleted", "description": "Hypothetical deletion of the entire supplied reference window, not a claim about a whole chromosome",
         "edits": [{"start": 0, "end": 900, "sequence": ""}]},
    ],
    "hypotheses": [
        {"id": "intended_biallelic", "alleles": ["intended", "intended"]},
        {"id": "intended_reference", "alleles": ["intended", "reference"]},
        {"id": "hidden_primer_deletion", "alleles": ["intended", "primer_site_deletion"]},
        {"id": "hidden_window_deletion", "alleles": ["intended", "window_deleted"]},
        {"id": "interior_deletion", "alleles": ["intended", "interior_deletion"]},
    ],
    "expected_hypothesis": "intended_biallelic",
    "assays": [{"id": "inner", "left_primer": {"start": 200, "end": 220},
                "right_primer": {"start": 680, "end": 700}, "readout": "full_insert", "max_product_bp": 1000}],
    "candidates": [
        {"id": "outer", "left_primer": {"start": 50, "end": 70}, "right_primer": {"start": 830, "end": 850},
         "readout": "full_insert", "max_product_bp": 1500, "cost_units": 2},
        {"id": "inner_repeat", "left_primer": {"start": 200, "end": 220}, "right_primer": {"start": 680, "end": 700},
         "readout": "full_insert", "max_product_bp": 1000, "cost_units": 1},
    ],
    "deletion_scan": {"start_min": 0, "start_max": 450, "end_min": 451, "end_max": 900, "step": 10},
}
for filename in ("demo.json", "paired_end.json"):
    if filename == "paired_end.json":
        data["assays"][0].update(readout="paired_end", read_bases=80)
    manifest = Manifest.model_validate(data)
    text = json.dumps(manifest.model_dump(mode="json"), indent=2) + "\n"
    for directory in (root / "examples", root / "src/editwitness/data"):
        (directory / filename).write_text(text, encoding="utf-8")
(root / "examples/synthetic.fasta").write_text(
    ">Synthetic teaching locus; not biological data\n" + "\n".join(sequence[i:i+80] for i in range(0, len(sequence), 80)) + "\n",
    encoding="utf-8",
)
