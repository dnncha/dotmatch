"""Offline synthetic example of joint inference and gene-control calibration."""

import csv
import json
from pathlib import Path
from types import SimpleNamespace

from . import __version__


def run_joint_demo(directory):
    """Publish a complete example only in a new directory.

    The generator declares neutral genes before fitting. These simulated
    controls are solely for exercising the software, not biological evidence.
    """
    from .calibration import run_calibration
    from .cli import input_record, synthetic_screen
    from .joint import run_joint

    directory = Path(directory)
    if directory.exists():
        raise ValueError("Joint demo output already exists; choose a new directory")
    directory.mkdir(parents=True, exist_ok=False)
    counts, design = synthetic_screen(directory, genes=48, guides=4, seed=1729)
    mapping = directory / "samples.tsv"
    with mapping.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["Sample", "Condition", "Control"])
        for condition, samples in (("BASE", ("baseline1", "baseline2")),
                                   ("CONTROL", ("control1", "control2")),
                                   ("TREATED", ("treated1", "treated2"))):
            writer.writerows((sample, condition, "BASE") for sample in samples)
    controls = directory / "controls.txt"
    controls.write_text("".join(f"GENE{gene}\n" for gene in range(48) if gene % 6 >= 2))
    prefix = directory / "screen"
    summary = run_joint(SimpleNamespace(counts=counts, sample_map=mapping,
                                        control_gene=controls, output_prefix=prefix))
    details = Path(str(prefix) + ".joint-details.json")
    run_calibration(SimpleNamespace(fit_details=details, control_gene=controls, output_prefix=prefix,
                                    score="effect", adjust="by", family="global", min_controls=20, max_guides=40))
    manifest = {"schema_version": 1, "status": "complete", "tool": "CRISPRWorks Fit joint demo",
                "version": __version__, "seed": 1729, "genes": 48, "guides_per_gene": 4,
                "control_genes": 32, "conditions": ["CONTROL", "TREATED"],
                "scope": "Synthetic software demonstration; controls declared by the generator before fitting; not biological validation",
                "inputs": {name: input_record(str(path)) for name, path in
                           (("counts", counts), ("sample_map", mapping), ("control_genes", controls), ("mle_design", design))},
                "outputs": {name: input_record(str(path)) for name, path in
                            (("joint_summary", summary), ("joint_details", details),
                             ("calibration_table", Path(str(prefix) + ".calibration.tsv")),
                             ("calibration_manifest", Path(str(prefix) + ".calibration.json")))}}
    (directory / "demo.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return summary
