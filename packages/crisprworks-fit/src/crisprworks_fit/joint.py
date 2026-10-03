"""Joint-screen inference with a JACKS-derived guide-efficacy model.

Posterior variance smoothing is adapted from JACKS preprocessing by Felicity
Allen and Leopold Parts (2018, MIT; see licenses/JACKS-MIT.txt), commit
dd5c4be5e83baa5ee79a589b0a8c3a2ac3f7d6ad. Changes: strict tabular validation,
replicated treatment requirement, and omission of subsampling/interpolation.
The smoothing and default variational updates match the upstream algorithm.
"""

import csv
import hashlib
import json
import io
from .publication import unique_staging_path, reject_input_collisions, publish_table_bundle, read_text_input
from pathlib import Path
import platform
import time

import numpy as np

from ._joint_kernel import inferJACKSGene

JACKS_COMMIT = "dd5c4be5e83baa5ee79a589b0a8c3a2ac3f7d6ad"


def posterior_sd(data):
    """Upstream mean/variance window smoothing, with single-baseline NaNs."""
    n = data.shape[0]
    if data.shape[1] == 1:
        return np.full(n, np.nan)
    window = min(max(int(n / 100.), 30), 800)
    ordered = np.array(sorted(zip(data.mean(axis=1), np.nanstd(data, axis=1)**2,
                                  range(n))))
    values = ordered[:, 1]
    smooth = np.zeros(n)
    for k in range(window + 1):
        smooth[k] = np.mean(values[:2 * window + 1])
    for k in range(window + 1, n - window - 1):
        smooth[k] = smooth[k - 1] + (values[k + window + 1] - values[k - window]) / (2 * window + 1)
    for k in range(n - window - 1, n):
        smooth[k] = np.mean(values[n - 2 * window - 1:])
    for i in range(1, n):
        if not np.isnan(smooth[n - i - 1]):
            smooth[n - i - 1] = np.nanmax(smooth[n - i - 1:n - i + 1])
    result = np.zeros(n)
    result[ordered[:, 2].astype(int)] = smooth**.5
    return result


def read_screen(count_file, sample_map, input_records=None):
    """Read a guide-by-sample count TSV and Sample/Condition/Control TSV."""
    count_text, count_record = read_text_input(count_file)
    with io.StringIO(count_text, newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader, [])
        if len(header) < 3 or len(set(header)) != len(header) or any(not s for s in header):
            raise ValueError("Count table requires distinct guide, gene and sample columns")
        rows = list(reader)
    if not rows or any(len(row) != len(header) for row in rows):
        raise ValueError("Count table has missing or inconsistent rows")
    guides, genes = [r[0] for r in rows], [r[1] for r in rows]
    if len(set(guides)) != len(guides) or any(not s.strip() for s in guides + genes):
        raise ValueError("Guide IDs must be unique and guide/gene IDs nonempty")
    counts = np.asarray([r[2:] for r in rows], dtype=float)
    if not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError("Counts must be finite and nonnegative")
    if len(guides) < 64:
        raise ValueError("Joint variance smoothing requires at least 64 guides")
    map_text, map_record = read_text_input(sample_map)
    with io.StringIO(map_text, newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["Sample", "Condition", "Control"]:
            raise ValueError("Sample map header must be Sample, Condition, Control (tab-separated)")
        metadata = list(reader)
    if input_records is not None:
        input_records.update(counts=count_record, sample_map=map_record)
    if not metadata or any(set(r) != {"Sample", "Condition", "Control"} or
                           any(not isinstance(v, str) or not v.strip() for v in r.values())
                           for r in metadata):
        raise ValueError("Sample map contains missing or inconsistent values")
    samples = [r["Sample"] for r in metadata]
    if len(set(samples)) != len(samples) or set(samples) != set(header[2:]):
        raise ValueError("Sample map must map each count sample exactly once")
    lookup = {r["Sample"]: r for r in metadata}
    control_by_condition = {}
    for row in metadata:
        condition, control = row["Condition"], row["Control"]
        if condition in control_by_condition and control_by_condition[condition] != control:
            raise ValueError("Replicates of a condition must use the same baseline condition")
        control_by_condition[condition] = control
    conditions = sorted(control_by_condition)
    for condition, control in control_by_condition.items():
        if control not in control_by_condition or control_by_condition[control] != control:
            raise ValueError("Each control must name a self-mapped baseline condition")
        if condition != control and sum(r["Condition"] == condition for r in metadata) < 2:
            raise ValueError("Each treatment condition requires at least two replicates")
    treatments = [s for s in conditions if control_by_condition[s] != s]
    if not treatments:
        raise ValueError("At least one treatment condition is required")
    columns = [[i for i, sample in enumerate(header[2:]) if lookup[sample]["Condition"] == condition]
               for condition in conditions]
    order = np.argsort(np.asarray(guides))
    return (counts[order], np.asarray(guides)[order], np.asarray(genes)[order],
            conditions, columns, control_by_condition, treatments)


def fit_joint(count_file, sample_map, control_genes=None, *, max_iterations=50, tolerance=.1, input_records=None):
    """Fit effects in log2 units with guide efficacies shared across conditions.

    Efficacies are continuous relative weights centered near one, not editing
    probabilities. Posterior z scores are descriptive; calibrated gene-control
    tails must be requested separately. No Gaussian posterior FDR is implied.
    """
    if not isinstance(max_iterations, int) or isinstance(max_iterations, bool) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer")
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive")
    counts, guides, genes, conditions, columns, controls, treatments = read_screen(count_file, sample_map, input_records)
    gene_names = sorted(set(genes))
    training = set(control_genes or [])
    if control_genes is not None and not training:
        raise ValueError("Negative-control gene file is empty")
    missing = training - set(gene_names)
    if missing:
        raise ValueError("Normalization control genes absent from counts: " + ", ".join(sorted(missing)[:5]))
    logcounts = np.log2(counts + 32.)
    norm_rows = np.isin(genes, list(training)) if training else np.ones(len(genes), dtype=bool)
    logcounts -= np.tile(np.nanmedian(logcounts[norm_rows], axis=0), (len(genes), 1))
    condensed = np.zeros((len(genes), len(conditions), 2))
    for j, indexes in enumerate(columns):
        values = logcounts[:, indexes]
        condensed[:, j, 0] = values.mean(axis=1)
        condensed[:, j, 1] = posterior_sd(values)
    treatment_indexes = [conditions.index(s) for s in treatments]
    baseline_indexes = [conditions.index(controls[s]) for s in treatments]
    test, baseline = condensed[:, treatment_indexes], condensed[:, baseline_indexes]
    results = {}
    for gene in gene_names:
        indexes = np.flatnonzero(genes == gene)
        diagnostics = {}
        _, _, x, x2, w, w2 = inferJACKSGene(test[indexes, :, 0], test[indexes, :, 1].copy(),
                                         baseline[indexes, :, 0], baseline[indexes, :, 1].copy(),
                                         max_iterations, tol=tolerance, diagnostics=diagnostics)
        effect_var, efficacy_var = w2 - w**2, x2 - x**2
        if (not all(np.isfinite(a).all() for a in (x, x2, w, w2)) or
                (effect_var <= 0).any() or (efficacy_var < 0).any()):
            raise ValueError(f"Joint inference failed for {gene}; no results published")
        results[gene] = {"guides": len(indexes), "guide_ids": guides[indexes].tolist(),
                         "fit_diagnostics": diagnostics,
                         "effect_estimate": w.tolist(), "effect_sd": np.sqrt(effect_var).tolist(),
                         "effect_zscore": (w / np.sqrt(effect_var)).tolist(),
                         "guide_efficacy": x.tolist(), "guide_efficacy_sd": np.sqrt(efficacy_var).tolist()}
    return {"schema_version": 1, "model": "JACKS-derived joint guide-efficacy inference",
            "source_commit": JACKS_COMMIT, "effect_units": "log2 relative abundance",
            "conditions": treatments, "baseline_by_condition": {s: controls[s] for s in treatments},
            "normalization": "control-guide median log counts" if training else "all-guide median log counts",
            "training_controls": sorted(training), "count_pseudocount": 32,
            "max_iterations": max_iterations, "absolute_bound_tolerance": tolerance,
            "nonconverged_genes": sum(not result["fit_diagnostics"]["converged"] for result in results.values()),
            "genes": results}


def run_joint(args):
    started = time.perf_counter()
    inputs = {}
    controls = None
    if args.control_gene is not None:
        text, inputs["control_gene"] = read_text_input(args.control_gene)
        controls = [s.strip() for s in text.splitlines() if s.strip()]
    details = fit_joint(args.counts, args.sample_map, controls,
                        max_iterations=getattr(args, "max_iterations", 50),
                        tolerance=getattr(args, "tolerance", .1), input_records=inputs)
    details["inputs"] = inputs
    details["elapsed_seconds"] = time.perf_counter() - started
    from . import __version__
    details["versions"] = {"crisprworks_fit": __version__, "python": platform.python_version(), "numpy": np.__version__}
    details["implementation_sha256"] = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                         for name in ("joint.py", "_joint_kernel.py")}
    output = Path(str(args.output_prefix) + ".joint-details.json")
    summary = Path(str(args.output_prefix) + ".joint.tsv")
    reject_input_collisions((output, summary), (args.counts, args.sample_map, args.control_gene))
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = unique_staging_path(summary)
    try:
        with temporary.open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(["gene", "condition", "guides", "effect", "posterior_sd", "posterior_z"])
            for gene, values in details["genes"].items():
                for j, condition in enumerate(details["conditions"]):
                    writer.writerow([gene, condition, values["guides"], values["effect_estimate"][j],
                                     values["effect_sd"][j], values["effect_zscore"][j]])
        publish_table_bundle(temporary, summary, output, details)
    finally:
        temporary.unlink(missing_ok=True)
    return summary
