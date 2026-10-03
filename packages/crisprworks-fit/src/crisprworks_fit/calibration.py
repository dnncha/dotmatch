"""Gene-level negative-control calibration of frozen fitted effects.

Unlike guide shuffling, this null retains within-gene guide correlation.
Validity requires independently chosen controls whose null scores are
exchangeable with tested null genes within each guide-count stratum.
"""

import csv
import json
import sys
from .publication import unique_staging_path, reject_input_collisions, publish_table_bundle, verify_output_record, read_text_input
from pathlib import Path

import numpy as np

from .inference import finite_permutation_tails


def adjust_fdr(pvalues, method="by"):
    """BH or BY adjusted p-values, preserving every declared hypothesis.

    BY handles arbitrary dependence when marginal p-values are valid. BH
    requires independence or suitable positive dependence. Neither repairs
    nonrepresentative controls. Invalid p-values fail rather than disappear.
    """
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1 or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError("FDR input must be a finite 1D array of p-values in [0, 1]")
    if method not in ("bh", "by"):
        raise ValueError("FDR method must be bh or by")
    if not p.size:
        return p.copy()
    order = np.argsort(p, kind="stable")
    factor = np.sum(1 / np.arange(1, p.size + 1)) if method == "by" else 1.
    sorted_q = np.minimum.accumulate((p[order] * p.size * factor / np.arange(1, p.size + 1))[::-1])[::-1]
    q = np.empty_like(p)
    q[order] = np.minimum(1., sorted_q)
    return q


def calibrate_genes(genes, control_genes, *, score="effect", method="by", family="global",
                    min_controls=20, max_guides=40, alpha=.05):
    """Calibrate full-precision Fit records without refitting their effects.

    Training controls are excluded from the declared testing family and have
    p=q=1 in outputs. Each guide-count group uses only its own
    controls. Undersupported strata, skipped fits and nonfinite scores receive
    p=1. The default adjusts all test-gene-by-condition hypotheses together, separately
    for each of the three predeclared alternatives (lower, upper, two-sided).
    Selecting a direction after seeing results requires the two-sided family.
    """
    if score not in ("effect", "beta", "z") or method not in ("bh", "by") or family not in ("global", "condition"):
        raise ValueError("Unsupported score, FDR method or testing family")
    if not isinstance(min_controls, int) or min_controls < 2:
        raise ValueError("min_controls must be an integer >=2")
    if not isinstance(max_guides, int) or max_guides < 2:
        raise ValueError("max_guides must be an integer >=2")
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) or not np.isfinite(alpha) or not 0 < alpha <= 1:
        raise ValueError("alpha must be finite and in (0, 1]")
    if not genes:
        raise ValueError("At least one gene is required")
    names = sorted(genes)
    controls = set(control_genes)
    if not controls:
        raise ValueError("At least one control gene is required")
    missing = controls - set(names)
    if missing:
        raise ValueError("Control genes absent from fit: " + ", ".join(sorted(missing)[:5]))
    guide_counts, rows = [], []
    if len({"effect_estimate" in gene for gene in genes.values()}) != 1:
        raise ValueError("Genes mix incompatible effect models")
    for name in names:
        gene = genes[name]
        if "fit_diagnostics" in gene:
            diagnostics = gene["fit_diagnostics"]
            if not isinstance(diagnostics, dict) or not isinstance(diagnostics.get("converged"), bool):
                raise ValueError(f"Malformed fit diagnostics for {name}")
            if not diagnostics["converged"]:
                raise ValueError(f"Nonconverged fit for {name}; refit before calibration")
            iterations = diagnostics.get("iterations")
            change = diagnostics.get("final_bound_change")
            if (not isinstance(iterations, int) or isinstance(iterations, bool) or iterations < 1
                    or diagnostics.get("termination_reason") != "bound_tolerance"
                    or not isinstance(change, (int, float)) or isinstance(change, bool)
                    or not np.isfinite(change) or change < 0):
                raise ValueError(f"Malformed fit diagnostics for {name}")
        n = gene["guides"]
        if not isinstance(n, int) or isinstance(n, bool) or n < 1:
            raise ValueError(f"Invalid guide count for {name}")
        joint = "effect_estimate" in gene
        beta = np.asarray(gene["effect_estimate" if joint else "beta_estimate"], dtype=float)
        z = np.asarray(gene["effect_zscore" if joint else "beta_zscore"], dtype=float)
        if beta.ndim != 1 or z.ndim != 1 or not z.size or beta.size != (z.size if joint else n + z.size):
            raise ValueError(f"Malformed effect dimensions for {name}")
        rows.append((beta if joint else beta[n:]) if score in ("effect", "beta") else z)
        guide_counts.append(n)
    if len({row.size for row in rows}) != 1:
        raise ValueError("Genes have different condition counts")
    scores = np.stack(rows)
    sizes = np.asarray(guide_counts)
    training = np.asarray([name in controls for name in names])
    p_two, p_upper, p_lower = (np.ones_like(scores) for _ in range(3))
    statuses = np.full(scores.shape, "tested", dtype=object)
    groups = []
    for n in sorted(set(guide_counts)):
        members = sizes == n
        null = scores[members & training]
        group = {"guides": int(n), "genes": int(members.sum()), "controls": len(null),
                 "minimum_directional_p": 1. / (len(null) + 1),
                 "nonfinite_null_by_condition": (~np.isfinite(null)).sum(axis=0).tolist()}
        groups.append(group)
        if not joint and n >= max_guides:
            statuses[members] = "skipped_fit"
        elif len(null) < min_controls:
            statuses[members] = "insufficient_controls"
        else:
            tails = finite_permutation_tails(null, scores[members])
            for p, values in zip((p_two, p_upper, p_lower), tails):
                p[members] = values
            for condition in range(scores.shape[1]):
                if not np.isfinite(null[:, condition]).all():
                    statuses[members, condition] = "invalid_control_fit"
                invalid = members & ~np.isfinite(scores[:, condition])
                statuses[invalid, condition] = "invalid_fit"
    for p in (p_two, p_upper, p_lower):
        p[training] = 1.
    statuses[training] = "training_control"
    adjusted = []
    testing = ~training
    # Optimistic simultaneous lower bounds preserve the same family, including
    # unsupported/invalid hypotheses at p=1. They do not predict realized power.
    floor = np.ones_like(scores)
    for group in groups:
        eligible = (sizes == group["guides"])[:, None] & (statuses == "tested")
        floor[eligible] = group["minimum_directional_p"]
    adjusted_floor = []
    for index, p in enumerate((p_two, p_upper, p_lower, np.minimum(1., 2 * floor), floor)):
        q = np.ones_like(p)
        if family == "global":
            q[testing] = adjust_fdr(p[testing].ravel(), method).reshape(p[testing].shape)
        else:
            for column in range(p.shape[1]):
                q[testing, column] = adjust_fdr(p[testing, column], method)
        (adjusted if index < 3 else adjusted_floor).append(q)
    directional, two_sided = adjusted_floor[1][testing], adjusted_floor[0][testing]
    resolution = {
        "alpha": float(alpha),
        "minimum_attainable_directional_q": float(directional.min()) if directional.size else None,
        "minimum_attainable_two_sided_q": float(two_sided.min()) if two_sided.size else None,
        "directional_hypotheses_with_attainable_q": int((directional <= alpha).sum()),
        "two_sided_hypotheses_with_attainable_q": int((two_sided <= alpha).sum()),
        "scope": "optimistic lower bounds for the declared family and current eligible strata; not expected power or a biological error-control guarantee",
    }
    result = []
    for i, name in enumerate(names):
        for condition in range(scores.shape[1]):
            result.append({"gene": name, "condition": condition + 1, "guides": guide_counts[i],
                           "score": float(scores[i, condition]) if np.isfinite(scores[i, condition]) else None,
                           "status": statuses[i, condition], "p_two": float(p_two[i, condition]),
                           "p_positive": float(p_upper[i, condition]), "p_negative": float(p_lower[i, condition]),
                           "q_two": float(adjusted[0][i, condition]), "q_positive": float(adjusted[1][i, condition]),
                           "q_negative": float(adjusted[2][i, condition])})
    return {"schema_version": 1, "method": "guide-count-stratified gene-control empirical tails",
            "score": score, "fdr_method": method, "testing_family": family,
            "minimum_control_genes": min_controls, "maximum_guides_exclusive": None if joint else max_guides,
            "family_size": int(testing.sum() * scores.shape[1] if family == "global" else testing.sum()),
            "training_controls_excluded_from_family": int(training.sum()),
            "genes": len(names), "conditions": scores.shape[1], "groups": groups, "rows": result,
            "resolution": resolution,
            "scope": "validity requires independently chosen control genes exchangeable with null test genes within strata after normalization and fitting; not guaranteed biological FDR"}


def run_calibration(args):
    text, details_record = read_text_input(args.fit_details)
    control_text, control_record = read_text_input(args.control_gene)
    inputs = {"fit_details": details_record, "control_gene": control_record}
    data = json.loads(text)
    if not isinstance(data, dict) or data.get("schema_version") != 1 or not isinstance(data.get("genes"), dict):
        raise ValueError("Expected schema-version-1 Fit details")
    verify_output_record(data)
    controls = [name.strip() for name in control_text.splitlines() if name.strip()]
    result = calibrate_genes(data["genes"], controls, score=args.score, method=args.adjust,
                             family=args.family, min_controls=args.min_controls, max_guides=args.max_guides,
                             alpha=getattr(args, "fdr_alpha", .05))
    if any("fit_diagnostics" in gene for gene in data["genes"].values()) and (
            "max_iterations" in data or "absolute_bound_tolerance" in data):
        limit, tolerance = data.get("max_iterations"), data.get("absolute_bound_tolerance")
        if (not isinstance(limit, int) or isinstance(limit, bool) or limit < 1
                or not isinstance(tolerance, (int, float)) or isinstance(tolerance, bool)
                or not np.isfinite(tolerance) or tolerance <= 0):
            raise ValueError("Malformed joint stopping settings")
        for name, gene in data["genes"].items():
            diagnostics = gene.get("fit_diagnostics")
            if diagnostics is None or diagnostics["iterations"] > limit or diagnostics["final_bound_change"] >= tolerance:
                raise ValueError(f"Fit diagnostics contradict stopping settings for {name}")
    if "conditions" in data:
        labels = data["conditions"]
        if (not isinstance(labels, list) or len(labels) != result["conditions"] or
                any(not isinstance(label, str) or not label.strip() for label in labels) or
                len(set(labels)) != len(labels)):
            raise ValueError("Condition labels must be distinct nonempty strings matching fitted effects")
        result["condition_labels"] = labels
        for row in result["rows"]:
            row["condition_label"] = labels[row["condition"] - 1]
    joint = "effect_estimate" in next(iter(data["genes"].values()))
    result["effect_model"] = data.get("model", "JACKS-derived joint guide-efficacy inference" if joint else "MAGeCK2 MLE")
    result["effect_units"] = data.get("effect_units", "log2 relative abundance" if joint else "natural-log beta coefficient")
    from . import __version__
    result["versions"] = {"crisprworks_fit": __version__, "numpy": np.__version__}
    if result["family_size"]:
        for alternative in ("directional", "two_sided"):
            if not result["resolution"][f"{alternative}_hypotheses_with_attainable_q"]:
                print(f"crisprworks-fit: warning: current control resolution and eligible strata prevent any "
                      f"{alternative} discovery at q <= {result['resolution']['alpha']:g}; "
                      "no calls would not establish absence of effects", file=sys.stderr)
    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    output = Path(str(prefix) + ".calibration.tsv")
    manifest = Path(str(prefix) + ".calibration.json")
    companion = data.get("output", {}).get("path")
    reject_input_collisions((output, manifest), (args.fit_details, args.control_gene, companion))
    # Stage both files before publication; a manifest identifies completion.
    temporary = unique_staging_path(output)
    try:
        with temporary.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(result["rows"][0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(result.pop("rows"))
        result["inputs"] = inputs
        if companion is not None:
            result["inputs"]["companion_output"] = data["output"].copy()
        publish_table_bundle(temporary, output, manifest, result)
    finally:
        temporary.unlink(missing_ok=True)
    return output
