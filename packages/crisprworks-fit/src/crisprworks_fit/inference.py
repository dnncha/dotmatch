"""Conservative finite Monte Carlo tails for the experimental MLE workflow.

This changes inference, not the fitted effects or the pooled-guide null model.
The add-one correction requires exchangeability for a randomization-test
interpretation; MAGeCK's pooled pseudo-genes do not establish that assumption.
"""

import numpy as np


def finite_permutation_tails(null_betas, observed_betas):
    """Return two-sided, upper and lower add-one tails, counting ties.

    For B null draws, p = (1 + number at least as extreme) / (B + 1).
    Two-sided p-values double the smaller directional tail and are capped at
    one. A nonfinite observation, or any nonfinite null in its column, receives
    p=1 in all tails: a failed fit must never become a discovery. Null draws
    remain in the denominator; we do not silently discard failed draws.
    """
    null = np.asarray(null_betas, dtype=float)
    observed = np.asarray(observed_betas, dtype=float)
    if null.ndim != 2 or observed.ndim != 2 or null.shape[1] != observed.shape[1]:
        raise ValueError("Null and observed betas must be 2D with matching columns")
    if null.shape[0] == 0:
        raise ValueError("At least one permutation beta is required")
    upper = np.ones(observed.shape, dtype=float)
    lower = np.ones(observed.shape, dtype=float)
    denominator = null.shape[0] + 1
    for column in range(null.shape[1]):
        if not np.isfinite(null[:, column]).all():
            continue
        values = np.sort(null[:, column])
        valid = np.isfinite(observed[:, column])
        scores = observed[valid, column]
        lower[valid, column] = (1 + np.searchsorted(values, scores, side="right")) / denominator
        upper[valid, column] = (1 + len(values) - np.searchsorted(values, scores, side="left")) / denominator
    return np.minimum(1.0, 2 * np.minimum(upper, lower)), upper, lower


def assign_finite_permutation_pvalues(null_betas, genedict, *, max_guides=None, diagnostics=None):
    genes = list(genedict.values())
    if not genes:
        return
    observed = np.asarray([gene.beta_estimate[gene.nb_count.shape[1]:] for gene in genes])
    two_sided, upper, lower = finite_permutation_tails(null_betas, observed)
    skipped = np.array([
        max_guides is not None and gene.nb_count.shape[1] >= max_guides for gene in genes
    ])
    for values in (two_sided, upper, lower):
        values[skipped] = 1.0
    for i, gene in enumerate(genes):
        gene.beta_permute_pval = two_sided[i]
        gene.beta_permute_pval_pos = upper[i]
        gene.beta_permute_pval_neg = lower[i]
    if diagnostics is not None:
        null = np.asarray(null_betas)
        diagnostics.append({
            "genes": len(genes), "guide_counts": sorted({gene.nb_count.shape[1] for gene in genes}),
            "null_draws": null.shape[0], "conditions": null.shape[1],
            "minimum_directional_p": 1.0 / (null.shape[0] + 1),
            "minimum_two_sided_p": min(1.0, 2.0 / (null.shape[0] + 1)),
            "nonfinite_null_draws_by_condition": (~np.isfinite(null)).sum(axis=0).tolist(),
            "nonfinite_observed_by_condition": (~np.isfinite(observed)).sum(axis=0).tolist(),
            "skipped_genes": int(skipped.sum()),
        })
