"""MAGeCK2-compatible EM and permutation kernels.

The EM control flow and statistical conventions derive from MAGeCK2,
Copyright (c) 2021 Wei Li, BSD-3-Clause; see LICENSE.
"""

import numpy as np
from scipy.special import gammaln, xlog1py

try:
    from . import _native
except ImportError:
    _native = None

_engine = "auto"


def resolved_engine(requested="auto"):
    if requested not in ("auto", "native", "numpy"):
        raise ValueError("Kernel must be auto, native or numpy")
    if requested == "native" and _native is None:
        raise RuntimeError("Native Fit extension unavailable; build the package or select --kernel numpy")
    return "native" if requested != "numpy" and _native is not None else "numpy"


def em_whileloop(sk, beta_init_mat_0, size_vec, wfrac_list_0,
                 alpha_dispersion, alpha_val, estimateeff, updateeff,
                 removeoutliers, debug):
    if debug or resolved_engine(_engine) == "numpy":
        return em_whileloop_numpy(sk, beta_init_mat_0, size_vec, wfrac_list_0,
                                 alpha_dispersion, alpha_val, estimateeff, updateeff,
                                 removeoutliers, debug)
    from mageck2.mledesignmat import DesignMatCache
    n = sk.nb_count.shape[1]
    other = sk.design_mat.shape[0] - 1
    conditions = sk.design_mat.shape[1] - 1
    design = np.require(np.asarray(DesignMatCache.get_record(n)[2]), dtype=np.float64, requirements=["C", "A"])
    m, p = design.shape
    def buffer(values):
        return np.require(np.asarray(values).reshape(-1), dtype=np.float64, requirements=["C", "A"])
    dispersion = np.require(np.broadcast_to(np.asarray(alpha_dispersion), (m,)), dtype=np.float64, requirements=["C", "A"])
    try:
        result = _native.loop(design, buffer(sk.sgrna_kvalue), buffer(size_vec),
                              buffer(beta_init_mat_0), buffer(wfrac_list_0), dispersion,
                              n, conditions, other, alpha_val, estimateeff, updateeff)
    except ArithmeticError as error:
        raise np.linalg.LinAlgError(str(error)) from error
    efficiency, beta, mean, residual, weights, gram, regularized = (
        np.frombuffer(value, dtype=np.float64) for value in result[1:]
    )
    mean = mean.reshape(m, 1)
    residual = residual.reshape(m, 1)
    covariance = uncertainty(design, weights, gram.reshape(p, p),
                             regularized.reshape(p, p), residual / mean, n)
    return (result[0], efficiency.copy(), np.matrix(beta.reshape(p, 1)),
            np.matrix(covariance), np.matrix(mean), np.matrix(residual))


def negative_binomial_loglikelihood(counts, mean, dispersion):
    """Evaluate the reference NB log PMF without scipy.stats argument parsing.

    Retain upstream's variance arithmetic, rounded counts, distribution support
    and NaN-to-zero convention. In particular, do not simplify r to 1/alpha:
    that changes cancellation behavior for small means and dispersions.
    """
    counts = np.asarray(counts)
    mean = np.asarray(mean)
    if counts.shape[0] != mean.shape[0]:
        raise ValueError("Count table dimension is not the same as mu vector dimension.")
    with np.errstate(all="ignore"):
        k = np.round(counts)
        mean_squared = mean * mean
        variance = mean + np.asarray(dispersion) * mean_squared
        p = mean / variance
        r = mean_squared / (variance - mean)
        k, r, p = np.broadcast_arrays(k, r, p)
        valid_parameters = (r > 0) & (p > 0) & (p <= 1)
        supported = (k >= 0) & (k <= np.inf) & (k == np.floor(k))
        values = (gammaln(r + k) - gammaln(k + 1) - gammaln(r)
                  + r * np.log(p) + xlog1py(k, -p))
        result = np.where(supported, values, -np.inf)
        result = np.where(valid_parameters & ~np.isnan(k), result, np.nan)
        return np.where(np.isnan(result), 0, result)


def weighted_fit(design, weights, response, ridge):
    """Solve ridge IRLS without constructing a diagonal weight matrix."""
    gram = design.T @ (weights[:, None] * design)
    regularized = gram.copy()
    regularized.flat[:: regularized.shape[0] + 1] += ridge
    beta = np.linalg.solve(regularized, design.T @ (weights[:, None] * response))
    return beta, gram, regularized


def uncertainty(design, weights, gram, regularized, residual_ratio, n_guides):
    """MAGeCK-NEST sandwich covariance and effective residual degrees of freedom.

    If G = X'WX and A = G + ridge*I, the hat matrix has trace tr(A^-1 G)
    and squared Frobenius norm tr((A^-1 G)^2). This avoids its m-by-m storage.
    """
    inverse = np.linalg.solve(regularized, np.eye(regularized.shape[0]))
    influence = inverse @ gram
    degrees = design.shape[0] - (
        2 * np.trace(influence) - np.einsum("ij,ji->", influence, influence)
    )
    # Only the condition block is consumed downstream.
    condition_rows = inverse[n_guides:, :]
    weighted_projection = (condition_rows @ design.T) * weights
    covariance = weighted_projection @ weighted_projection.T
    scale = np.sum(residual_ratio**2) / degrees
    return scale * covariance


def em_whileloop_numpy(sk, beta_init_mat_0, size_vec, wfrac_list_0,
                 alpha_dispersion, alpha_val, estimateeff, updateeff,
                 removeoutliers, debug):
    """Replacement for the pinned MAGeCK2 EM loop, retaining its conventions.

    Upstream's active solution does not use removeoutliers. Return values retain
    np.matrix types at the adapter boundary for its existing downstream code.
    """
    from mageck2.mledesignmat import DesignMatCache

    beta = np.asarray(beta_init_mat_0).copy()
    efficiency = np.asarray(wfrac_list_0).copy()
    counts = np.asarray(sk.sgrna_kvalue)
    n_guides = sk.nb_count.shape[1]
    n_other = sk.design_mat.shape[0] - 1
    n_conditions = sk.design_mat.shape[1] - 1
    design = np.asarray(DesignMatCache.get_record(n_guides)[2])
    sizes = np.asarray(size_vec)
    dispersion = np.asarray(alpha_dispersion)
    likelihood_dispersion = (
        alpha_dispersion if dispersion.ndim == 0 else dispersion.reshape(-1, 1)
    )
    status = 0
    iteration = 1
    while True:
        log_mean = design @ beta
        mean = sizes * np.exp(log_mean)
        residual = counts - mean
        if estimateeff:
            coefficient = np.concatenate((
                np.ones(n_guides), np.tile(efficiency, n_other),
                np.tile(1 - efficiency, n_other),
            )).reshape(-1, 1)
            response = residual * coefficient / mean + log_mean
        else:
            response = residual / mean + log_mean

        # In upstream this likelihood is otherwise a discarded diagnostic.
        if (estimateeff and updateeff) or debug:
            likelihood = negative_binomial_loglikelihood(
                counts, mean, likelihood_dispersion,
            )
        if estimateeff and updateeff:
            baseline = beta.copy()
            baseline[n_guides:n_guides + n_conditions] = 0
            baseline_mean = sizes * np.exp(design @ baseline)
            baseline_likelihood = negative_binomial_loglikelihood(
                counts, baseline_mean, likelihood_dispersion,
            )
            difference = np.asarray(baseline_likelihood - likelihood)
            difference = difference[n_guides:n_guides * (n_other + 1)]
            difference = difference.reshape(n_other, n_guides).T
            difference = np.minimum(difference, 100)
            efficiency = (1 / (1 + np.exp(np.min(difference, axis=1))))
            if debug:
                print("frac:" + " ".join(format(x, ".3g") for x in efficiency))

        weights = 1 / (1 / mean.ravel() + dispersion)
        if estimateeff:
            extended_efficiency = np.concatenate((
                np.tile(efficiency, n_other + 1),
                np.tile(1 - efficiency, n_other),
            ))
            extended_efficiency = np.maximum(extended_efficiency, 1e-2)
            weights[n_guides:] *= extended_efficiency[n_guides:]
        weights = np.maximum(weights, np.max(weights) / 100)
        if debug:
            print("w:" + " ".join(str(x) for x in weights))
        proposal, gram, regularized = weighted_fit(design, weights, response, alpha_val)
        difference = proposal - beta
        # Match upstream's rollback and maximum-coefficient behavior.
        if np.isnan(np.sum(proposal)) or np.max(np.abs(proposal)) > 100:
            status = 2
            break
        beta = proposal
        iteration += 1
        difference_size = np.sum(difference[n_guides:]**2)
        beta_size = np.sum(proposal[n_guides:]**2)
        if abs(beta_size) < 1e-9:
            beta_size = 1.0
        fraction = difference_size / beta_size
        if debug:
            print(f"Iteration {iteration}, updated beta:" +
                  " ".join(format(x, ".3g") for x in proposal.ravel()))
            print(f"log likelihood:{np.sum(likelihood)}, frac:{fraction:.3g},abs:{beta_size:.3g}")
        if fraction < 1e-9:
            break
        if iteration > 1000:
            status = 1
            break

    covariance = uncertainty(
        design, weights, gram, regularized, residual / mean, n_guides,
    )
    return (status, efficiency, np.matrix(beta), np.matrix(covariance),
            np.matrix(mean), np.matrix(residual))


def permutation_tails(null_betas, observed_betas):
    """Strict empirical tails, matching upstream's > and < (including ties).

    NaNs compare false upstream: retain the original denominator while removing
    them from each sorted column. A NaN observation has both tails zero.
    """
    null = np.asarray(null_betas)
    observed = np.asarray(observed_betas)
    if null.ndim != 2 or observed.ndim != 2 or null.shape[1] != observed.shape[1]:
        raise ValueError("Null and observed betas must be 2D with matching columns")
    if null.shape[0] == 0:
        raise ValueError("At least one permutation beta is required")
    upper = np.empty(observed.shape, dtype=float)
    lower = np.empty(observed.shape, dtype=float)
    for column in range(null.shape[1]):
        values = null[:, column]
        values = np.sort(values[~np.isnan(values)])
        scores = observed[:, column]
        lower[:, column] = np.searchsorted(values, scores, side="left") / null.shape[0]
        upper[:, column] = (len(values) - np.searchsorted(values, scores, side="right")) / null.shape[0]
        missing = np.isnan(scores)
        lower[missing, column] = 0
        upper[missing, column] = 0
    return 2 * np.minimum(upper, lower), upper, lower


def assign_p_value_from_permuted_beta(betazeros, genedict):
    genes = list(genedict.values())
    if not genes:
        return
    observed = np.asarray([gene.beta_estimate[gene.nb_count.shape[1]:] for gene in genes])
    two_sided, upper, lower = permutation_tails(betazeros, observed)
    for i, gene in enumerate(genes):
        gene.beta_permute_pval = two_sided[i]
        gene.beta_permute_pval_pos = upper[i]
        gene.beta_permute_pval_neg = lower[i]
