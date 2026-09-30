"""Scientific regression checks against pinned, unmodified MAGeCK2 kernels."""

import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
import warnings

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from threadpoolctl import threadpool_limits

from mageck2 import mleem, mlemultiprocessing
from mageck2.mleclassdef import SimCaseSimple, gene_fdr_correction
from mageck2.mledesignmat import DesignMatCache
from crisprworks_fit import kernels
from crisprworks_fit.backend import accelerated, verify_upstream, Runner
from crisprworks_fit.cli import synthetic_screen

warnings.filterwarnings("ignore", category=PendingDeprecationWarning)


class VarianceModel:
    def get_lm_var(self, values, returnalpha=False):
        return 0.01 + 0.02 / np.sqrt(np.asarray(values) + 1)


def gene_case(seed=1, guides=6, design=None):
    if design is None:
        design = [[1, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 0], [1, 0, 1], [1, 0, 1]]
    rng = np.random.default_rng(seed)
    gene = SimCaseSimple()
    gene.prefix = f"GENE{seed}"
    gene.design_mat = np.matrix(design, dtype=float)
    beta = rng.uniform(-1.3, 0.8, gene.design_mat.shape[1] - 1)
    means = rng.uniform(4, 900, guides)
    sample_effects = np.exp(np.asarray(gene.design_mat)[:, 1:] @ beta)
    expected = sample_effects[:, None] * means[None, :]
    gene.nb_count = np.matrix(rng.negative_binomial(20, 20 / (20 + expected)) + 1.0)
    gene.w_estimate = rng.uniform(0.1, 0.9, guides)
    gene.beta_estimate = []
    gene.sgrnaid = [f"{gene.prefix}_{i}" for i in range(guides)]
    return gene


class LinearAlgebraTests(unittest.TestCase):
    def test_covariance_and_degrees_match_dense_hat_matrix(self):
        rng = np.random.default_rng(71)
        for rows, columns in ((12, 3), (66, 8), (400, 13)):
            x = rng.normal(size=(rows, columns))
            w = rng.uniform(0.01, 100, rows)
            z = rng.normal(size=(rows, 1))
            ratio = rng.normal(size=(rows, 1))
            for ridge in (0.0, 1e-6, 0.01, 3.0):
                with self.subTest(rows=rows, ridge=ridge):
                    diagonal = np.diag(w)
                    inverse = np.linalg.inv(x.T @ diagonal @ x + ridge * np.eye(columns))
                    reference_beta = inverse @ x.T @ diagonal @ z
                    reference_cov = inverse @ x.T @ diagonal @ diagonal @ x @ inverse
                    hat = np.sqrt(diagonal) @ x @ inverse @ x.T @ np.sqrt(diagonal)
                    degrees = rows - np.trace(2 * hat - hat @ hat.T)
                    reference_cov *= np.sum(ratio**2) / degrees
                    beta, gram, regularized = kernels.weighted_fit(x, w, z, ridge)
                    covariance = kernels.uncertainty(x, w, gram, regularized, ratio, 2)
                    assert_allclose(beta, reference_beta, rtol=1e-10, atol=1e-12)
                    assert_allclose(covariance, reference_cov[2:, 2:], rtol=1e-10, atol=1e-12)

    def test_singular_solve_does_not_silently_change_model(self):
        with self.assertRaises(np.linalg.LinAlgError):
            kernels.weighted_fit(np.ones((6, 2)), np.ones(6), np.ones((6, 1)), 0)


class InferenceTests(unittest.TestCase):
    def test_fits_match_reference_across_designs_and_efficiency_modes(self):
        designs = (
            [[1, 0], [1, 0], [1, 1], [1, 1]],
            [[1, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 0], [1, 0, 1], [1, 0, 1]],
            [[1, 0, 0], [1, 0, 0], [1, 1, 0.2], [1, 1, 0.8], [1, 0, 0.5], [1, 0, 1.0]],
        )
        with threadpool_limits(limits=1):
            for design in designs:
                for guides in (1, 4, 12):
                    for seed in (3, 11):
                        for estimate, update in ((False, False), (True, False), (True, True)):
                            with self.subTest(design=design, guides=guides, seed=seed,
                                              estimate=estimate, update=update):
                                original = gene_case(seed, guides, design)
                                reference = copy.deepcopy(original)
                                candidate = copy.deepcopy(original)
                                settings = dict(debug=False, estimateeff=estimate, updateeff=update,
                                                meanvarmodel=VarianceModel(),
                                                size_factor=[0.8 + i * 0.1 for i in range(len(design))],
                                                removeoutliers=True, logem=False)
                                DesignMatCache.cache = {}
                                mleem.iteratenbem(reference, **settings)
                                with accelerated():
                                    mleem.iteratenbem(candidate, **settings)
                                for field in ("beta_estimate", "beta_zscore", "beta_pval",
                                              "beta_pval_pos", "beta_pval_neg", "w_estimate",
                                              "mu_estimate", "sgrna_residule"):
                                    assert_allclose(getattr(candidate, field), getattr(reference, field),
                                                    rtol=1e-7, atol=1e-8, err_msg=field)
                                # Compare standard errors derived from beta/z separately.
                                count = guides
                                reference_se = reference.beta_estimate[count:] / reference.beta_zscore
                                candidate_se = candidate.beta_estimate[count:] / candidate.beta_zscore
                                assert_allclose(candidate_se, reference_se, rtol=1e-7, atol=1e-8)

    def test_restart_matches(self):
        reference = gene_case()
        DesignMatCache.cache = {}
        mleem.iteratenbem(reference, debug=False, estimateeff=True)
        candidate = copy.deepcopy(reference)
        mleem.iteratenbem(reference, debug=False, estimateeff=True, restart=True)
        with accelerated():
            mleem.iteratenbem(candidate, debug=False, estimateeff=True, restart=True)
        assert_allclose(candidate.beta_estimate, reference.beta_estimate, rtol=1e-7, atol=1e-8)
        assert_allclose(candidate.beta_pval, reference.beta_pval, rtol=1e-7, atol=1e-8)

    def test_nested_context_and_failure_restore_patches_and_cache(self):
        original = mleem.em_whileloop
        cache = DesignMatCache.cache
        with self.assertRaisesRegex(RuntimeError, "Nested"):
            with accelerated():
                with accelerated():
                    pass
        self.assertIs(mleem.em_whileloop, original)
        self.assertIs(DesignMatCache.cache, cache)
        verify_upstream()

    def test_untested_upstream_is_rejected(self):
        original = mleem.em_whileloop
        try:
            mleem.em_whileloop = kernels.em_whileloop
            with self.assertRaisesRegex(RuntimeError, "Untested"):
                verify_upstream()
        finally:
            mleem.em_whileloop = original

    def test_skip_threshold_matches_upstream(self):
        genes = {"gene": gene_case(guides=6)}
        Runner()(genes, SimpleNamespace(max_sgrnapergene_permutation=6), nproc=1)
        assert_array_equal(genes["gene"].beta_estimate, np.zeros(8))
        assert_array_equal(genes["gene"].beta_pval, np.ones(2))


class PermutationTests(unittest.TestCase):
    def test_strict_tails_preserve_ties_nan_and_infinity(self):
        null = np.array([[0, 1], [0, 2], [1, 2], [np.nan, np.inf], [-np.inf, 3]])
        observed = np.array([[0, 2], [1, np.inf], [np.nan, 3], [-np.inf, 1]])
        two, upper, lower = kernels.permutation_tails(null, observed)
        for index, beta in enumerate(observed):
            expected_upper = np.sum(null > beta, axis=0) / len(null)
            expected_lower = np.sum(null < beta, axis=0) / len(null)
            assert_array_equal(upper[index], expected_upper)
            assert_array_equal(lower[index], expected_lower)
            assert_array_equal(two[index], 2 * np.minimum(expected_upper, expected_lower))
        with self.assertRaises(ValueError):
            kernels.permutation_tails(np.empty((0, 2)), observed)

    def test_assignment_and_fdr_match(self):
        rng = np.random.default_rng(1)
        null = rng.integers(-3, 4, size=(2000, 2)).astype(float)
        genes = {str(i): gene_case(i) for i in range(60)}
        for gene in genes.values():
            gene.beta_estimate = rng.normal(size=8)
            gene.beta_pval = rng.uniform(size=2)
            gene.beta_pval_pos = gene.beta_pval
            gene.beta_pval_neg = gene.beta_pval
        expected = copy.deepcopy(genes)
        mlemultiprocessing.assign_p_value_from_permuted_beta(null, expected)
        kernels.assign_p_value_from_permuted_beta(null, genes)
        gene_fdr_correction(expected, "fdr")
        gene_fdr_correction(genes, "fdr")
        for key in genes:
            for field in ("beta_permute_pval", "beta_permute_pval_neg", "beta_permute_pval_pos",
                          "beta_permute_pval_fdr", "beta_permute_pval_neg_fdr", "beta_permute_pval_pos_fdr"):
                assert_array_equal(getattr(genes[key], field), getattr(expected[key], field))


class WorkflowTests(unittest.TestCase):
    def test_public_fixture_fit_efficiency_and_mean_variance_match(self):
        fixture = Path(__file__).parent / "data" / "mageck2-counts.tsv"
        content = fixture.read_bytes()
        blob = hashlib.sha1(b"blob " + str(len(content)).encode() + b"\0" + content).hexdigest()
        self.assertEqual(blob, "633f72aaeb99a0b2e7bd40e4742056f0d94340ed")
        with tempfile.TemporaryDirectory() as directory:
            details = []
            for backend in ("reference", "accelerated"):
                prefix = Path(directory) / backend
                completed = subprocess.run([
                    sys.executable, "-m", "crisprworks_fit", "mle", "-k", str(fixture),
                    "-d", str(fixture.parent / "design.tsv"), "-n", str(prefix),
                    "--backend", backend, "--update-efficiency", "--genes-varmodeling", "20",
                    "--write-fit-details", "--seed", "42",
                ], capture_output=True, text=True, timeout=90)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                details.append(json.loads(Path(str(prefix) + ".fit-details.json").read_text())["genes"])
            self.assertEqual(details[0].keys(), details[1].keys())
            for gene_id, expected in details[0].items():
                for field, values in expected.items():
                    assert_allclose(details[1][gene_id][field], values, rtol=1e-7, atol=1e-8,
                                    err_msg=f"{gene_id}: {field}")

    def test_cli_reference_accelerated_and_spawn_workers_match(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            counts, design = synthetic_screen(root / "inputs", genes=16, guides=4)
            summaries = []
            for backend, workers in (("reference", 1), ("accelerated", 1), ("accelerated", 2)):
                prefix = root / f"{backend}-{workers}"
                result = subprocess.run([
                    sys.executable, "-m", "crisprworks_fit", "mle", "-k", str(counts),
                    "-d", str(design), "-n", str(prefix), "--backend", backend,
                    "--threads", str(workers), "--seed", "42", "--permutation-round", "2",
                ], capture_output=True, text=True, timeout=90)
                self.assertEqual(result.returncode, 0, result.stderr)
                summaries.append(Path(str(prefix) + ".gene_summary.txt").read_bytes())
                manifest = json.loads(Path(str(prefix) + ".crisprworks.json").read_text())
                self.assertEqual(manifest["status"], "complete")
                self.assertEqual(manifest["inputs"]["count_table"]["sha256"],
                                 hashlib.sha256(counts.read_bytes()).hexdigest())
            self.assertEqual(summaries[0], summaries[1])
            self.assertEqual(summaries[1], summaries[2])

    def test_invalid_design_and_zero_rounds_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            counts, design = synthetic_screen(root)
            for extra in (("--permutation-round", "0"), ("-d", "1,1;1,0")):
                result = subprocess.run([
                    sys.executable, "-m", "crisprworks_fit", "mle", "-k", str(counts),
                    "-d", str(design), "-n", str(root / "failed"), *extra,
                ], capture_output=True, text=True, timeout=30)
                self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
