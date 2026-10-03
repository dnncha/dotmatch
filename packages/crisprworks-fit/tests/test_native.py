"""Compiled engine contracts and direct comparison with the NumPy accelerator."""
import copy
import unittest
import warnings
from unittest.mock import patch
import numpy as np
from numpy.testing import assert_allclose
from mageck2 import mleem
from crisprworks_fit import kernels
from crisprworks_fit.backend import accelerated
from test_fit import gene_case, VarianceModel


class NativeTests(unittest.TestCase):
    def test_selection_and_state_restoration(self):
        previous = kernels._engine
        with accelerated(kernel="numpy"):
            self.assertEqual(kernels._engine, "numpy")
        self.assertEqual(kernels._engine, previous)
        with patch.object(kernels, "_native", None):
            self.assertEqual(kernels.resolved_engine("auto"), "numpy")
            with self.assertRaisesRegex(RuntimeError, "unavailable"):
                kernels.resolved_engine("native")
        with self.assertRaises(ValueError):
            kernels.resolved_engine("other")

    @unittest.skipIf(kernels._native is None, "Native extension not installed")
    def test_native_matches_numpy_for_matrix_designs_and_low_counts(self):
        for guides in (1, 4, 12, 24):
            for seed in (2, 17):
                original = gene_case(seed, guides)
                original.design_mat = np.matrix(original.design_mat)
                original.nb_count += 1
                original.nb_count[-1, 0] = 1
                original.nb_count[-2, 0] = 2
                settings = dict(debug=False, logem=False, estimateeff=True,
                                updateeff=True, meanvarmodel=VarianceModel())
                expected, actual = copy.deepcopy(original), copy.deepcopy(original)
                with accelerated(kernel="numpy"):
                    mleem.iteratenbem(expected, **settings)
                with accelerated(kernel="native"):
                    mleem.iteratenbem(actual, **settings)
                for field in ("beta_estimate", "beta_zscore", "beta_pval", "w_estimate",
                              "mu_estimate", "sgrna_residule"):
                    assert_allclose(getattr(actual, field), getattr(expected, field),
                                    rtol=1e-7, atol=1e-8, err_msg=field)

    @unittest.skipIf(kernels._native is None, "Native extension not installed")
    def test_invalid_buffers_are_rejected_before_computation(self):
        arrays = [np.ones((3, 2)), np.ones(3), np.ones(3), np.ones(2), np.ones(1), np.ones(3)]
        for index in range(6):
            bad = list(arrays)
            bad[index] = np.ones(1, dtype=np.float32)
            with self.assertRaises(ValueError):
                kernels._native.loop(*bad, 1, 1, 1, .01, True, True)
        unaligned = np.ndarray((3,), dtype=np.float64, buffer=bytearray(25), offset=1)
        arrays[1] = unaligned
        with self.assertRaises(ValueError):
            kernels._native.loop(*arrays, 1, 1, 1, .01, True, True)
        with self.assertRaises(ValueError):
            kernels._native.loop(*arrays, -1, 1, 1, .01, True, True)

    @unittest.skipIf(kernels._native is None, "Native extension not installed")
    def test_nonfinite_initialization_matches_numpy_failure_conventions(self):
        original = gene_case(guides=4)
        original.nb_count[0, 0] = 0
        expected, actual = copy.deepcopy(original), copy.deepcopy(original)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            for kernel, gene in (("numpy", expected), ("native", actual)):
                with accelerated(kernel=kernel):
                    mleem.iteratenbem(gene, debug=False, logem=False, estimateeff=True, updateeff=True)
        for field in ("beta_estimate", "beta_zscore", "beta_pval", "w_estimate",
                      "mu_estimate", "sgrna_residule"):
            assert_allclose(getattr(actual, field), getattr(expected, field),
                            rtol=1e-7, atol=1e-8, equal_nan=True)

    @unittest.skipIf(kernels._native is None, "Native extension not installed")
    def test_native_singular_solve_rejects_without_extra_regularization(self):
        arrays = [np.ones((3, 2)), np.full(3, 10.), np.ones(3),
                  np.ones(2), np.ones(1), np.full(3, .1)]
        with self.assertRaisesRegex(ArithmeticError, "Singular"):
            kernels._native.loop(*arrays, 1, 1, 1, 0., False, False)
