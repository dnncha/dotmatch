"""Statistical safeguards, separate from legacy MAGeCK numerical parity."""

import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
from numpy.testing import assert_array_equal, assert_allclose
from mageck2 import mlemultiprocessing

from crisprworks_fit.backend import accelerated, permutation_calibration
from crisprworks_fit.cli import synthetic_screen
from crisprworks_fit.inference import finite_permutation_tails, assign_finite_permutation_pvalues


class FiniteTailTests(unittest.TestCase):
    def test_ties_extremes_and_failed_fits_cannot_be_zero_p(self):
        null = np.array([[0., -1., 0.], [0., 0., np.nan], [0., 1., 2.]])
        observed = np.array([[0., -2., -100.], [0., 2., 100.], [np.nan, np.inf, 0.]])
        two, upper, lower = finite_permutation_tails(null, observed)
        assert_array_equal(two[:, 0], [1., 1., 1.])
        assert_array_equal(two[:, 1], [.5, .5, 1.])
        assert_array_equal(two[:, 2], [1., 1., 1.])
        assert_array_equal(lower[:, 1], [.25, 1., 1.])
        assert_array_equal(upper[:, 1], [1., .25, 1.])
        for bad in (np.nan, np.inf, -np.inf):
            tails = finite_permutation_tails([[bad]], [[0.]])
            for tail in tails:
                assert_array_equal(tail, [[1.]])

    def test_exact_enumeration_is_superuniform_for_exchangeable_draws_with_ties(self):
        # Enumerate every equally likely observation + three null draws from a
        # discrete null. This checks Type I error, not just our formula.
        pvalues = [[], [], []]
        for sample in itertools.product((-1., 0., 1.), repeat=4):
            tails = finite_permutation_tails(np.array(sample[1:])[:, None], [[sample[0]]])
            for values, tail in zip(pvalues, tails):
                values.append(tail.item())
        for values in pvalues:
            for alpha in (0., .05, .25, .5, .75, 1.):
                self.assertLessEqual(np.mean(np.array(values) <= alpha), alpha + 1e-15)

    def test_vectorized_tails_match_independent_counting_oracle(self):
        rng = np.random.default_rng(88)
        null = rng.integers(-4, 5, size=(73, 3))
        observed = rng.integers(-6, 7, size=(41, 3))
        two, upper, lower = finite_permutation_tails(null, observed)
        for i, row in enumerate(observed):
            for j, value in enumerate(row):
                lo = (1 + sum(draw <= value for draw in null[:, j])) / 74
                hi = (1 + sum(draw >= value for draw in null[:, j])) / 74
                self.assertEqual(lower[i, j], lo)
                self.assertEqual(upper[i, j], hi)
                self.assertEqual(two[i, j], min(1., 2 * min(lo, hi)))
        with self.assertRaises(ValueError):
            finite_permutation_tails(np.empty((0, 3)), observed)
        with self.assertRaises(ValueError):
            finite_permutation_tails(null, np.ones((2, 2)))

    def test_skipped_genes_and_resolution_are_recorded(self):
        genes = {str(n): SimpleNamespace(nb_count=np.ones((4, n)),
                                       beta_estimate=np.r_[np.zeros(n), -100.]) for n in (4, 40)}
        diagnostics = []
        assign_finite_permutation_pvalues(np.arange(9.)[:, None], genes,
                                         max_guides=40, diagnostics=diagnostics)
        assert_array_equal(genes['4'].beta_permute_pval_neg, [.1])
        assert_array_equal(genes['40'].beta_permute_pval_neg, [1.])
        self.assertEqual(diagnostics[0]['skipped_genes'], 1)
        self.assertEqual(diagnostics[0]['minimum_directional_p'], .1)

    def test_context_restores_both_reference_and_accelerated_after_failure(self):
        original = mlemultiprocessing.assign_p_value_from_permuted_beta
        with self.assertRaisesRegex(RuntimeError, 'injected'):
            with permutation_calibration('finite'):
                self.assertIsNot(mlemultiprocessing.assign_p_value_from_permuted_beta, original)
                with self.assertRaisesRegex(RuntimeError, 'Nested'):
                    with permutation_calibration('finite'):
                        pass
                raise RuntimeError('injected')
        self.assertIs(mlemultiprocessing.assign_p_value_from_permuted_beta, original)
        with accelerated():
            accelerated_function = mlemultiprocessing.assign_p_value_from_permuted_beta
            with permutation_calibration('finite'):
                self.assertIsNot(mlemultiprocessing.assign_p_value_from_permuted_beta, accelerated_function)
            self.assertIs(mlemultiprocessing.assign_p_value_from_permuted_beta, accelerated_function)
        self.assertIs(mlemultiprocessing.assign_p_value_from_permuted_beta, original)

    def test_cli_records_inference_and_preserves_estimates(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            counts, design = synthetic_screen(root, genes=12, guides=4)
            details = {}
            for backend, mode in (('accelerated', 'legacy'), ('accelerated', 'finite'), ('reference', 'finite')):
                prefix = root / f'{backend}-{mode}'
                command = [sys.executable, '-m', 'crisprworks_fit', 'mle', '-k', str(counts),
                           '-d', str(design), '-n', str(prefix), '--backend', backend,
                           '--permutation-pvalues', mode, '--permutation-round', '2',
                           '--genes-varmodeling', '0', '--seed', '42', '--write-fit-details']
                completed = subprocess.run(command, capture_output=True, text=True, timeout=60)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                details[backend, mode] = json.loads(Path(str(prefix) + '.fit-details.json').read_text())['genes']
                manifest = json.loads(Path(str(prefix) + '.crisprworks.json').read_text())
                self.assertEqual(manifest['options']['permutation_pvalues'], mode)
                if mode == 'finite':
                    self.assertEqual(manifest['permutation_calibration']['groups'][0]['null_draws'], 24)
            for gene, legacy in details['accelerated', 'legacy'].items():
                for key in ('beta_estimate', 'beta_zscore', 'beta_pval', 'w_estimate'):
                    assert_array_equal(details['accelerated', 'finite'][gene][key], legacy[key])
                for key, values in details['accelerated', 'finite'][gene].items():
                    if key != 'guides':
                        assert_allclose(details['reference', 'finite'][gene][key], values, rtol=1e-7, atol=1e-8)
                for key in ('beta_permute_pval', 'beta_permute_pval_pos', 'beta_permute_pval_neg'):
                    self.assertTrue(np.all(np.array(details['accelerated', 'finite'][gene][key]) >= 1 / 25))
            spawned_prefix = root / 'finite-spawned'
            spawned = subprocess.run([
                sys.executable, '-m', 'crisprworks_fit', 'mle', '-k', str(counts), '-d', str(design),
                '-n', str(spawned_prefix), '--threads', '2', '--permutation-pvalues', 'finite',
                '--permutation-round', '2', '--genes-varmodeling', '0', '--seed', '42', '--write-fit-details',
            ], capture_output=True, text=True, timeout=60)
            self.assertEqual(spawned.returncode, 0, spawned.stderr)
            spawned_details = json.loads(Path(str(spawned_prefix) + '.fit-details.json').read_text())['genes']
            self.assertEqual(spawned_details, details['accelerated', 'finite'])

    def test_control_null_is_recorded_and_ungrouped_controls_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            counts, design = synthetic_screen(root, genes=12, guides=4)
            controls = root / 'controls.txt'
            controls.write_text('GENE2\nGENE3\nGENE4\nGENE5\n')
            prefix = root / 'control'
            command = [sys.executable, '-m', 'crisprworks_fit', 'mle', '-k', str(counts),
                       '-d', str(design), '-n', str(prefix), '--permutation-pvalues', 'finite',
                       '--permutation-round', '2', '--genes-varmodeling', '0', '--norm-method', 'control',
                       '--control-gene', str(controls)]
            completed = subprocess.run(command, capture_output=True, text=True, timeout=60)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            manifest = json.loads(Path(str(prefix) + '.crisprworks.json').read_text())
            self.assertEqual(manifest['permutation_calibration']['null_source'], 'control guides')
            self.assertIsNotNone(manifest['inputs']['control_gene']['sha256'])
            rejected = subprocess.run(command + ['--no-permutation-by-group'], capture_output=True,
                                      text=True, timeout=60)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn('requires grouped permutations', rejected.stderr)


if __name__ == '__main__':
    unittest.main()
