"""Gene-control calibration and declared multiple-testing families."""

import copy
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from numpy.testing import assert_allclose
from scipy.stats import false_discovery_control

from crisprworks_fit.calibration import adjust_fdr, calibrate_genes


def gene(effects, guides=4):
    return {'guides': guides, 'beta_estimate': [0.] * guides + list(effects), 'beta_zscore': list(effects)}


class CalibrationTests(unittest.TestCase):
    def test_resolution_reports_impossible_discovery_despite_extreme_effects(self):
        genes = {f'control{i}': gene([float(i)]) for i in range(20)}
        controls = list(genes)
        genes.update({f'test{i}': gene([-1000.]) for i in range(500)})
        result = calibrate_genes(genes, controls)
        resolution = result['resolution']
        expected = false_discovery_control(np.full(500, 1 / 21), method='by')[0]
        self.assertAlmostEqual(resolution['minimum_attainable_directional_q'], expected)
        self.assertEqual(resolution['directional_hypotheses_with_attainable_q'], 0)
        self.assertEqual(resolution['two_sided_hypotheses_with_attainable_q'], 0)
        self.assertGreater(expected, .05)

    def test_resolution_uses_the_declared_family_without_changing_qvalues(self):
        genes = {'a': gene([0., 0.]), 'b': gene([1., 1.]), 'test': gene([-100., -100.])}
        for method in ('bh', 'by'):
            for family in ('global', 'condition'):
                before = calibrate_genes(genes, ['a', 'b'], min_controls=2, method=method, family=family)
                after = calibrate_genes(genes, ['a', 'b'], min_controls=2, method=method,
                                        family=family, alpha=.9)
                self.assertEqual(before['rows'], after['rows'])
                row = next(row for row in before['rows'] if row['gene'] == 'test')
                self.assertAlmostEqual(before['resolution']['minimum_attainable_directional_q'], row['q_negative'])
                self.assertAlmostEqual(before['resolution']['minimum_attainable_two_sided_q'], row['q_two'])
                self.assertEqual(after['resolution']['directional_hypotheses_with_attainable_q'], 2)

    def test_resolution_keeps_unsupported_hypotheses_and_handles_empty_families(self):
        genes = {'a': gene([0.]), 'b': gene([1.]), 'test': gene([-100.], guides=5)}
        result = calibrate_genes(genes, ['a', 'b'], min_controls=2)
        self.assertEqual(result['resolution']['minimum_attainable_directional_q'], 1.)
        empty = calibrate_genes({'a': gene([0.]), 'b': gene([1.])}, ['a', 'b'], min_controls=2)
        self.assertIsNone(empty['resolution']['minimum_attainable_directional_q'])
        for alpha in (0, -1, 1.1, float('nan'), True):
            with self.assertRaises(ValueError):
                calibrate_genes(genes, ['a', 'b'], min_controls=2, alpha=alpha)

    def test_explicit_nonconvergence_and_malformed_diagnostics_are_rejected(self):
        for diagnostics in ({'converged': False}, {'converged': 'false'}, {}, None,
                            {'converged': True},
                            {'converged': True, 'iterations': 0, 'termination_reason': 'bound_tolerance', 'final_bound_change': .01},
                            {'converged': True, 'iterations': 1, 'termination_reason': 'iteration_limit', 'final_bound_change': .01},
                            {'converged': True, 'iterations': 1, 'termination_reason': 'bound_tolerance', 'final_bound_change': float('nan')}):
            genes = {'a': gene([0.]), 'b': gene([1.]), 'test': gene([-10.])}
            genes['test']['fit_diagnostics'] = diagnostics
            with self.assertRaisesRegex(ValueError, 'fit diagnostics|Nonconverged fit'):
                calibrate_genes(genes, ['a', 'b'], min_controls=2)

    def test_successful_joint_fits_are_not_excluded_by_mle_permutation_limit(self):
        for size in (40, 80):
            genes = {name: {'guides': size, 'effect_estimate': [effect],
                            'effect_zscore': [effect]}
                     for name, effect in (('a', 0.), ('b', 1.), ('test', -10.))}
            result = calibrate_genes(genes, ['a', 'b'], min_controls=2)
            row = next(row for row in result['rows'] if row['gene'] == 'test')
            self.assertEqual(row['status'], 'tested')
            self.assertEqual(row['p_negative'], 1 / 3)
            self.assertIsNone(result['maximum_guides_exclusive'])

    def test_bh_and_by_match_independent_scipy_oracle_with_ties(self):
        rng = np.random.default_rng(1729)
        for size in (1, 5, 40, 1000):
            values = rng.integers(0, 100, size=size) / 100
            for method in ('bh', 'by'):
                assert_allclose(adjust_fdr(values, method), false_discovery_control(values, method=method),
                                rtol=2e-15, atol=2e-15)
        for bad in ([np.nan], [-.1], [1.1], [[.1]]):
            with self.assertRaises(ValueError):
                adjust_fdr(bad)

    def test_strata_and_training_controls_cannot_be_called(self):
        genes = {'a': gene([0.]), 'b': gene([1.]), 'test': gene([-10.]),
                 'u': gene([100.], guides=5), 'v': gene([101.], guides=5),
                 'other': gene([102.], guides=5)}
        untouched = copy.deepcopy(genes)
        result = calibrate_genes(genes, ['a', 'b', 'u', 'v'], min_controls=2, method='bh')
        rows = {row['gene']: row for row in result['rows']}
        self.assertEqual(rows['test']['p_negative'], 1 / 3)
        self.assertEqual(rows['other']['p_negative'], 1.)
        self.assertEqual(rows['other']['p_positive'], 1 / 3)
        for control in ('a', 'b', 'u', 'v'):
            self.assertEqual(rows[control]['p_negative'], 1.)
            self.assertEqual(rows[control]['q_negative'], 1.)
            self.assertEqual(rows[control]['status'], 'training_control')
        self.assertEqual(genes, untouched)

    def test_missing_support_failed_scores_and_skipped_fits_are_explicit(self):
        genes = {'a': gene([0., 0.]), 'b': gene([1., np.nan]), 'test': gene([-10., -10.]),
                 'invalid': gene([np.nan, -1.]), 'unsupported': gene([-100., -100.], guides=5),
                 'skipped': gene([-100., -100.], guides=40)}
        result = calibrate_genes(genes, ['a', 'b'], min_controls=2)
        rows = {(row['gene'], row['condition']): row for row in result['rows']}
        self.assertEqual(rows['test', 1]['p_negative'], 1 / 3)
        self.assertEqual(rows['test', 2]['status'], 'invalid_control_fit')
        self.assertEqual(rows['test', 2]['p_negative'], 1.)
        self.assertEqual(rows['invalid', 1]['status'], 'invalid_fit')
        self.assertIsNone(rows['invalid', 1]['score'])
        self.assertEqual(rows['unsupported', 1]['status'], 'insufficient_controls')
        self.assertEqual(rows['skipped', 1]['status'], 'skipped_fit')
        self.assertEqual(result['family_size'], 8)
        with self.assertRaisesRegex(ValueError, 'absent'):
            calibrate_genes(genes, ['missing'])
        mixed = copy.deepcopy(genes)
        mixed['a'] = {'guides': 4, 'effect_estimate': [0., 0.], 'effect_zscore': [0., 0.]}
        with self.assertRaisesRegex(ValueError, 'incompatible'):
            calibrate_genes(mixed, ['a', 'b'])
        genes['test']['beta_estimate'] = [-1.]
        with self.assertRaisesRegex(ValueError, 'dimensions'):
            calibrate_genes(genes, ['a', 'b'])

    def test_global_family_adjusts_all_conditions_and_ties_are_conservative(self):
        genes = {'a': gene([0., 0.]), 'b': gene([0., 0.]), 'test': gene([0., -10.])}
        result = calibrate_genes(genes, ['a', 'b'], min_controls=2, method='by')
        row = next(row for row in result['rows'] if row['gene'] == 'test' and row['condition'] == 1)
        self.assertEqual(row['p_negative'], 1.)
        self.assertEqual(row['p_two'], 1.)
        self.assertEqual(row['q_negative'], 1.)
        condition = calibrate_genes(genes, ['a', 'b'], min_controls=2, family='condition')
        self.assertEqual(condition['family_size'], 1)

    def test_cli_roundtrip_hashes_and_bad_input_does_not_replace_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            details = root / 'screen.fit-details.json'
            details.write_text(json.dumps({'schema_version': 1,
                                          'genes': {'a': gene([0.]), 'b': gene([1.]), 'test': gene([-10.])}}))
            controls = root / 'controls.txt'
            controls.write_text('a\nb\n')
            prefix = root / 'calibrated'
            command = [sys.executable, '-m', 'crisprworks_fit', 'calibrate', '--fit-details', str(details),
                       '--control-gene', str(controls), '-n', str(prefix), '--min-controls', '2']
            completed = subprocess.run(command, text=True, capture_output=True, timeout=30)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = json.loads((root / 'calibrated.calibration.json').read_text())
            self.assertEqual(result['fdr_method'], 'by')
            self.assertEqual(result['testing_family'], 'global')
            output = root / 'calibrated.calibration.tsv'
            original = output.read_bytes()
            self.assertEqual(result['output']['sha256'], hashlib.sha256(original).hexdigest())
            rows = list(csv.DictReader(output.open(), delimiter='\t'))
            self.assertEqual(len(rows), 3)
            controls.write_text('missing\n')
            failed = subprocess.run(command, text=True, capture_output=True, timeout=30)
            self.assertNotEqual(failed.returncode, 0)
            self.assertEqual(output.read_bytes(), original)

    def test_condition_names_survive_calibration_and_invalid_labels_preserve_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            details, controls = root / 'details.json', root / 'controls.txt'
            data = {'schema_version': 1, 'model': 'MAGeCK2 MLE', 'effect_units': 'natural-log beta coefficient',
                    'conditions': ['drug', 'recovery'],
                    'genes': {'a': gene([0., 0.]), 'b': gene([1., 1.]), 'test': gene([-2., 2.])}}
            details.write_text(json.dumps(data))
            controls.write_text('a\nb\n')
            command = [sys.executable, '-m', 'crisprworks_fit', 'calibrate', '--fit-details', str(details),
                       '--control-gene', str(controls), '--min-controls', '2', '-n', str(root / 'out')]
            result = subprocess.run(command, text=True, capture_output=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            output = root / 'out.calibration.tsv'
            original = output.read_bytes()
            with output.open() as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            self.assertEqual({r['condition_label'] for r in rows}, {'drug', 'recovery'})
            for labels in ('ab', ['drug'], ['drug', 'drug'], ['', 'recovery'], [None, 'recovery']):
                data['conditions'] = labels
                details.write_text(json.dumps(data))
                bad = subprocess.run(command, text=True, capture_output=True, timeout=30)
                self.assertNotEqual(bad.returncode, 0)
                self.assertEqual(output.read_bytes(), original)


if __name__ == '__main__':
    unittest.main()
