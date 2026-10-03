"""Joint-screen inference agrees with an independent pinned upstream oracle."""

import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np
from numpy.testing import assert_allclose

from crisprworks_fit.joint import fit_joint
from crisprworks_fit.joint import run_joint


FIXTURE = Path(__file__).parent / 'fixtures' / 'joint_upstream.json'


def write_table(path, rows):
    with path.open('w', newline='') as stream:
        csv.writer(stream, delimiter='\t').writerows(rows)


class JointTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fixture = json.loads(FIXTURE.read_text())
        self.counts, self.mapping = self.root / 'counts.tsv', self.root / 'map.tsv'
        write_table(self.counts, [self.fixture['header']] + self.fixture['rows'])
        write_table(self.mapping, self.fixture['sample_map'])

    def test_full_pipeline_matches_unmodified_upstream(self):
        result = fit_joint(self.counts, self.mapping, self.fixture['controls'])
        self.assertEqual(result['conditions'], self.fixture['conditions'])
        self.assertEqual(result['count_pseudocount'], 32)
        for gene, expected in self.fixture['expected'].items():
            for field, values in expected.items():
                assert_allclose(result['genes'][gene][field], values, atol=1e-11, rtol=1e-11,
                                err_msg=f'{gene} {field}')
        self.assertEqual(result['effect_units'], 'log2 relative abundance')

    def test_diagnostics_distinguish_stopping_from_iteration_exhaustion(self):
        stopped = fit_joint(self.counts, self.mapping, self.fixture['controls'])
        capped = fit_joint(self.counts, self.mapping, self.fixture['controls'], max_iterations=1,
                          tolerance=1e-12)
        self.assertEqual(stopped['nonconverged_genes'], 0)
        self.assertGreater(capped['nonconverged_genes'], 0)
        for values in stopped['genes'].values():
            diagnostics = values['fit_diagnostics']
            self.assertTrue(diagnostics['converged'])
            self.assertEqual(diagnostics['termination_reason'], 'bound_tolerance')
            self.assertLess(diagnostics['final_bound_change'], .1)
            self.assertGreaterEqual(diagnostics['iterations'], 1)
        for values in capped['genes'].values():
            diagnostics = values['fit_diagnostics']
            self.assertEqual(diagnostics['iterations'], 1)
            self.assertEqual(diagnostics['termination_reason'], 'iteration_limit')

    def test_stopping_settings_are_validated(self):
        for kwargs in ({'max_iterations': 0}, {'max_iterations': True},
                       {'tolerance': 0}, {'tolerance': float('nan')}):
            with self.assertRaises(ValueError):
                fit_joint(self.counts, self.mapping, **kwargs)

    def test_input_changes_during_fitting_preserve_previous_completed_bundle(self):
        args = SimpleNamespace(counts=self.counts, sample_map=self.mapping,
                               control_gene=None, output_prefix=self.root / 'fit')
        summary = run_joint(args)
        manifest = self.root / 'fit.joint-details.json'
        original = summary.read_bytes(), manifest.read_bytes()
        def changing_fit(*arguments, **kwargs):
            result = fit_joint(*arguments, **kwargs)
            self.counts.write_bytes(self.counts.read_bytes() + b'\n')
            return result
        with patch('crisprworks_fit.joint.fit_joint', side_effect=changing_fit):
            with self.assertRaisesRegex(ValueError, 'Input changed'):
                run_joint(args)
        self.assertEqual((summary.read_bytes(), manifest.read_bytes()), original)

    def test_order_invariance_and_all_guide_normalization(self):
        before = fit_joint(self.counts, self.mapping, self.fixture['controls'])
        write_table(self.counts, [self.fixture['header']] + list(reversed(self.fixture['rows'])))
        write_table(self.mapping, [self.fixture['sample_map'][0]] + list(reversed(self.fixture['sample_map'][1:])))
        after = fit_joint(self.counts, self.mapping, self.fixture['controls'])
        self.assertEqual(before, after)
        self.assertEqual(fit_joint(self.counts, self.mapping)['normalization'], 'all-guide median log counts')

    def test_invalid_counts_and_controls_fail_before_fitting(self):
        with self.assertRaisesRegex(ValueError, 'absent'):
            fit_joint(self.counts, self.mapping, ['missing'])
        with self.assertRaisesRegex(ValueError, 'empty'):
            fit_joint(self.counts, self.mapping, [])
        for value in ('nan', 'inf', '-1', '1+2'):
            rows = json.loads(json.dumps(self.fixture['rows']))
            rows[0][2] = value
            write_table(self.counts, [self.fixture['header']] + rows)
            with self.assertRaises(ValueError):
                fit_joint(self.counts, self.mapping)
        write_table(self.counts, [self.fixture['header']] + self.fixture['rows'][:60])
        with self.assertRaisesRegex(ValueError, '64 guides'):
            fit_joint(self.counts, self.mapping)

    def test_invalid_sample_mapping_is_rejected(self):
        rows = json.loads(json.dumps(self.fixture['sample_map']))
        rows[1][2] = 'A'
        write_table(self.mapping, rows)
        with self.assertRaisesRegex(ValueError, 'self-mapped'):
            fit_joint(self.counts, self.mapping)
        rows = json.loads(json.dumps(self.fixture['sample_map']))
        rows[5][1] = 'SINGLETON'
        write_table(self.mapping, rows)
        with self.assertRaisesRegex(ValueError, 'two replicates'):
            fit_joint(self.counts, self.mapping)
        write_table(self.mapping, self.fixture['sample_map'] + [self.fixture['sample_map'][1]])
        with self.assertRaisesRegex(ValueError, 'exactly once'):
            fit_joint(self.counts, self.mapping)

    def test_cli_joint_and_gene_calibration_preserve_scale_and_family(self):
        controls = self.root / 'controls.txt'
        controls.write_text('\n'.join(self.fixture['controls']) + '\n')
        command = [sys.executable, '-m', 'crisprworks_fit', 'joint', '-k', str(self.counts),
                   '--sample-map', str(self.mapping), '--control-gene', str(controls), '-n', str(self.root / 'fit')]
        completed = subprocess.run(command, text=True, capture_output=True, timeout=30)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        details = self.root / 'fit.joint-details.json'
        self.assertTrue(details.exists())
        calibration = subprocess.run([sys.executable, '-m', 'crisprworks_fit', 'calibrate', '--fit-details', str(details),
                                      '--control-gene', str(controls), '--min-controls', '2', '-n', str(self.root / 'cal')],
                                     text=True, capture_output=True, timeout=30)
        self.assertEqual(calibration.returncode, 0, calibration.stderr)
        manifest = json.loads((self.root / 'cal.calibration.json').read_text())
        self.assertEqual(manifest['family_size'], (24 - 6) * 2)
        self.assertEqual(manifest['condition_labels'], ['A', 'B'])
        self.assertEqual(manifest['effect_units'], 'log2 relative abundance')
        original = (self.root / 'fit.joint.tsv').read_bytes()
        (self.root / 'fit.joint.tsv').write_text('tampered summary\n')
        before_calibration = (self.root / 'cal.calibration.tsv').read_bytes()
        rejected = subprocess.run([sys.executable, '-m', 'crisprworks_fit', 'calibrate',
            '--fit-details', str(details), '--control-gene', str(controls), '--min-controls', '2',
            '-n', str(self.root / 'cal')], text=True, capture_output=True, timeout=30)
        self.assertNotEqual(rejected.returncode, 0)
        self.assertIn('hash mismatch', rejected.stderr)
        self.assertEqual((self.root / 'cal.calibration.tsv').read_bytes(), before_calibration)
        (self.root / 'fit.joint.tsv').write_bytes(original)
        controls.write_text('missing\n')
        bad = subprocess.run(command, text=True, capture_output=True, timeout=30)
        self.assertNotEqual(bad.returncode, 0)
        self.assertEqual((self.root / 'fit.joint.tsv').read_bytes(), original)


if __name__ == '__main__':
    unittest.main()
