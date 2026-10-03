"""Batch preflight must reject collisions before launching any screen."""

import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("fit_batch", Path(__file__).parents[1] / "examples" / "batch.py")
batch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(batch)


class BatchTests(unittest.TestCase):
    def test_paths_commands_and_preflight_collisions(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "counts with spaces.tsv").touch()
            (root / "design.tsv").touch()
            manifest = root / "screens.json"
            screen = {"name": "screen-a", "count_table": "counts with spaces.tsv", "design_matrix": "design.tsv"}
            manifest.write_text(json.dumps([screen]))
            tasks = batch.commands(manifest, root / "results", 2, 42, 10)
            self.assertIn(str((root / "counts with spaces.tsv").resolve()), tasks[0][1])
            manifest.write_text(json.dumps([screen, screen]))
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                batch.commands(manifest, root / "results", 2, 42, 10)
            manifest.write_text(json.dumps([{**screen, "name": "../escape"}]))
            with self.assertRaisesRegex(ValueError, "Screen names"):
                batch.commands(manifest, root / "results", 2, 42, 10)
            manifest.write_text(json.dumps([screen]))
            (root / "results" / "screen-a").mkdir(parents=True)
            with self.assertRaisesRegex(ValueError, "exists"):
                batch.commands(manifest, root / "results", 2, 42, 10)
            self.assertEqual(len(batch.commands(manifest, root / "results", 2, 42, 10, True)), 1)

    def test_finite_control_inference_and_missing_controls_are_preflighted(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for name in ('counts.tsv', 'design.tsv', 'controls.txt'):
                (root / name).touch()
            manifest = root / 'screens.json'
            manifest.write_text(json.dumps([{'name': 'a', 'count_table': 'counts.tsv',
                                            'design_matrix': 'design.tsv', 'control_gene': 'controls.txt'}]))
            command = batch.commands(manifest, root / 'results', 1, 42, 2,
                                     permutation_pvalues='finite')[0][1]
            self.assertIn('--control-gene', command)
            self.assertEqual(command[command.index('--norm-method') + 1], 'control')
            self.assertEqual(command[command.index('--permutation-pvalues') + 1], 'finite')
            self.assertNotIn('--write-fit-details', command)
            full = batch.commands(manifest, root / 'results', 1, 42, 2,
                                  permutation_pvalues='finite', write_fit_details=True, update_efficiency=True)[0][1]
            self.assertIn('--write-fit-details', full)
            self.assertIn('--update-efficiency', full)
            (root / 'controls.txt').unlink()
            with self.assertRaisesRegex(ValueError, 'Missing control-gene'):
                batch.commands(manifest, root / 'results', 1, 42, 2, permutation_pvalues='finite')
