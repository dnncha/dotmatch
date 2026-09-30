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
            self.assertIn(str(root / "counts with spaces.tsv"), tasks[0][1])
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
