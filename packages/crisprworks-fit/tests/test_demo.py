"""The documented joint demo exercises the installed workflow end to end."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class JointDemoTests(unittest.TestCase):
    def test_demo_calibrates_declared_controls_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "joint-demo"
            command = [sys.executable, "-m", "crisprworks_fit", "demo", "--model", "joint", "--out-dir", str(output)]
            result = subprocess.run(command, text=True, capture_output=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = json.loads((output / "demo.json").read_text())
            self.assertEqual(manifest["status"], "complete")
            for record in manifest["inputs"].values():
                self.assertEqual(record["sha256"], hashlib.sha256(Path(record["path"]).read_bytes()).hexdigest())
            for record in manifest["outputs"].values():
                self.assertEqual(record["sha256"], hashlib.sha256(Path(record["path"]).read_bytes()).hexdigest())
            calibration = json.loads((output / "screen.calibration.json").read_text())
            self.assertEqual(calibration["testing_family"], "global")
            self.assertEqual(calibration["fdr_method"], "by")
            self.assertEqual(calibration["family_size"], 32)
            self.assertEqual(calibration["condition_labels"], ["CONTROL", "TREATED"])
            with (output / "screen.calibration.tsv").open() as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(len(rows), 96)
            training = [r for r in rows if r["status"] == "training_control"]
            self.assertEqual(len(training), 64)
            self.assertTrue(all(float(r["q_two"]) == 1 for r in training))
            self.assertTrue(all(float(r["p_two"]) > 0 for r in rows))
            original = (output / "demo.json").read_bytes()
            failed = subprocess.run(command, text=True, capture_output=True, timeout=30)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("already exists", failed.stderr)
            self.assertEqual((output / "demo.json").read_bytes(), original)

    def test_joint_demo_rejects_mle_execution_options(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "joint-demo"
            result = subprocess.run([sys.executable, "-m", "crisprworks_fit", "demo", "--model", "joint",
                                     "--threads", "2", "--out-dir", str(output)],
                                    text=True, capture_output=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("MLE demo behavior", result.stderr)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
