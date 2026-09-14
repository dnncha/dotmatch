"""Fail-closed archive and trace regression tests; no network or raw inputs."""
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import unittest
import summarize


class SummaryTests(unittest.TestCase):
    def make_artifact(self, root):
        data = root / 'evidence.txt'
        data.write_text('evidence\n')
        manifest = {'status': 'complete', 'files': {'evidence.txt': {
            'bytes': data.stat().st_size, 'sha256': summarize.sha(data)}}}
        (root / 'MANIFEST.json').write_text(json.dumps(manifest))
        return manifest

    def test_valid_and_corrupted_artifact(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.make_artifact(root)
            self.assertEqual(summarize.verify(root)['files_verified'], 1)
            (root / 'evidence.txt').write_text('Evidence\n')
            with self.assertRaisesRegex(ValueError, 'corrupt'):
                summarize.verify(root)

    def test_incomplete_and_missing_artifact(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = self.make_artifact(root)
            manifest['status'] = 'failed'
            (root / 'MANIFEST.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'incomplete'):
                summarize.verify(root)
            self.make_artifact(root)
            (root / 'evidence.txt').unlink()
            with self.assertRaisesRegex(ValueError, 'missing'):
                summarize.verify(root)

    def test_path_escape_is_rejected(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = self.make_artifact(root)
            manifest['files']['../escaped.txt'] = manifest['files']['evidence.txt']
            (root / 'MANIFEST.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'unsafe'):
                summarize.verify(root)

    def test_trace_gene_and_overlap_semantics(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / 'traces.tsv'
            path.write_text('record_ordinal\toffset_and_guide\tdistinct_guides\n'
                            '1\t21:a;22:b\t2\n2\t21:a;23:c\t2\n')
            result = summarize.annotate_trace(path, {'a':'G','b':'G','c':'H'}, {(1,'a','b')})
            self.assertEqual(result['multi_offset_records'], 2)
            self.assertEqual(result['same_gene'], 1)
            self.assertEqual(result['cross_gene'], 1)
            self.assertEqual(result['records_with_exact_overlap_pair'], 1)
            self.assertEqual(result['observed_pairs'], 2)

    def test_multiple_calls_at_one_offset_rejected(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / 'traces.tsv'
            path.write_text('record_ordinal\toffset_and_guide\tdistinct_guides\n1\t21:a;21:b\t2\n')
            with self.assertRaisesRegex(ValueError, 'invalid'):
                summarize.annotate_trace(path, {'a':'G','b':'H'}, set())


if __name__ == '__main__':
    unittest.main(verbosity=2)
