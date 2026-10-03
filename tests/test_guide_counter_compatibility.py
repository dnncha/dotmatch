#!/usr/bin/env python3
"""Independent byte-oracle regressions for guide-counter compatibility."""

import csv
import gzip
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get('DOTMATCH_BIN', ROOT / 'dotmatch')).resolve()


def oracle(targets, reads, k, sample_size, fraction):
    targets = [seq.upper() for seq in targets]
    def match(window):
        if set(window) - set('ACGT'):
            return None
        distances = [sum(a != b for a, b in zip(window, target)) for target in targets]
        best = min(distances)
        hits = [i for i, distance in enumerate(distances) if distance == best]
        return hits[0] if best <= k and len(hits) == 1 else None
    length = len(targets[0])
    scores = {}
    for read in reads[:sample_size]:
        for offset in range(min(500, len(read) - length + 1)):
            if match(read[offset:offset + length]) is not None:
                scores[offset] = scores.get(offset, 0) + 1
    matched = sum(scores.values())
    offsets = [offset for offset, n in scores.items() if n / matched >= fraction]
    counts = [0] * len(targets)
    for read in reads:
        for offset in offsets:
            if offset + length <= len(read):
                target = match(read[offset:offset + length])
                if target is not None:
                    counts[target] += 1
    return counts


class Compatibility(unittest.TestCase):
    def check_counts(self, targets, reads, *, k=1, sample_size=100000, fraction=.0025):
        expected = oracle(targets, reads, k, sample_size, fraction)
        for compressed in (False, True):
            with self.subTest(compressed=compressed), tempfile.TemporaryDirectory() as temporary:
                folder = Path(temporary)
                library = folder / 'library.tsv'
                library.write_text('guide\tsequence\tgene\n' + ''.join(
                    f'g{i}\t{target}\tGENE{i}\n' for i, target in enumerate(targets)))
                fastq = folder / ('reads.fastq.gz' if compressed else 'reads.fastq')
                payload = ''.join(f'@r{i}\n{seq}\n+\n{"I" * len(seq)}\n' for i, seq in enumerate(reads))
                if compressed:
                    with gzip.open(fastq, 'wt') as stream:
                        stream.write(payload)
                else:
                    fastq.write_text(payload)
                command = [str(BIN), 'guide-counter', 'count', '--input', str(fastq),
                           '--samples', 'sample', '--library', str(library), '--output', str(folder / 'out'),
                           '--offset-sample-size', str(sample_size), '--offset-min-fraction', str(fraction)]
                if k == 0:
                    command.append('--exact-match')
                subprocess.run(command, capture_output=True, text=True, check=True)
                with (folder / 'out.counts.txt').open() as stream:
                    actual = [int(row['sample']) for row in csv.DictReader(stream, delimiter='\t')]
                self.assertEqual(actual, expected)
                with (folder / 'out.stats.txt').open() as stream:
                    stats = next(csv.DictReader(stream, delimiter='\t'))
                self.assertEqual(int(stats['total_reads']), len(reads))
                self.assertEqual(int(stats['mapped_reads']), sum(expected))
                self.assertEqual(int(stats['zero_read_guides']), expected.count(0))

    def test_one_mismatch_can_establish_an_offset(self):
        self.check_counts(['ACGTCAGT'], ['NNACGTCAGA', 'NNACGTCAGT'], sample_size=1)

    def test_n_and_lowercase_windows_are_rejected(self):
        self.check_counts(['acgtcagt'], ['ACGTCAGT', 'ACGTCAGN', 'acgtcagt'])

    def test_each_selected_window_is_counted_even_for_the_same_guide(self):
        self.check_counts(['ACGTCAGT', 'TTAACCGG'],
                          ['ACGTCAGTNNACGTCAGT', 'ACGTCAGTNN TTAACCGG'.replace(' ', '')])

    def test_offsets_use_the_matched_window_denominator(self):
        self.check_counts(['ACGTCAGT', 'TTAACCGG'],
                          ['ACGTCAGTNNNNNNNN', 'NNNNNNNNTTAACCGG'] + ['N' * 16] * 8,
                          fraction=.4)

    def test_no_fallback_when_no_offset_passes_the_threshold(self):
        self.check_counts(['ACGTCAGT', 'TTAACCGG'],
                          ['ACGTCAGTNNNNNNNN', 'NNNNNNNNTTAACCGG'], fraction=.6)

    def test_fraction_precision_at_offset_selection_boundary(self):
        for fraction in (math.nextafter(.5, 0.), .5, math.nextafter(.5, 1.), .5000000001):
            for k in (0, 1):
                with self.subTest(fraction=fraction, k=k):
                    self.check_counts(['ACGTCAGT', 'TTAACCGG'],
                                      ['ACGTCAGTNNNNNNNN', 'NNNNNNNNTTAACCGG'],
                                      fraction=fraction, k=k)

    def test_no_fallback_when_the_sample_has_no_matches(self):
        self.check_counts(['ACGTCAGT'], ['N' * 8, 'ACGTCAGT'], sample_size=1)

    def test_exact_priority_and_equal_distance_ties(self):
        self.check_counts(['ACGTCAGT', 'ACGTCAGC'], ['ACGTCAGT', 'ACGTCAGC', 'ACGTCAGA'])

    def test_exact_mode_does_not_use_mismatches_for_offset_detection(self):
        self.check_counts(['ACGTCAGT'], ['NNACGTCAGA', 'NNACGTCAGT'], k=0, sample_size=1)

    def test_packed_length_boundaries_and_truncated_reads(self):
        for length in (1, 19, 20, 32):
            target = ('ACGT' * 8)[:length]
            with self.subTest(length=length):
                self.check_counts([target], ['NN' + target, target[:-1] or 'N', 'N' * (length + 2)])

    def test_empty_fastq(self):
        self.check_counts(['ACGTCAGT'], [])

    def test_duplicate_sequences_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            library = folder / 'library.tsv'
            library.write_text('guide\tsequence\tgene\ng0\tACGT\tG0\ng1\tacgt\tG1\n')
            fastq = folder / 'reads.fastq'
            fastq.write_text('@r\nACGT\n+\nIIII\n')
            completed = subprocess.run([str(BIN), 'guide-counter', 'count', '--input', str(fastq),
                                        '--library', str(library), '--output', str(folder / 'out')],
                                       capture_output=True, text=True)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn('unique guide sequences', completed.stderr)
            self.assertFalse((folder / 'out.counts.txt').exists())


if __name__ == '__main__':
    unittest.main()
