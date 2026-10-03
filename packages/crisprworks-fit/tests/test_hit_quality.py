"""Benchmark metrics must handle ties and exclude test genes from training."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np

spec = importlib.util.spec_from_file_location(
    'hit_quality', Path(__file__).parents[1] / 'benchmarks' / 'hit_quality.py')
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


class HitQualityTests(unittest.TestCase):
    def test_known_rankings_and_ties(self):
        labels = [1, 1, 0, 0]
        perfect = benchmark.ranking_metrics(labels, [4, 3, 2, 1])
        self.assertEqual(perfect['average_precision'], 1.)
        self.assertEqual(perfect['roc_auc'], 1.)
        self.assertEqual(perfect['recall_at_95pct_precision'], 1.)
        tied = benchmark.ranking_metrics(labels, [1, 1, 1, 1])
        self.assertEqual(tied['average_precision'], .5)
        self.assertEqual(tied['roc_auc'], .5)
        self.assertEqual(tied['recall_at_95pct_precision'], 0.)
        reversed_result = benchmark.ranking_metrics(labels, [1, 2, 3, 4])
        self.assertAlmostEqual(reversed_result['average_precision'], (1 / 3 + 1 / 2) / 2)
        self.assertEqual(reversed_result['roc_auc'], 0.)
        with self.assertRaises(ValueError):
            benchmark.ranking_metrics(labels, [1, 2, np.nan, 4])

    def test_metrics_do_not_depend_on_order_within_ties(self):
        rng = np.random.default_rng(44)
        labels = rng.integers(0, 2, 100)
        scores = rng.integers(-3, 4, 100)
        expected = benchmark.ranking_metrics(labels, scores)
        for _ in range(20):
            order = rng.permutation(100)
            self.assertEqual(benchmark.ranking_metrics(labels[order], scores[order]), expected)

    def test_folds_are_disjoint_stratified_and_reproducible(self):
        essentials = [f'E{i}' for i in range(13)]
        nonessentials = [f'N{i}' for i in range(17)]
        folds = benchmark.heldout_folds(essentials, nonessentials)
        self.assertEqual(folds, benchmark.heldout_folds(essentials[::-1], nonessentials[::-1]))
        self.assertEqual(set(folds), set(essentials + nonessentials))
        for fold in range(5):
            evaluation = {gene for gene in folds if folds[gene] == fold}
            training = {gene for gene in folds if folds[gene] != fold}
            self.assertFalse(evaluation & training)
            self.assertTrue(evaluation & set(essentials))
            self.assertTrue(evaluation & set(nonessentials))
        with self.assertRaises(ValueError):
            benchmark.heldout_folds(essentials, essentials)


if __name__ == '__main__':
    unittest.main()
