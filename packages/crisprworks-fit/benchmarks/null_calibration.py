"""Known-truth score simulation for within-gene correlated guide artifacts.

This isolates null calibration. It is not an end-to-end sequencing/count-fit
benchmark or evidence of biological accuracy. Both null methods use the same
negative controls and exactly the same declared test-gene family.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from crisprworks_fit.calibration import adjust_fdr
from crisprworks_fit.inference import finite_permutation_tails


def simulate(seed, *, correlated, alternatives, tests=600, controls=4000, guides=4):
    rng = np.random.default_rng(seed)
    artifacts = rng.normal(0, 1. if correlated else 0., size=tests + controls)
    measurements = artifacts[:, None] + rng.normal(0, .35, size=(tests + controls, guides))
    truth = np.zeros(tests, dtype=bool)
    if alternatives:
        truth[:tests // 5] = True
        measurements[:tests // 5] -= 2.5
    observed = measurements[:tests].mean(axis=1)
    control_values = measurements[tests:]
    # Guide resampling breaks the latent within-gene covariance while the
    # whole-gene mean distribution preserves it.
    shuffled = rng.choice(control_values.ravel(), size=(12000, guides)).mean(axis=1)
    gene_p = finite_permutation_tails(control_values.mean(axis=1)[:, None], observed[:, None])[2][:, 0]
    guide_p = finite_permutation_tails(shuffled[:, None], observed[:, None])[2][:, 0]
    output = []
    for method, p in (('pooled_guide_bh', guide_p), ('whole_gene_bh', gene_p), ('whole_gene_by', gene_p)):
        q = adjust_fdr(p, 'by' if method.endswith('_by') else 'bh')
        calls = q <= .05
        false = int((calls & ~truth).sum())
        true = int((calls & truth).sum())
        output.append({'seed': seed, 'correlated': int(correlated), 'alternatives': int(alternatives),
                       'method': method, 'test_genes': tests, 'control_genes': controls, 'guides': guides,
                       'true_calls': true, 'false_calls': false, 'total_calls': int(calls.sum()),
                       'false_discovery_proportion': false / max(1, int(calls.sum())),
                       'power': true / truth.sum() if truth.any() else '',
                       'any_false_call': int(false > 0)})
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=200)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('repeats must be positive')
    rows = []
    for correlated in (False, True):
        for alternatives in (False, True):
            for repeat in range(args.repeats):
                rows.extend(simulate(20000 + repeat, correlated=correlated, alternatives=alternatives))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for correlated in (0, 1):
        for alternatives in (0, 1):
            for method in ('pooled_guide_bh', 'whole_gene_bh', 'whole_gene_by'):
                part = [r for r in rows if r['correlated'] == correlated and r['alternatives'] == alternatives and r['method'] == method]
                print(json.dumps({'correlated': correlated, 'alternatives': alternatives, 'method': method,
                                  'mean_FDP': float(np.mean([r['false_discovery_proportion'] for r in part])),
                                  'any_false_call_probability': float(np.mean([r['any_false_call'] for r in part])),
                                  'mean_power': float(np.mean([r['power'] for r in part])) if alternatives else None}))


if __name__ == '__main__':
    main()
