"""Regenerate score-level null-calibration evidence from replicate-level CSV."""

import argparse
import csv
from pathlib import Path

import numpy as np


def report(source, output):
    with source.open() as stream:
        rows = list(csv.DictReader(stream))
    summaries = []
    text = ['# Correlated-guide null calibration', '',
            'This known-truth simulation isolates calibration of frozen scores. It does not simulate read counts, fit MAGeCK/JACKS models or measure biological accuracy.',
            'Each of 200 seeds supplies 600 test genes, 4,000 independently sampled control genes and four guides per gene.',
            'Guide residual SD is 0.35. The correlated scenario adds one shared N(0,1) artifact to all guides of each gene.',
            'Mixed scenarios shift 20% of test genes by -2.5. All callers use the same controls and test family at directional q<=0.05.',
            'The pooled null resamples individual control guides into 12,000 four-guide averages; the gene null retains the 4,000 control gene averages.',
            'The reported mean false discovery proportion estimates FDR in this specified simulation. Intervals are normal Monte Carlo intervals for the replicate mean, not a guarantee across unknown biological scenarios.', '',
            '| Shared artifact | True effects | Method | Mean FDP [MC 95% interval] | Mean power | P(any false call) |',
            '| --- | --- | --- | --- | --- | --- |']
    for correlated in (0, 1):
        for alternatives in (0, 1):
            for method in ('pooled_guide_bh', 'whole_gene_bh', 'whole_gene_by'):
                part = [r for r in rows if int(r['correlated']) == correlated and
                        int(r['alternatives']) == alternatives and r['method'] == method]
                if len(part) < 2:
                    raise ValueError('At least two replicates per scenario are required')
                fdp = np.array([float(r['false_discovery_proportion']) for r in part])
                mean = float(fdp.mean())
                width = 1.96 * float(fdp.std(ddof=1)) / len(part)**.5
                power = float(np.mean([float(r['power']) for r in part])) if alternatives else ''
                any_false = float(np.mean([int(r['any_false_call']) for r in part]))
                s = {'correlated': correlated, 'alternatives': alternatives, 'method': method,
                     'repeats': len(part), 'mean_fdp': mean, 'mc_ci_low': max(0., mean-width),
                     'mc_ci_high': min(1., mean+width), 'mean_power': power, 'any_false_probability': any_false}
                summaries.append(s)
                power_text = '—' if power == '' else f'{power:.3f}'
                text.append(f"| {bool(correlated)} | {bool(alternatives)} | {method} | {mean:.3f} [{s['mc_ci_low']:.3f}, {s['mc_ci_high']:.3f}] | {power_text} | {any_false:.3f} |")
    text += ['', 'Under shared artifacts with true effects, pooled guide averaging understates null variance: mean FDP rises above the nominal 0.05 level.',
             'The whole-gene null preserves that correlation and reduces false discoveries, with a substantial power cost. BY is more conservative still.',
             'Under null-only scenarios, finite p-value resolution can prevent all calls; zero observed false calls does not prove zero population risk.',
             'A shared finite control set induces dependence between empirical test p-values. BH validity requires appropriate independence or positive dependence; this simulation alone does not establish that for all assays.',
             'BY accommodates arbitrary dependence of valid marginal p-values. Both methods still require exchangeable, representative control genes within guide-count strata.', '']
    output.mkdir(parents=True, exist_ok=True)
    with (output / 'null_metrics.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    (output / 'NULL-CALIBRATION.md').write_text('\n'.join(text))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scores', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    report(args.scores, args.out)
