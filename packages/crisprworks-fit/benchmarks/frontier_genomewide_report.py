import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from crisprworks_fit.calibration import calibrate_genes
from frontier_genomewide import LIBRARIES, panel_cells
from hit_quality import ranking_metrics

DESCRIPTION = 'Score the genome-wide DepMap frontier panel and summarise across independent cell lines.'

METHODS = ('fit', 'fit_gene_bh', 'chronos', 'jacks', 'joint', 'joint_gene_bh')
CALLERS = {'fit': 'fit_q', 'fit_gene_bh': 'fit_gene_bh_q', 'chronos': 'chronos_q', 'joint_gene_bh': 'joint_gene_bh_q'}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def converged_family(genes, controls):
    unconverged = sorted(g for g, v in genes.items()
                         if 'fit_diagnostics' in v and not v['fit_diagnostics']['converged'])
    kept = {g: v for g, v in genes.items() if g not in unconverged}
    return kept, [c for c in controls if c not in unconverged], unconverged


def collect(panel):
    cells = panel_cells(panel)
    rows = []
    collect.unconverged = []
    for library, guide_count in LIBRARIES.items():
        folder = panel / library
        labels = pd.read_csv(folder / 'labels.csv')
        for fold in range(5):
            controls = (folder / f'controls-{fold}.txt').read_text().splitlines()
            heldout = labels[labels.fold == fold]
            if set(heldout.gene) & set(controls):
                raise ValueError('Held-out genes leaked into training controls')
            chronos = pd.read_csv(folder / 'chronos' / str(fold) / 'effects.csv', index_col=0)
            chronos_q = pd.read_csv(folder / 'chronos' / str(fold) / 'q.csv', index_col=0)
            jacks = pd.read_csv(folder / 'jacks' / str(fold) / 'effects.csv', index_col=0)
            joint = json.loads((folder / 'joint' / str(fold) / 'screen.joint-details.json').read_text())
            family, family_controls, unconverged = converged_family(joint['genes'], controls)
            collect.unconverged += [{'library': library, 'fold': fold, 'gene': g} for g in unconverged]
            joint_bh = {(r['gene'], r['condition']): r
                        for r in calibrate_genes(family, family_controls, method='bh')['rows']}
            for cell in cells:
                genes = json.loads((folder / 'fit' / f'{fold}-{cell}' / 'screen.fit-details.json').read_text())['genes']
                fit_bh = {r['gene']: r for r in calibrate_genes(genes, controls, method='bh')['rows']}
                column = joint['conditions'].index(cell)
                for target in heldout.itertuples():
                    gene = target.gene
                    joint_effect = -joint['genes'][gene]['effect_estimate'][column]
                    rows.append({'library': library, 'cell_line': cell, 'gene': gene,
                                 'essential': int(target.essential), 'fold': fold,
                                 'fit': -genes[gene]['beta_estimate'][guide_count],
                                 'fit_q': genes[gene]['beta_permute_pval_neg_fdr'][0],
                                 'fit_gene_bh': -fit_bh[gene]['score'],
                                 'fit_gene_bh_q': fit_bh[gene]['q_negative'],
                                 'chronos': -float(chronos.loc[cell, gene]),
                                 'chronos_q': float(chronos_q.loc[cell, gene]),
                                 'jacks': -float(jacks.loc[cell, gene]),
                                 'joint': joint_effect,
                                 'joint_gene_bh': joint_effect,
                                 'joint_gene_bh_q': joint_bh[gene, column + 1]['q_negative']
                                 if (gene, column + 1) in joint_bh else 1.0})
    for row in rows:
        if not all(np.isfinite(row[m]) for m in METHODS):
            raise ValueError(f'Nonfinite score: {row}')
    return rows


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def per_cell(rows):
    frame = pd.DataFrame(rows)
    out = []
    for (library, cell), cohort in frame.groupby(['library', 'cell_line']):
        labels = cohort.essential.to_numpy(int)
        for method in METHODS:
            scores = cohort[method].to_numpy(float)
            row = {'library': library, 'cell_line': cell, 'method': method,
                   'essentials': int(labels.sum()), 'nonessentials': int((labels == 0).sum()),
                   **ranking_metrics(labels, scores)}
            if method in CALLERS:
                calls = (cohort[CALLERS[method]].to_numpy(float) <= .05) & (scores > 0)
                row['essential_call_rate'] = float(calls[labels == 1].mean())
                row['nonessential_call_rate'] = float(calls[labels == 0].mean())
            out.append(row)
    return pd.DataFrame(out)


def cell_bootstrap(values, draws=10000, seed=1729):
    rng = np.random.default_rng(seed)
    means = rng.choice(values, (draws, len(values))).mean(axis=1)
    return float(np.quantile(means, .025)), float(np.quantile(means, .975))


def summarise(metrics, panel_info, output):
    ap = metrics.pivot_table(index=['library', 'cell_line'], columns='method', values='average_precision')
    lines = ['# Genome-wide frontier: DepMap 24Q4', '',
             f"{len(panel_info['cell_lines'])} cell lines screened in both Avana and KY, selected as the first sorted eligible ModelIDs without looking at effects "
             f"({panel_info['eligible_shared_cells']} were eligible).",
             'Every UsedByChronos single-alignment gene with the modal guide count is fitted: ' +
             ', '.join(f"{lib} {v['genes']} genes / {v['guides']} guides / {v['samples']} samples" for lib, v in panel_info['libraries'].items()) + '.',
             'Evaluation genes are Achilles common-essential and nonessential controls, held out from negative-control training in five folds.', '',
             '## Average precision across cell lines', '',
             '| Library | Method | Mean AP | Median AP | Min AP | Essential call rate | Nonessential call rate |',
             '| --- | --- | --- | --- | --- | --- | --- |']
    for library in LIBRARIES:
        sub = metrics[metrics.library == library]
        for method in METHODS:
            m = sub[sub.method == method]
            tp = f"{m.essential_call_rate.mean():.3f}" if method in CALLERS else '—'
            fp = f"{m.nonessential_call_rate.mean():.3f}" if method in CALLERS else '—'
            lines.append(f"| {library} | {method} | {m.average_precision.mean():.4f} | {m.average_precision.median():.4f} | "
                         f"{m.average_precision.min():.4f} | {tp} | {fp} |")
    lines += ['', '## Paired differences across independent cell lines', '',
              'Each cell line is one paired observation; intervals are 10,000 cell-level bootstrap draws (seed 1729) of the mean AP difference.', '',
              '| Library | Candidate | Comparator | Mean AP difference | 95% interval | Cells where candidate is higher |',
              '| --- | --- | --- | --- | --- | --- |']
    paired = []
    for library in LIBRARIES:
        table = ap.loc[library]
        for candidate, comparator in (('joint', 'chronos'), ('fit', 'chronos'), ('joint', 'fit'), ('jacks', 'chronos')):
            delta = (table[candidate] - table[comparator]).to_numpy()
            low, high = cell_bootstrap(delta)
            wins = int((delta > 0).sum())
            paired.append({'library': library, 'candidate': candidate, 'comparator': comparator,
                           'mean_ap_difference': float(delta.mean()), 'ci_low': low, 'ci_high': high,
                           'cells_higher': wins, 'cells': len(delta)})
            lines.append(f"| {library} | {candidate} | {comparator} | {delta.mean():+.4f} | [{low:+.4f}, {high:+.4f}] | {wins}/{len(delta)} |")
    lines += ['', '## Methods and limits', '',
              'Chronos (commit d85e55f, TensorFlow 2.20) trains 1001 epochs with default regularization; its calls use the official empirical p-value and two-stage BH.',
              'Chronos did not use the GPU meaningfully at this scale (0% utilisation, one CPU thread, about 29 minutes per fold); two folds ran with an L4 visible and eight on CPU only.',
              'JACKS (commit dd5c4be) and Fit joint learn guide efficacy jointly across all panel cell lines in each library; Fit MLE fits each cell separately.',
              'Fit gene BH calibrates frozen Fit effects with whole control genes, excluding training controls from the testing family.',
              'Call rates use negative effects and directional q<=0.05; nonessential call rates are class-specific false-positive frequencies, not genome-wide FDR.',
              'No copy-number correction is applied to any method. The Achilles reference sets derive partly from DepMap screens, so absolute AP is optimistic for every method; paired differences are the quantity of interest.',
              'No hyperparameters were selected on evaluation annotations.',
              f"Fit joint hit its 50-iteration limit for {len(panel_info['unconverged_joint_genes'])} gene-fold fits "
              f"({', '.join(sorted({u['gene'] for u in panel_info['unconverged_joint_genes']})) or 'none'}); "
              'those genes are left out of the joint gene-BH testing family and are never called, but keep their effects for ranking.', '']
    (output / 'FRONTIER-GENOMEWIDE.md').write_text('\n'.join(lines))
    write_csv(output / 'paired.csv', paired)
    return lines


def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--panel', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = collect(args.panel)
    write_csv(args.out / 'scores.csv', rows)
    metrics = per_cell(rows)
    metrics.to_csv(args.out / 'metrics.csv', index=False)
    panel_info = json.loads((args.panel / 'panel.json').read_text())
    panel_info['scores_sha256'] = digest(args.out / 'scores.csv')
    panel_info['unconverged_joint_genes'] = collect.unconverged
    (args.out / 'provenance.json').write_text(json.dumps(panel_info, indent=2) + '\n')
    print('\n'.join(summarise(metrics, panel_info, args.out)))


if __name__ == '__main__':
    main()
