"""Collect held-out frontier comparisons without fitting on evaluation labels."""

import argparse
import ast
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from crisprworks_fit.calibration import calibrate_genes
from frontier import CELL_LINES, GUIDES
from hit_quality import ranking_metrics


METHODS = ('fit', 'fit_gene_bh', 'fit_gene_by', 'chronos', 'jacks', 'joint', 'joint_gene_bh', 'joint_gene_by')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def collect(panel):
    rows = []
    provenance = json.loads((panel / 'panel.json').read_text())
    provenance['workers'] = []
    provenance['artifacts'] = {}
    worker_source = Path(__file__).with_name('frontier.py').read_text()
    provenance['worker_function_sha256'] = {node.name: hashlib.sha256(
        ast.get_source_segment(worker_source, node).encode()).hexdigest()
        for node in ast.parse(worker_source).body if isinstance(node, ast.FunctionDef) and
        node.name in ('chronos_worker', 'jacks_worker', 'joint_worker')}
    for library in GUIDES:
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
            joint_path = folder / 'joint' / str(fold) / 'screen.joint-details.json'
            joint = json.loads(joint_path.read_text())
            joint_calibration = {method: {(r['gene'], r['condition']): r
                for r in calibrate_genes(joint['genes'], controls, method=method)['rows']}
                for method in ('bh', 'by')}
            provenance['artifacts'][str(joint_path.relative_to(panel))] = digest(joint_path)
            for method in ('chronos', 'jacks', 'joint', 'joint_numpy126'):
                output = folder / method / str(fold)
                provenance['workers'].append({'library': library, 'output_tag': method,
                    **json.loads((output / 'worker.json').read_text()),
                    **json.loads((output / 'command.json').read_text())})
                for path in output.glob('*.csv'):
                    provenance['artifacts'][str(path.relative_to(panel))] = digest(path)
            for cell in CELL_LINES:
                fit_folder = folder / 'fit' / f'{fold}-{cell}'
                details_path = fit_folder / 'screen.fit-details.json'
                genes = json.loads(details_path.read_text())['genes']
                provenance['artifacts'][str(details_path.relative_to(panel))] = digest(details_path)
                provenance['workers'].append({'library': library, 'cell_line': cell, 'method': 'fit', 'fold': fold,
                    **json.loads((fit_folder / 'command.json').read_text()),
                    'manifest': json.loads((fit_folder / 'screen.crisprworks.json').read_text())})
                calibration = {method: {r['gene']: r for r in calibrate_genes(genes, controls, method=method)['rows']}
                               for method in ('bh', 'by')}
                for _, target in heldout.iterrows():
                    gene = target.gene
                    fit = genes[gene]
                    if gene not in chronos.columns or gene not in jacks.columns:
                        raise ValueError('Comparator omitted a held-out reference gene')
                    rows.append({'library': library, 'cell_line': cell, 'gene': gene,
                        'essential': int(target.essential), 'fold': fold,
                        'fit': -fit['beta_estimate'][GUIDES[library]],
                        'fit_q': fit['beta_permute_pval_neg_fdr'][0],
                        'fit_gene_bh': -calibration['bh'][gene]['score'],
                        'fit_gene_bh_q': calibration['bh'][gene]['q_negative'],
                        'fit_gene_by': -calibration['by'][gene]['score'],
                        'fit_gene_by_q': calibration['by'][gene]['q_negative'],
                        'chronos': -float(chronos.loc[cell, gene]),
                        'chronos_q': float(chronos_q.loc[cell, gene]),
                        'jacks': -float(jacks.loc[cell, gene]),
                        'joint': -joint['genes'][gene]['effect_estimate'][joint['conditions'].index(cell)],
                        'joint_gene_bh': -joint['genes'][gene]['effect_estimate'][joint['conditions'].index(cell)],
                        'joint_gene_bh_q': joint_calibration['bh'][gene, joint['conditions'].index(cell)+1]['q_negative'],
                        'joint_gene_by': -joint['genes'][gene]['effect_estimate'][joint['conditions'].index(cell)],
                        'joint_gene_by_q': joint_calibration['by'][gene, joint['conditions'].index(cell)+1]['q_negative']})
    for row in rows:
        if not all(np.isfinite(row[method]) for method in METHODS):
            raise ValueError('Nonfinite comparator score')
    return rows, provenance


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def report(rows, output):
    summaries = []
    paired = []
    for library in GUIDES:
        for cell in CELL_LINES:
            cohort = [r for r in rows if r['library'] == library and r['cell_line'] == cell]
            labels = np.array([int(r['essential']) for r in cohort])
            for method in METHODS:
                scores = np.array([float(r[method]) for r in cohort])
                qkey = method + '_q'
                calls = np.array([float(r[qkey]) <= .05 and float(r[method]) > 0 for r in cohort]) if qkey in cohort[0] else None
                summary = {'library': library, 'cell_line': cell, 'method': method,
                           'reference_genes': len(cohort), 'essentials': int(labels.sum()),
                           'nonessentials': int((labels == 0).sum()), **ranking_metrics(labels, scores),
                           'essential_calls': int(calls[labels == 1].sum()) if calls is not None else '',
                           'nonessential_calls': int(calls[labels == 0].sum()) if calls is not None else '',
                           'essential_call_rate': float(calls[labels == 1].mean()) if calls is not None else '',
                           'nonessential_call_rate': float(calls[labels == 0].mean()) if calls is not None else ''}
                summaries.append(summary)
            # Paired, class-stratified conditional gene bootstrap. Fits and
            # training folds remain frozen; this is not across-screen uncertainty.
            rng = np.random.default_rng(1729)
            positive, negative = np.flatnonzero(labels == 1), np.flatnonzero(labels == 0)
            draws = [np.r_[rng.choice(positive, len(positive)), rng.choice(negative, len(negative))]
                     for _ in range(2000)]
            candidate = np.array([float(r['joint']) for r in cohort])
            for comparator in ('fit', 'chronos'):
                baseline = np.array([float(r[comparator]) for r in cohort])
                deltas = [ranking_metrics(labels[d], candidate[d])['average_precision'] -
                          ranking_metrics(labels[d], baseline[d])['average_precision'] for d in draws]
                paired.append({'library': library, 'cell_line': cell, 'candidate': 'joint',
                               'comparator': comparator, 'bootstraps': 2000, 'seed': 1729,
                               'ap_difference': ranking_metrics(labels, candidate)['average_precision'] -
                                                ranking_metrics(labels, baseline)['average_precision'],
                               'ap_difference_ci_low': float(np.quantile(deltas, .025)),
                               'ap_difference_ci_high': float(np.quantile(deltas, .975))})
    write_csv(output / 'metrics.csv', summaries)
    write_csv(output / 'paired.csv', paired)
    text = ['# Gene-level frontier comparison', '',
            'A frozen panel compares Fit with Chronos 2.3.15 and JACKS 0.2 on three cell lines and two libraries.',
            'Avana (Broad) and KY (Sanger) are curated Chronos vignette subsets enriched for reference controls; these are not full genome-wide screens.',
            'Cell IDs are the first three sorted shared IDs with >=2 replicates: T47D, L363 and NCI-H1915.',
            'Every evaluation gene is excluded from negative-control training in its external fold.', '',
            '| Library | Cell | Method | AP | Essential calls / class size | Nonessential calls / class size |',
            '| --- | --- | --- | --- | --- | --- |']
    for s in summaries:
        calls = '—' if s['essential_calls'] == '' else f"{s['essential_calls']}/{s['essentials']}"
        false = '—' if s['nonessential_calls'] == '' else f"{s['nonessential_calls']}/{s['nonessentials']}"
        text.append(f"| {s['library']} | {s['cell_line']} | {s['method']} | {s['average_precision']:.4f} | {calls} | {false} |")
    text += ['', '## Paired ranking differences', '',
             'Positive AP differences favor production Fit joint. Percentile intervals use 2,000 paired, class-stratified gene bootstrap draws (seed 1729), conditional on frozen fits and folds.',
             'They do not measure uncertainty across independent screens, labels or refitted training sets; the six comparisons are correlated.', '',
             '| Library | Cell | Comparator | AP difference | Conditional 95% interval |',
             '| --- | --- | --- | --- | --- |']
    for p in paired:
        text.append(f"| {p['library']} | {p['cell_line']} | {p['comparator']} | {p['ap_difference']:+.4f} | [{p['ap_difference_ci_low']:+.4f}, {p['ap_difference_ci_high']:+.4f}] |")
    text += ['', '## Methods and limits', '',
             'Fit fits each cell separately; Chronos and JACKS learn guide efficacy jointly across the three cells within each library.',
             'All methods see the same finite-read complete-guide cohort (four guides Avana; five guides KY) and the same training controls.',
             'Chronos trains for a predeclared 1001 epochs with default model regularization and its official empirical p-value / two-stage BH caller.',
             'JACKS retains unmodified upstream inference and control-guide normalization with pseudocount 32; only effect ranking is evaluated here.',
             'Fit joint adapts that tested JACKS model. Every fitted effect is checked against independent unmodified JACKS at 1e-10 absolute/relative tolerance in the same NumPy 1.26.4 environment.',
             'Reported Fit joint runs use the production NumPy 2.3.5 dependency set. NumPy version changes alter variance-window ordering at near-tied means in some guides; cross-version differences are recorded rather than described as exact parity.',
             'Fit uses updated efficiency, control normalization, two pooled-guide permutation rounds, and median/variance conventions inherited from MAGeCK2.',
             'The gene-level Fit callers calibrate those frozen effects using whole control genes within guide-count strata, excluding training controls from the testing family.',
             'Fit MLE families cover one cell per fit; joint families cover all three conditions together. These call counts therefore compare complete analysis choices, not only effect estimators.',
             'BH assumes independent or suitable positively dependent valid p-values; BY tolerates arbitrary dependence of valid p-values, at a substantial power cost.',
             'Calls use negative effects and directional q<=0.05. Nonessential call rates are class-specific false-positive frequencies, not measured genome-wide FDR.',
             'Fit gene callers have identical AP to Fit because calibration does not refit effects. No hyperparameters were selected on the evaluation annotations.',
             '', 'This panel does not establish overall SOTA: it is enriched for controls, has only three cell lines, has no copy-number correction, and does not independently validate biological truth.',
             'The annotation collections and source screen populations may overlap; gene holdout cannot remove that collection-level bias.', '']
    (output / 'FRONTIER.md').write_text('\n'.join(text))
    return summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--panel', type=Path)
    parser.add_argument('--scores', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.scores:
        rows = list(csv.DictReader(args.scores.open()))
    else:
        if args.panel is None:
            parser.error('panel or scores is required')
        rows, provenance = collect(args.panel)
        write_csv(args.out / 'scores.csv', rows)
        provenance['scores_sha256'] = digest(args.out / 'scores.csv')
        paths = [Path(__file__), Path(__file__).with_name('frontier.py'),
                 Path(__file__).parents[1] / 'src' / 'crisprworks_fit' / 'calibration.py',
                 Path(__file__).parents[1] / 'src' / 'crisprworks_fit' / 'joint.py',
                 Path(__file__).parents[1] / 'src' / 'crisprworks_fit' / '_joint_kernel.py']
        provenance['implementation_hashes'] = {str(p.relative_to(Path(__file__).parents[1])): digest(p) for p in paths}
        (args.out / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    for row in report(rows, args.out):
        print(f"{row['library']} {row['cell_line']} {row['method']} AP={row['average_precision']:.4f} TP={row['essential_calls']} FP={row['nonessential_calls']}")


if __name__ == '__main__':
    main()
