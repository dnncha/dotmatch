"""Plot paired ranking differences and conditional bootstrap intervals."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def plot(source, output):
    with source.open() as stream:
        rows = list(csv.DictReader(stream))
    cells = {'ACH-000147': 'T47D', 'ACH-000183': 'L363', 'ACH-000434': 'NCI-H1915'}
    groups = [(library, cell) for library in ('Avana', 'KY') for cell in cells]
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'svg.hashsalt': 'crisprworks-fit-frontier'})
    fig, ax = plt.subplots(figsize=(8.3, 5.2))
    for method, offset, color, label in [('fit', -.13, '#2563a0', 'Joint minus Fit MLE'),
                                         ('chronos', .13, '#ba6326', 'Joint minus Chronos')]:
        subset = [next(r for r in rows if r['library'] == library and r['cell_line'] == cell and
                       r['comparator'] == method) for library, cell in groups]
        means = np.array([float(r['ap_difference']) for r in subset])
        low = np.array([float(r['ap_difference_ci_low']) for r in subset])
        high = np.array([float(r['ap_difference_ci_high']) for r in subset])
        ax.errorbar(means, np.arange(len(groups)) + offset, xerr=[means-low, high-means], fmt='o',
                    markersize=5, capsize=3, linewidth=1.5, color=color, label=label)
    ax.axvline(0, color='#64748b', linestyle='--', linewidth=1)
    ax.axhline(2.5, color='#e2e8f0', linewidth=1)
    ax.set_yticks(range(len(groups)), [f'{lib} / {cells[cell]}' for lib, cell in groups])
    ax.invert_yaxis()
    ax.set_xlabel('Average precision difference (positive favors Fit joint)')
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.grid(axis='x', color='#e2e8f0', linewidth=.6)
    ax.set_axisbelow(True)
    fig.legend(*ax.get_legend_handles_labels(), loc='upper left', bbox_to_anchor=(.04, .87),
               frameon=False, fontsize=9, ncol=2)
    fig.suptitle('Joint guide-efficacy inference: paired ranking comparison', fontsize=13, x=.04, ha='left', y=.98)
    fig.text(.04, .91, 'Three cell lines × two curated library subsets; 2,000 paired gene bootstrap draws', fontsize=9, color='#475569')
    fig.text(.04, .02, '95% intervals are conditional on frozen fits and reference annotations; comparisons share genes and cells.',
             fontsize=8, color='#475569')
    fig.tight_layout(rect=(0, .06, 1, .89))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, metadata={'Date': None})
    fig.savefig(output.with_suffix('.pdf'), metadata={'CreationDate': None, 'ModDate': None})
    fig.savefig(output.with_suffix('.png'), dpi=140)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--paired', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    plot(args.paired, args.out)
