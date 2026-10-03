"""Create a small numerical oracle using unmodified, pinned JACKS.

Run in the JACKS comparator environment (NumPy 1.26.4 / SciPy 1.11.4).
The fixture stores input counts as well as upstream full-precision outputs.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np


def generate(source, output):
    commit = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    if commit != 'dd5c4be5e83baa5ee79a589b0a8c3a2ac3f7d6ad':
        raise ValueError('Unexpected JACKS source commit')
    from jacks.preprocess import loadDataAndPreprocess, collateTestControlSamples
    from jacks.infer import inferJACKS, LOG
    import scipy
    import logging
    LOG.setLevel(logging.WARNING)
    rng = np.random.default_rng(20261001)
    names = [f'gene-{i:02}' for i in range(24)]
    controls = names[:6]
    samples = [('baseline-a', 'BASE-A', 'BASE-A'), ('baseline-b1', 'BASE-B', 'BASE-B'),
               ('baseline-b2', 'BASE-B', 'BASE-B'), ('a1', 'A', 'BASE-A'), ('a2', 'A', 'BASE-A'),
               ('b1', 'B', 'BASE-B'), ('b2', 'B', 'BASE-B'), ('b3', 'B', 'BASE-B')]
    rows = []
    for i, gene in enumerate(names):
        effects = np.zeros(2) if gene in controls else rng.normal(-1., 1.2, 2)
        for j, efficacy in enumerate((.4, .8, 1.2, 1.6)):
            abundance = rng.lognormal(7., .35)
            mu = np.array([abundance] * 3 + [abundance * 2**(effects[0] * efficacy)] * 2 +
                          [abundance * 2**(effects[1] * efficacy)] * 3)
            counts = rng.poisson(mu * rng.lognormal(0., .08, len(samples))).tolist()
            if i == 23 and j == 3:
                counts[3] = 0
            rows.append([f'guide-{i:02}-{j}', gene] + counts)
    rows = [rows[i] for i in rng.permutation(len(rows))]
    with tempfile.TemporaryDirectory() as temp:
        counts = Path(temp) / 'counts.tsv'
        with counts.open('w', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t')
            writer.writerow(['sgRNA', 'Gene'] + [r[0] for r in samples])
            writer.writerows(rows)
        gene_spec = {r[0]: r[1] for r in rows}
        ctrl_spec = {r[1]: r[2] for r in samples}
        data, meta, sample_ids, genes, indexes = loadDataAndPreprocess(
            {str(counts): [(r[1], r[0]) for r in samples]}, gene_spec, ctrl_spec,
            ctrl_geneset=set(controls), normtype='ctrl_guides', prior=32)
        test, baseline, test_ids = collateTestControlSamples(data, sample_ids, ctrl_spec)
        results = inferJACKS(indexes, test, baseline)
    expected = {gene: {'effect_estimate': r[4].tolist(), 'effect_sd': np.sqrt(r[5]-r[4]**2).tolist(),
                       'guide_efficacy': r[2].tolist(), 'guide_efficacy_sd': np.sqrt(r[3]-r[2]**2).tolist()}
                for gene, r in results.items()}
    result = {'schema_version': 1, 'source_commit': commit, 'numpy': np.__version__, 'scipy': scipy.__version__,
              'source_hashes': {name: hashlib.sha256((source / name).read_bytes()).hexdigest()
                                for name in ('jacks/jacks/infer.py', 'jacks/jacks/preprocess.py')},
              'seed': 20261001, 'generator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'header': ['sgRNA', 'Gene'] + [r[0] for r in samples], 'rows': rows,
              'sample_map': [['Sample', 'Condition', 'Control']] + [list(r) for r in samples],
              'controls': controls, 'conditions': [sample_ids[i] for i in test_ids], 'expected': expected}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jacks-source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    generate(args.jacks_source, args.out)
