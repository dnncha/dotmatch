"""Frozen Broad/Sanger panel with current Chronos and unmodified JACKS.

Comparator workers run in separate environments. The three cell lines are
the first sorted shared IDs with >=2 replicates, selected without looking at
effects. This curated vignette panel is not a genome-wide accuracy estimate.
"""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

CHRONOS_COMMIT = "d85e55f4caabaea240d58c8355894613e5410b5a"
JACKS_COMMIT = "dd5c4be5e83baa5ee79a589b0a8c3a2ac3f7d6ad"
CELL_LINES = ("ACH-000147", "ACH-000183", "ACH-000434")
GUIDES = {"Avana": 4, "KY": 5}
SOURCE_HASHES = {
    "AchillesCommonEssentialControls.csv": "494ab3798c083bb6dd77d0f82b9f6ba4a26261df84704dd117847e3a9533ec51",
    "AchillesNonessentialControls.csv": "662e90c7d3f7ee02bdd64a0a517845281f37935e35d7934251e9a9e95c3cb7cd",
    "AvanaReadcounts.hdf5": "185808e90d4a7690e5ac9092cbad516b657831cadfc2a61c54516800d92744fa",
    "AvanaGuideMap.csv": "50d58da481f4f1634c06e89ae189e8fe7bfb4ea3b0a20de3eac8b6c4a6cdfe3f",
    "AvanaSequenceMap.csv": "1f9caefbe998822d3decf114dd48427dc5d34a2de1a3dea102635c00c2b4cef1",
    "KYReadcounts.hdf5": "ea287dbe44678018d02cca71ce14df9a17f322234618596f74bac50ee57835bc",
    "KYGuideMap.csv": "23d942b8af1542eca28759ecebac6d15dcd93876e9c3f0589bd647356abafa32",
    "KYSequenceMap.csv": "bf0df34353db70fdfa827ebedeffad1ea6fe3b6443617c236a0aeeadb3f0b606",
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_hdf(path):
    import h5py
    import pandas as pd
    with h5py.File(path) as stream:
        return pd.DataFrame(stream['data'][:],
                            index=[s.decode() for s in stream['dim_0'][:]],
                            columns=[s.decode() for s in stream['dim_1'][:]])


def prepare(source, output):
    import pandas as pd
    from hit_quality import heldout_folds
    data_dir = source / 'Data' / 'SampleData'
    expected_commit = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    if expected_commit != CHRONOS_COMMIT:
        raise ValueError('Chronos source commit differs from the recorded version')
    for name, expected in SOURCE_HASHES.items():
        if digest(data_dir / name) != expected:
            raise ValueError('Frozen source data changed: ' + name)
    provenance = {'schema_version': 1, 'chronos_commit': CHRONOS_COMMIT, 'jacks_commit': JACKS_COMMIT,
                  'cell_lines': CELL_LINES, 'seed': 1729, 'folds': 5, 'chronos_epochs': 1001,
                  'selection': 'first three sorted shared cell IDs with >=2 replicates; complete modal-guide-count genes and finite shared reads',
                  'scope': 'curated vignette subset enriched for control annotations; not full-library or unbiased genome-wide accuracy',
                  'sources': {}, 'libraries': {}}
    classes = {}
    for name, value in [('AchillesCommonEssentialControls.csv', 1), ('AchillesNonessentialControls.csv', 0)]:
        path = data_dir / name
        provenance['sources'][name] = digest(path)
        for gene in pd.read_csv(path).Gene:
            if gene in classes and classes[gene] != value:
                raise ValueError('Conflicting reference annotation')
            classes[gene] = value
    maps = {lib: pd.read_csv(data_dir / (lib + 'SequenceMap.csv')) for lib in GUIDES}
    eligible = {lib: set(m[m.cell_line_name != 'pDNA'].groupby('cell_line_name').size().loc[lambda n: n >= 2].index)
                for lib, m in maps.items()}
    if tuple(sorted(eligible['Avana'] & eligible['KY'])[:3]) != CELL_LINES:
        raise ValueError('Cell selection changed')
    for library, guide_count in GUIDES.items():
        folder = output / library
        folder.mkdir(parents=True, exist_ok=True)
        for suffix in ('Readcounts.hdf5', 'GuideMap.csv', 'SequenceMap.csv'):
            name = library + suffix
            provenance['sources'][name] = digest(data_dir / name)
        counts = read_hdf(data_dir / (library + 'Readcounts.hdf5'))
        mapping = pd.read_csv(data_dir / (library + 'GuideMap.csv'))[['sgrna', 'gene']]
        seq = maps[library]
        treatments = seq[seq.cell_line_name.isin(CELL_LINES)]
        baselines = seq[(seq.cell_line_name == 'pDNA') & seq.pDNA_batch.isin(treatments.pDNA_batch)]
        selected_seq = pd.concat([baselines, treatments]).reset_index(drop=True)
        sizes = mapping.groupby('gene').size()
        mapping = mapping[mapping.gene.isin(sizes[sizes == guide_count].index)]
        counts = counts.loc[selected_seq.sequence_ID, mapping.sgrna]
        finite_guides = counts.columns[np.isfinite(counts.values).all(axis=0)]
        mapping = mapping[mapping.sgrna.isin(finite_guides)]
        remaining = mapping.groupby('gene').size()
        mapping = mapping[mapping.gene.isin(remaining[remaining == guide_count].index)]
        counts = counts[mapping.sgrna]
        if not (counts.values >= 0).all():
            raise ValueError('Negative raw count')
        genes = set(mapping.gene)
        folds = heldout_folds([g for g in genes if classes.get(g) == 1],
                              [g for g in genes if classes.get(g) == 0])
        mapping.to_csv(folder / 'guide_map.csv', index=False)
        selected_seq.to_csv(folder / 'sequence_map.csv', index=False)
        counts.to_csv(folder / 'readcounts.csv')
        table = mapping.rename(columns={'sgrna': 'sgRNA', 'gene': 'Gene'}).set_index('sgRNA')
        table = table.join(counts.T)
        table.to_csv(folder / 'counts.tsv', sep='\t')
        pd.DataFrame([{'gene': gene, 'essential': classes[gene], 'fold': folds[gene]} for gene in sorted(folds)]).to_csv(folder / 'labels.csv', index=False)
        provenance['libraries'][library] = {'genes': len(genes), 'guides': len(mapping), 'guide_count': guide_count,
                                            'evaluated_genes': len(folds), 'counts_sha256': digest(folder / 'counts.tsv')}
        for fold in range(5):
            controls = sorted(g for g in folds if folds[g] != fold and classes[g] == 0)
            (folder / f'controls-{fold}.txt').write_text('\n'.join(controls) + '\n')
        for cell in CELL_LINES:
            cell_map = treatments[treatments.cell_line_name == cell]
            batches = set(cell_map.pDNA_batch)
            matched = baselines[baselines.pDNA_batch.isin(batches)]
            if len(batches) != 1 or len(matched) != 1:
                raise ValueError('Panel requires one matched pDNA sample per cell')
            samples = matched.sequence_ID.tolist() + cell_map.sequence_ID.tolist()
            table[['Gene'] + samples].to_csv(folder / f'counts-{cell}.tsv', sep='\t')
            (folder / f'design-{cell}.tsv').write_text('Samples\tbaseline\ttreatment\n' + ''.join(
                f'{sample}\t1\t{int(i > 0)}\n' for i, sample in enumerate(samples)))
    (output / 'panel.json').write_text(json.dumps(provenance, indent=2) + '\n')


def chronos_worker(folder, fold, output):
    import pandas as pd
    import tensorflow as tf
    from chronos import Chronos
    from chronos.hit_calling import get_pvalue_dependent, get_fdr_from_pvalues
    tf.compat.v1.set_random_seed(1729)
    np.random.seed(1729)
    counts = pd.read_csv(folder / 'readcounts.csv', index_col=0)
    mapping = pd.read_csv(folder / 'guide_map.csv')
    seq = pd.read_csv(folder / 'sequence_map.csv')
    controls = (folder / f'controls-{fold}.txt').read_text().splitlines()
    guides = mapping.loc[mapping.gene.isin(controls), 'sgrna'].tolist()
    model = Chronos(readcounts={'panel': counts}, guide_gene_map={'panel': mapping},
                   sequence_map={'panel': seq}, negative_control_sgrnas={'panel': guides},
                   print_to=str(output / 'training.log'))
    initial_cost = float(model.cost)
    model.train(nepochs=1001)
    effects = model.gene_effect
    if not np.isfinite(effects.values).all():
        raise ValueError('Nonfinite Chronos effect')
    effects.to_csv(output / 'effects.csv')
    p = get_pvalue_dependent(effects, negative_controls=controls)
    p.to_csv(output / 'p.csv')
    get_fdr_from_pvalues(p).to_csv(output / 'q.csv')
    from importlib.metadata import version
    (output / 'worker.json').write_text(json.dumps({'method': 'chronos', 'fold': fold,
        'epochs': model.epoch, 'initial_cost': initial_cost, 'final_cost': float(model.cost),
        'version': version('crispr_chronos'), 'tensorflow': tf.__version__, 'numpy': np.__version__,
        'inference': 'official empirical_pvalue and default FDR_TSBH', 'negative_control_genes': len(controls)}, indent=2) + '\n')
    model.sess.close()


def jacks_worker(folder, fold, output):
    import logging
    from jacks.preprocess import loadDataAndPreprocess, collateTestControlSamples
    from jacks.infer import inferJACKS, LOG
    import pandas as pd
    LOG.setLevel(logging.WARNING)
    seq = pd.read_csv(folder / 'sequence_map.csv')
    mapping = pd.read_csv(folder / 'guide_map.csv')
    controls = set((folder / f'controls-{fold}.txt').read_text().splitlines())
    sample_spec = []
    ctrl_spec = {}
    for cell in CELL_LINES:
        cm = seq[seq.cell_line_name == cell]
        batch = cm.pDNA_batch.iloc[0]
        baseline = seq[(seq.cell_line_name == 'pDNA') & (seq.pDNA_batch == batch)].sequence_ID.iloc[0]
        control_id = 'CONTROL-' + str(batch)
        if control_id not in ctrl_spec:
            sample_spec.append((control_id, baseline))
            ctrl_spec[control_id] = control_id
        sample_spec.extend((cell, sample) for sample in cm.sequence_ID)
        ctrl_spec[cell] = control_id
    gene_spec = dict(zip(mapping.sgrna, mapping.gene))
    data, meta, sample_ids, genes, gene_index = loadDataAndPreprocess(
        {str(folder / 'counts.tsv'): sample_spec}, gene_spec, ctrl_spec=ctrl_spec,
        normtype='ctrl_guides', ctrl_geneset=controls, prior=32)
    testdata, ctrldata, test_idxs = collateTestControlSamples(data, sample_ids, ctrl_spec)
    results = inferJACKS(gene_index, testdata, ctrldata)
    effects = pd.DataFrame({gene: result[4] for gene, result in results.items()},
                           index=[sample_ids[i] for i in test_idxs])
    if not np.isfinite(effects.values).all():
        raise ValueError('Nonfinite JACKS effect')
    effects.to_csv(output / 'effects.csv')
    (output / 'worker.json').write_text(json.dumps({'method': 'jacks', 'fold': fold,
        'numpy': np.__version__, 'max_iterations': 50, 'count_prior': 32,
        'normalization': 'control-guide median log counts', 'negative_control_genes': len(controls),
        'scope': 'joint three-cell effects; no pretrained efficacy reference or pseudogenes'}, indent=2) + '\n')


def run_job(command, output):
    output.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    with (output / 'stdout.log').open('w') as stream:
        completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, timeout=1800)
    (output / 'command.json').write_text(json.dumps({'command': command, 'returncode': completed.returncode,
                                                   'wall_seconds': time.perf_counter() - start}, indent=2) + '\n')
    if completed.returncode:
        raise RuntimeError(f'Comparator failed; see {output / "stdout.log"}')


def joint_worker(folder, fold, output):
    import pandas as pd
    from crisprworks_fit.joint import fit_joint
    seq = pd.read_csv(folder / 'sequence_map.csv')
    samples, seen = [], set()
    for cell in CELL_LINES:
        cm = seq[seq.cell_line_name == cell]
        batch = cm.pDNA_batch.iloc[0]
        baseline = seq[(seq.cell_line_name == 'pDNA') & (seq.pDNA_batch == batch)].sequence_ID.iloc[0]
        control = 'CONTROL-' + str(batch)
        if control not in seen:
            samples.append([baseline, control, control])
            seen.add(control)
        samples.extend([sample, cell, control] for sample in cm.sequence_ID)
    sample_map = output / 'sample-map.tsv'
    with sample_map.open('w', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['Sample', 'Condition', 'Control'])
        writer.writerows(samples)
    controls = (folder / f'controls-{fold}.txt').read_text().splitlines()
    result = fit_joint(folder / 'counts.tsv', sample_map, controls)
    (output / 'screen.joint-details.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    effects = pd.DataFrame({g: v['effect_estimate'] for g, v in result['genes'].items()}, index=result['conditions'])
    effects.to_csv(output / 'effects.csv')
    oracle = pd.read_csv(folder / 'jacks' / str(fold) / 'effects.csv', index_col=0).loc[effects.index, effects.columns]
    oracle_numpy = json.loads((folder / 'jacks' / str(fold) / 'worker.json').read_text())['numpy']
    same_numpy = np.__version__ == oracle_numpy
    if same_numpy:
        np.testing.assert_allclose(effects.values, oracle.values, rtol=1e-10, atol=1e-10)
    (output / 'worker.json').write_text(json.dumps({'method': 'joint', 'fold': fold, 'numpy': np.__version__,
        'oracle': 'unmodified JACKS effects, same frozen counts and controls',
        'oracle_numpy': oracle_numpy, 'same_numpy_strict_parity': same_numpy,
        'max_abs_effect_error': float(np.max(np.abs(effects.values-oracle.values))),
        'source_commit': JACKS_COMMIT}, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--method', choices=('prepare', 'chronos', 'jacks', 'fit', 'joint'), required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--panel', type=Path, required=True)
    parser.add_argument('--library', choices=tuple(GUIDES))
    parser.add_argument('--fold', type=int, choices=range(5))
    parser.add_argument('--out', type=Path)
    parser.add_argument('--python', type=Path, help='Run all jobs of this method with this worker Python')
    parser.add_argument('--output-tag', choices=('joint', 'joint_numpy126'), help='Separate same-NumPy oracle validation from production joint results')
    args = parser.parse_args()
    for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'TF_NUM_INTRAOP_THREADS', 'TF_NUM_INTEROP_THREADS'):
        os.environ[name] = '1'
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
    if args.method == 'prepare':
        if not args.source:
            parser.error('prepare requires --source')
        args.panel.mkdir(parents=True, exist_ok=True)
        prepare(args.source.resolve(), args.panel.resolve())
    elif args.python:
        for library in GUIDES:
            for fold in range(5):
                folder = args.panel / library
                if args.method == 'fit':
                    for cell in CELL_LINES:
                        output = folder / 'fit' / f'{fold}-{cell}'
                        command = [str(args.python), '-m', 'crisprworks_fit', 'mle',
                                   '-k', str(folder / f'counts-{cell}.tsv'), '-d', str(folder / f'design-{cell}.tsv'),
                                   '-n', str(output / 'screen'), '--kernel', 'native', '--norm-method', 'control',
                                   '--control-gene', str(folder / f'controls-{fold}.txt'), '--permutation-pvalues', 'finite',
                                   '--seed', '1729', '--permutation-round', '2', '--genes-varmodeling', '1000',
                                   '--update-efficiency', '--write-fit-details']
                        run_job(command, output)
                        print(f'fit {library} fold {fold} {cell}', flush=True)
                else:
                    output = folder / (args.output_tag or args.method) / str(fold)
                    command = [str(args.python), str(Path(__file__).resolve()), '--method', args.method,
                               '--panel', str(args.panel.resolve()), '--library', library, '--fold', str(fold), '--out', str(output)]
                    run_job(command, output)
                    print(f'{args.method} {library} fold {fold}', flush=True)
    else:
        if args.library is None or args.fold is None or args.out is None:
            parser.error('worker requires library, fold, out')
        args.out.mkdir(parents=True, exist_ok=True)
        function = {'chronos': chronos_worker, 'jacks': jacks_worker, 'joint': joint_worker}.get(args.method)
        if function is None:
            parser.error('fit requires --python orchestration')
        function(args.panel / args.library, args.fold, args.out)


if __name__ == '__main__':
    main()
