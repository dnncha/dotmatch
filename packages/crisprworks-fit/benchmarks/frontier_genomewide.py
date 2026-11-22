import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import frontier

DESCRIPTION = "Genome-wide DepMap 24Q4 frontier: Fit, Fit joint, Chronos and JACKS across many cell lines. Uses the public DepMap 24Q4 Figshare release (article 27993248). Cells are the first N sorted ModelIDs screened in both Avana and KY with >=2 passing replicates and one pDNA batch, selected without looking at effects. Every UsedByChronos single-alignment gene with the library's modal guide count is fitted; evaluation uses Achilles reference genes with gene-held-out controls."

FIGSHARE_ARTICLE = 27993248
LIBRARIES = {"Avana": 4, "KY": 5}
SOURCE_FILES = ("AchillesCommonEssentialControls.csv", "AchillesNonessentialControls.csv",
                "AvanaGuideMap.csv", "AvanaRawReadcounts.csv", "KYGuideMap.csv", "KYRawReadcounts.csv",
                "ScreenSequenceMap.csv")


def digest(path):
    sha = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 24), b''):
            sha.update(block)
    return sha.hexdigest()


def eligible_cells(seq):
    usable = seq[(seq.ModelID != 'pDNA') & seq.PassesQC.astype(bool) & (seq.ScreenType == '2DS')
                 & ~seq.DrugTreated.fillna(False).astype(bool)]
    per_library = {}
    for library in LIBRARIES:
        rows = usable[usable.Library == library]
        replicates = rows.groupby('ModelID').size()
        batches = rows.groupby('ModelID').pDNABatch.nunique()
        screens = rows.groupby('ModelID').ScreenID.nunique()
        per_library[library] = set(replicates[(replicates >= 2) & (batches == 1) & (screens == 1)].index)
    return sorted(per_library['Avana'] & per_library['KY']), usable


def prepare(source, output, n_cells):
    import pandas as pd
    from hit_quality import heldout_folds
    provenance = {'schema_version': 1, 'release': 'DepMap 24Q4 Public', 'figshare_article': FIGSHARE_ARTICLE,
                  'chronos_commit': frontier.CHRONOS_COMMIT, 'jacks_commit': frontier.JACKS_COMMIT,
                  'seed': 1729, 'folds': 5, 'chronos_epochs': 1001,
                  'selection': f'first {n_cells} sorted ModelIDs in Avana and KY with >=2 passing 2D untreated replicates, one screen and one pDNA batch',
                  'scope': 'genome-wide: every UsedByChronos single-alignment gene with the modal guide count',
                  'sources': {name: digest(source / name) for name in SOURCE_FILES}, 'libraries': {}}
    classes = {}
    for name, value in [('AchillesCommonEssentialControls.csv', 1), ('AchillesNonessentialControls.csv', 0)]:
        for gene in pd.read_csv(source / name).Gene:
            if gene in classes and classes[gene] != value:
                raise ValueError('Conflicting reference annotation')
            classes[gene] = value
    seq_all = pd.read_csv(source / 'ScreenSequenceMap.csv')
    shared, usable = eligible_cells(seq_all)
    cells = tuple(shared[:n_cells])
    if len(cells) < n_cells:
        raise ValueError(f'Only {len(cells)} eligible cells')
    provenance['cell_lines'] = cells
    provenance['eligible_shared_cells'] = len(shared)
    for library, guide_count in LIBRARIES.items():
        folder = output / library
        folder.mkdir(parents=True, exist_ok=True)
        treatments = usable[(usable.Library == library) & usable.ModelID.isin(cells)]
        baselines = seq_all[(seq_all.ModelID == 'pDNA') & (seq_all.Library == library)
                            & seq_all.pDNABatch.isin(treatments.pDNABatch)]
        selected = pd.concat([baselines, treatments]).reset_index(drop=True)
        seq = pd.DataFrame({'sequence_ID': selected.SequenceID,
                            'cell_line_name': selected.ModelID,
                            'pDNA_batch': selected.pDNABatch,
                            'days': selected.Days})
        guide_map = pd.read_csv(source / f'{library}GuideMap.csv')
        guide_map = guide_map[guide_map.UsedByChronos.astype(bool) & (guide_map.nAlignments == 1) & guide_map.Gene.notna()]
        mapping = guide_map[['sgRNA', 'Gene']].rename(columns={'sgRNA': 'sgrna', 'Gene': 'gene'})
        reads = pd.read_csv(source / f'{library}RawReadcounts.csv', index_col=0,
                            usecols=lambda c: c == 'Unnamed: 0' or c in set(seq.sequence_ID))
        missing = set(seq.sequence_ID) - set(reads.columns)
        if missing:
            raise ValueError(f'Missing read columns: {sorted(missing)[:5]}')
        counts = reads.T.loc[seq.sequence_ID]
        mapping = mapping[mapping.sgrna.isin(counts.columns)]
        counts = counts[mapping.sgrna]
        finite = counts.columns[np.isfinite(counts.values).all(axis=0)]
        mapping = mapping[mapping.sgrna.isin(finite)]
        sizes = mapping.groupby('gene').size()
        mapping = mapping[mapping.gene.isin(sizes[sizes == guide_count].index)]
        counts = counts[mapping.sgrna]
        if not (counts.values >= 0).all():
            raise ValueError('Negative raw count')
        genes = set(mapping.gene)
        folds = heldout_folds([g for g in genes if classes.get(g) == 1],
                              [g for g in genes if classes.get(g) == 0])
        mapping.to_csv(folder / 'guide_map.csv', index=False)
        seq.to_csv(folder / 'sequence_map.csv', index=False)
        counts.to_csv(folder / 'readcounts.csv')
        table = mapping.rename(columns={'sgrna': 'sgRNA', 'gene': 'Gene'}).set_index('sgRNA').join(counts.T)
        table.to_csv(folder / 'counts.tsv', sep='\t')
        pd.DataFrame([{'gene': g, 'essential': classes[g], 'fold': folds[g]} for g in sorted(folds)]).to_csv(
            folder / 'labels.csv', index=False)
        for fold in range(5):
            controls = sorted(g for g in folds if folds[g] != fold and classes[g] == 0)
            (folder / f'controls-{fold}.txt').write_text('\n'.join(controls) + '\n')
        for cell in cells:
            cell_map = seq[seq.cell_line_name == cell]
            matched = seq[(seq.cell_line_name == 'pDNA') & seq.pDNA_batch.isin(set(cell_map.pDNA_batch))]
            if cell_map.pDNA_batch.nunique() != 1 or len(matched) != 1:
                raise ValueError(f'{library} {cell}: requires one matched pDNA sample')
            samples = matched.sequence_ID.tolist() + cell_map.sequence_ID.tolist()
            whitespace_free = {s: '_'.join(s.split()) for s in samples}
            table[['Gene'] + samples].rename(columns=whitespace_free).to_csv(folder / f'counts-{cell}.tsv', sep='\t')
            (folder / f'design-{cell}.tsv').write_text('Samples\tbaseline\ttreatment\n' + ''.join(
                f'{whitespace_free[s]}\t1\t{int(i > 0)}\n' for i, s in enumerate(samples)))
        provenance['libraries'][library] = {'genes': len(genes), 'guides': len(mapping), 'guide_count': guide_count,
                                            'samples': len(seq), 'evaluated_genes': len(folds),
                                            'counts_sha256': digest(folder / 'counts.tsv')}
        print(library, provenance['libraries'][library], flush=True)
    (output / 'panel.json').write_text(json.dumps(provenance, indent=2) + '\n')


def panel_cells(panel):
    return tuple(json.loads((panel / 'panel.json').read_text())['cell_lines'])


def run_job(command, output, timeout):
    output.mkdir(parents=True, exist_ok=True)
    if (output / 'command.json').exists() and json.loads((output / 'command.json').read_text())['returncode'] == 0:
        return 'cached'
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', TF_CPP_MIN_LOG_LEVEL='2')
    start = time.perf_counter()
    with (output / 'stdout.log').open('w') as stream:
        completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, timeout=timeout, env=env)
    (output / 'command.json').write_text(json.dumps({'command': command, 'returncode': completed.returncode,
                                                     'wall_seconds': time.perf_counter() - start}, indent=2) + '\n')
    return 'ok' if completed.returncode == 0 else f'failed {output}'


def jobs_for(method, panel, python, output_tag):
    cells = panel_cells(panel)
    for library in LIBRARIES:
        folder = panel / library
        for fold in range(5):
            if method == 'fit':
                for cell in cells:
                    output = folder / 'fit' / f'{fold}-{cell}'
                    yield [str(python), '-m', 'crisprworks_fit', 'mle',
                           '-k', str(folder / f'counts-{cell}.tsv'), '-d', str(folder / f'design-{cell}.tsv'),
                           '-n', str(output / 'screen'), '--kernel', 'native', '--norm-method', 'control',
                           '--control-gene', str(folder / f'controls-{fold}.txt'), '--permutation-pvalues', 'finite',
                           '--seed', '1729', '--permutation-round', '2', '--genes-varmodeling', '1000',
                           '--update-efficiency', '--write-fit-details'], output
            else:
                output = folder / (output_tag or method) / str(fold)
                yield [str(python), str(Path(__file__).resolve()), '--method', method, '--panel', str(panel),
                       '--library', library, '--fold', str(fold), '--out', str(output)], output


def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--method', choices=('prepare', 'chronos', 'jacks', 'fit', 'joint'), required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--panel', type=Path, required=True)
    parser.add_argument('--cells', type=int, default=24)
    parser.add_argument('--library', choices=tuple(LIBRARIES))
    parser.add_argument('--fold', type=int, choices=range(5))
    parser.add_argument('--out', type=Path)
    parser.add_argument('--python', type=Path)
    parser.add_argument('--output-tag', choices=('joint', 'joint_numpy126'))
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--timeout', type=float, default=6 * 3600)
    args = parser.parse_args()
    panel = args.panel.resolve()
    if args.method == 'prepare':
        if not args.source:
            parser.error('prepare requires --source')
        panel.mkdir(parents=True, exist_ok=True)
        prepare(args.source.resolve(), panel, args.cells)
    elif args.python:
        jobs = list(jobs_for(args.method, panel, args.python, args.output_tag))
        with ThreadPoolExecutor(args.workers) as pool:
            for (command, output), status in zip(jobs, pool.map(lambda j: run_job(j[0], j[1], args.timeout), jobs)):
                print(args.method, output.relative_to(panel), status, flush=True)
    else:
        if args.library is None or args.fold is None or args.out is None:
            parser.error('worker requires library, fold, out')
        frontier.CELL_LINES = panel_cells(panel)
        args.out.mkdir(parents=True, exist_ok=True)
        worker = {'chronos': frontier.chronos_worker, 'jacks': frontier.jacks_worker, 'joint': frontier.joint_worker}
        worker[args.method](panel / args.library, args.fold, args.out)


if __name__ == '__main__':
    main()
