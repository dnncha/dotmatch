#!/usr/bin/env python3
"""Verify archived study outputs and regenerate descriptive research summaries.

No raw sequencing input, new counts, p-values or biological hit calls are made.
Directories are the output directories extracted from the recorded Actions ZIPs.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def tsv(path):
    with Path(path).open(newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def verify(root):
    root = Path(root).resolve()
    manifest = read_json(root / 'MANIFEST.json')
    require(manifest['status'] == 'complete', f'incomplete source: {root}')
    require(manifest['files'], 'empty manifest')
    for rel, expected in manifest['files'].items():
        path = (root / rel).resolve()
        require(path.is_relative_to(root), 'unsafe manifest path')
        require(path.is_file(), f'missing artifact: {rel}')
        require(path.stat().st_size == expected['bytes'] and sha(path) == expected['sha256'], f'corrupt artifact: {rel}')
    return {'files_verified': len(manifest['files']), 'manifest_sha256': sha(root / 'MANIFEST.json')}


def annotate_trace(path, genes, edges):
    result = {'multi_offset_records': 0, 'same_gene': 0, 'cross_gene': 0,
              'records_with_exact_overlap_pair': 0, 'records_with_all_pairs_exact_overlap': 0,
              'observed_pairs': 0, 'pairs_with_exact_overlap': 0}
    for row in tsv(path):
        calls = sorted((int(event.split(':', 1)[0]), event.split(':', 1)[1]) for event in row['offset_and_guide'].split(';'))
        require(len(calls) > 1 and len({x[0] for x in calls}) == len(calls), 'invalid multi-offset trace')
        annotations = {genes[target] for _, target in calls}
        pairs = [(b[0] - a[0], a[1], b[1]) for a, b in itertools.combinations(calls, 2)]
        supported = sum(pair in edges for pair in pairs)
        result['multi_offset_records'] += 1
        result['same_gene'] += len(annotations) == 1
        result['cross_gene'] += len(annotations) > 1
        result['records_with_exact_overlap_pair'] += supported > 0
        result['records_with_all_pairs_exact_overlap'] += supported == len(pairs)
        result['observed_pairs'] += len(pairs)
        result['pairs_with_exact_overlap'] += supported
    return result


def summarize(args):
    roots = {'pilot': Path(args.pilot), 'yusa_prefix': Path(args.yusa_prefix),
             'brunello_prefix': Path(args.brunello_prefix), 'overlap': Path(args.overlap)}
    if args.full_yusa:
        roots['full_yusa'] = Path(args.full_yusa)
    verification = {label: verify(path) for label, path in roots.items()}
    pilot = read_json(roots['pilot'] / 'run.json')
    yusa = read_json(roots['yusa_prefix'] / 'reconstruction.json')
    brunello = read_json(roots['brunello_prefix'] / 'run.json')
    overlap = read_json(roots['overlap'] / 'run.json')
    result = {'schema': 'dotmatch.research.AR001.summary.v1', 'status': 'complete',
        'analysis_kind': 'descriptive and algorithmic validation; not biological inference',
        'summary_script_sha256': sha(__file__), 'verification': verification,
        'pilot': {'evaluation_reads': sum(s['records'] for s in pilot['samples']), 'samples': pilot['samples'],
                  'libraries': pilot['libraries']},
        'prefix': {}, 'overlap': overlap['libraries'], 'full_yusa': None,
        'read_count_warning': 'Pilot and historical prefixes overlap one another and, for Yusa, overlap full archives. Do not sum analysis counts as independent reads.'}
    require(result['pilot']['evaluation_reads'] == 600000 and len(pilot['samples']) == 6, 'pilot cohort mismatch')
    for sample in pilot['samples']:
        require(sample['independent_count_state_transition_changed_call_agreement'], 'pilot independent validation failed')
        for policy, outcomes in sample['states'].items():
            require(sum(outcomes.values()) == sample['records'], 'pilot read conservation failed')
    for name, data in [('yusa', yusa), ('brunello', brunello)]:
        require(data['status'] == 'complete' and data['historical_report_exactly_reproduced'], f'historical reproduction failed: {name}')
        require(pilot['libraries'][name]['provenance']['normalized_sha256'] == data['library']['normalized_sha256'] == overlap['libraries'][name]['source']['normalized_sha256'], f'reference differs across lanes: {name}')
        edge_rows = tsv(roots['overlap'] / f'{name}.exact_overlap_edges.tsv')
        edges = {(int(r['shift']), r['guide_a_at_earlier_offset'], r['guide_b_at_later_offset']) for r in edge_rows}
        require(len(edges) == len(edge_rows), 'duplicate exact overlap edge')
        count_path = roots[f'{name}_prefix'] / ('guide_counter_hamming.counts.txt' if name == 'yusa' else 'guide_counter.counts.txt')
        count_rows = tsv(count_path)
        genes = {r['guide']: r['gene'] for r in count_rows}
        require(len(genes) == len(count_rows), 'duplicate guide IDs')
        for row in edge_rows:
            require(genes[row['guide_a_at_earlier_offset']] == row['gene_a'] and genes[row['guide_b_at_later_offset']] == row['gene_b'], 'overlap gene annotation mismatch')
            require(int(row['same_gene']) == (row['gene_a'] == row['gene_b']), 'overlap same-gene flag mismatch')
        annotated, totals = {}, {k: 0 for k in ('records', 'matching_events', 'matched_read_records', 'extra_events', 'multi_offset_records', 'same_gene', 'cross_gene', 'records_with_exact_overlap_pair')}
        for sample, values in data['results'].items():
            require(values['reconstruction_count_disagreements'] == 0, 'event reconstruction failed')
            trace = annotate_trace(roots[f'{name}_prefix'] / f'{sample}.multi_offset_records.tsv', genes, edges)
            require(trace['multi_offset_records'] == values['multi_offset_same_guide_reads'] + values['multi_offset_different_guide_reads'], 'multi-record denominator mismatch')
            require(sum(int(r[sample]) for r in count_rows) == values['matching_events'], 'event matrix sum mismatch')
            require(sum(int(k) * v for k, v in values['event_multiplicity'].items()) == values['matching_events'], 'event multiplicity mismatch')
            require(sum(values['event_multiplicity'].values()) == values['records'], 'read multiplicity mismatch')
            require(values['matching_events'] - values['reads_with_at_least_one_event'] == values['excess_events_above_matched_reads'], 'extra event mismatch')
            annotated[sample] = {**values, **trace}
            totals['records'] += values['records']
            totals['matching_events'] += values['matching_events']
            totals['matched_read_records'] += values['reads_with_at_least_one_event']
            totals['extra_events'] += values['excess_events_above_matched_reads']
            for key in ('multi_offset_records', 'same_gene', 'cross_gene', 'records_with_exact_overlap_pair'):
                totals[key] += trace[key]
        result['prefix'][name] = {'historical_aggregate': data['historical_aggregate'], 'samples': annotated, 'totals': totals}
    if args.full_yusa:
        full = read_json(roots['full_yusa'] / 'run.json')
        require(full['status'] == 'complete' and len(full['results']) == 2, 'full archive experiment incomplete')
        require(full['library']['normalized_sha256'] == yusa['library']['normalized_sha256'], 'full/prefix reference mismatch')
        for label, sample in full['results'].items():
            source = full['samples'][label]
            require(source['archive_md5_locally_verified'] and source['local_md5'] == source['fastq_md5'], 'full MD5 verification absent')
            require(source['local_bytes'] == int(source['fastq_bytes']) and sample['records'] == int(source['read_count']), 'full metadata denominator mismatch')
            require(sample['independent_all_guide_count_agreement'], 'full independent count agreement absent')
            for outcomes in sample['states'].values():
                require(sum(outcomes.values()) == sample['records'], 'full state conservation failed')
            require(sample['multi_offset_same_gene_reads'] + sample['multi_offset_cross_gene_reads'] == sample['multi_offset_same_guide_reads'] + sample['multi_offset_different_guide_reads'], 'gene classification denominator mismatch')
        result['full_yusa'] = {'git_commit': full['git_commit'], 'samples': full['samples'], 'results': full['results'],
                               'records': sum(s['records'] for s in full['results'].values())}
        result['full_yusa']['totals'] = {key: sum(s[key] for s in full['results'].values()) for key in (
            'records', 'changed_records', 'matching_events', 'reads_with_any_event',
            'excess_events_above_matched_reads', 'multi_offset_same_gene_reads', 'multi_offset_cross_gene_reads',
            'multi_offset_same_guide_reads', 'multi_offset_different_guide_reads')}
        result['full_yusa']['totals']['fixed_policy_unique'] = {
            policy: sum(s['states'][policy]['unique'] for s in full['results'].values())
            for policy in ('exact', 'radius_k1', 'best_k1')}
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'results.json').write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
    prefix_header = ['library', 'sample', 'read_records', 'matching_events', 'matched_read_records', 'multi_offset_records', 'same_gene', 'cross_gene', 'exact_overlap_supported_records']
    with (out / 'prefix_accounting.tsv').open('w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
        writer.writerow(prefix_header)
        for library, data in result['prefix'].items():
            for sample, v in data['samples'].items():
                writer.writerow([library, sample, v['records'], v['matching_events'], v['reads_with_at_least_one_event'], v['multi_offset_records'], v['same_gene'], v['cross_gene'], v['records_with_exact_overlap_pair']])
    with (out / 'fixed_window_pilot.tsv').open('w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
        writer.writerow(['sample', 'offset', 'evaluated_records', 'exact_unique', 'radius_unique', 'best_unique', 'changed_records'])
        for s in result['pilot']['samples']:
            writer.writerow([s['sample'], s['target_start'], s['records'], *[s['states'][p]['unique'] for p in ('exact', 'radius_k1', 'best_k1')], s['changed_records']])
    print(json.dumps({'verified_files': sum(v['files_verified'] for v in verification.values()), 'pilot_reads': result['pilot']['evaluation_reads'], 'prefix_totals': {k:v['totals'] for k,v in result['prefix'].items()}, 'full_yusa_records': result['full_yusa']['records'] if result['full_yusa'] else None}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ('pilot', 'yusa-prefix', 'brunello-prefix', 'overlap', 'out'):
        parser.add_argument('--' + field, required=True)
    parser.add_argument('--full-yusa')
    summarize(parser.parse_args())
