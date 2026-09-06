#!/usr/bin/env python3
"""Independent, read-conserving audit of pinned public CRISPR counting results.

This research harness is not a production replacement counter or truth oracle.
It separates fixed-window policies, native multi-offset count events, and
constructed known-origin controls. Outputs are complete only after reconciliation.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
from functools import lru_cache
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import random
import shutil
import subprocess
import sys
import tempfile
import traceback
import unittest

POLICIES = ('exact', 'radius_k1', 'best_k1')
STATE = {-1: 'none', -2: 'ambiguous', -3: 'invalid'}
UPSTREAM_SHA = '96602fd0b9732204b530afb912ff679d48b0ba9e13d32c5eea67c10cbbdbf777'
ENGINE_COMMIT = '11d159fa1648365f2a4e96917b483c33aa5d9fe7'


def sha(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def dump(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.pending')
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + '\n')
    tmp.replace(path)


def tsv(path, header, rows):
    with Path(path).open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)


def one(root, pattern):
    matches = list(Path(root).rglob(pattern))
    if len(matches) != 1:
        raise ValueError(f'Expected exactly one {pattern}; found {len(matches)}')
    return matches[0]


def library(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    if not rows or any(not r.get(k) for r in rows for k in ('id', 'sequence', 'gene')):
        raise ValueError('Empty library or missing target fields')
    if len({r['id'] for r in rows}) != len(rows):
        raise ValueError('Duplicate target IDs')
    length = len(rows[0]['sequence'])
    if length < 2 or any(len(r['sequence']) != length or not set(r['sequence']) <= set('ACGT') for r in rows):
        raise ValueError('Expected uniform ACGT reference targets')
    return rows


def fastq(path):
    with gzip.open(path, 'rt', encoding='ascii', newline='') as f:
        ordinal = 0
        while True:
            head = f.readline()
            if not head:
                return
            seq, plus, qual = f.readline(), f.readline(), f.readline()
            ordinal += 1
            if not head.startswith('@') or not head[1:].split() or not plus.startswith('+') or not qual:
                raise ValueError(f'Malformed FASTQ record {ordinal}')
            seq, qual = seq.rstrip('\r\n'), qual.rstrip('\r\n')
            if not seq or len(seq) != len(qual) or any(c.isspace() for c in seq):
                raise ValueError(f'Invalid sequence/quality record {ordinal}')
            if min(qual) < '!' or max(qual) > '~':
                raise ValueError(f'Invalid quality at record {ordinal}')
            yield ordinal, head[1:].split()[0], seq.upper()


class Index:
    """At distance <=1 at least one of two disjoint halves is unchanged."""
    def __init__(self, sequences):
        self.sequences = tuple(sequences)
        self.length = len(self.sequences[0])
        self.cut = self.length // 2
        self.left, self.right, self.exact = defaultdict(list), defaultdict(list), defaultdict(list)
        for i, s in enumerate(self.sequences):
            if len(s) != self.length:
                raise ValueError('Mixed target lengths')
            self.left[s[:self.cut]].append(i)
            self.right[s[self.cut:]].append(i)
            self.exact[s].append(i)
        self.candidates = lru_cache(maxsize=131072)(self._candidates)

    def _candidates(self, window):
        if window is None or len(window) != self.length:
            return ()
        possible = set(self.left.get(window[:self.cut], ()))
        possible.update(self.right.get(window[self.cut:], ()))
        out = []
        for i in sorted(possible):
            distance = sum(a != b for a, b in zip(window, self.sequences[i]))
            if distance <= 1:
                out.append((i, distance))
        return tuple(out)

    def calls(self, window):
        if window is None:
            return (-3, -3, -3)
        candidates = self.candidates(window)
        exact = tuple(i for i, d in candidates if d == 0)
        e = exact[0] if len(exact) == 1 else -2 if exact else -1
        r = candidates[0][0] if len(candidates) == 1 else -2 if candidates else -1
        if candidates:
            distance = min(d for _, d in candidates)
            best = [i for i, d in candidates if d == distance]
            b = best[0] if len(best) == 1 else -2
        else:
            b = -1
        return e, r, b

    def upstream(self, window, exact_only=False):
        if window is None or len(window) != self.length or not set(window) <= set('ACGT'):
            return None
        exact = self.exact.get(window)
        if exact:
            return exact[-1]  # This explicitly reproduces upstream duplicate overwrite semantics.
        if exact_only:
            return None
        candidates = self.candidates(window)
        return candidates[0][0] if len(candidates) == 1 else None


def oracle(sequences, window):
    if window is None:
        return (-3, -3, -3)
    distances = [sum(a != b for a, b in zip(window, s)) for s in sequences]
    exact = [i for i, d in enumerate(distances) if d == 0]
    radius = [i for i, d in enumerate(distances) if d <= 1]
    best = [i for i in radius if distances[i] == min(distances)]
    choose = lambda xs: xs[0] if len(xs) == 1 else -2 if xs else -1
    return choose(exact), choose(radius), choose(best)


def mode(counter):
    if not counter:
        return 'N'
    top = max(counter.values())
    winners = [base for base, count in counter.items() if count == top]
    return winners[0] if len(winners) == 1 else 'N'


def discover(path, index, start, limit=100000, exact_only=False):
    offsets, left, right = Counter(), [], []
    exact_reads, records = 0, 0
    for ordinal, _, seq in fastq(path):
        if ordinal > limit:
            break
        records += 1
        for offset in range(max(0, len(seq) - index.length + 1)):
            if index.upstream(seq[offset:offset + index.length], exact_only) is not None:
                offsets[offset] += 1
        w = seq[start:start + index.length] if len(seq) >= start + index.length else None
        if w is not None and len(index.exact.get(w, ())) == 1:
            exact_reads += 1
            for counts, text in ((left, seq[:start]), (right, seq[start + index.length:])):
                while len(counts) < len(text):
                    counts.append(Counter())
                for pos, base in enumerate(text):
                    counts[pos][base] += 1
    events = sum(offsets.values())
    selected = sorted(k for k, n in offsets.items() if events and n / events >= 0.0025)
    return {'records': records, 'event_denominator': events, 'offset_events': dict(sorted(offsets.items())),
            'selected': selected, 'fixed_exact_discovery_reads': exact_reads,
            'left': ''.join(mode(c) for c in left), 'right': ''.join(mode(c) for c in right),
            'flank_support': {'left': left, 'right': right}}


def counted_events(index, seq, offsets, exact_only=False):
    out = []
    for offset in offsets:
        if offset + index.length <= len(seq):
            target = index.upstream(seq[offset:offset + index.length], exact_only)
            if target is not None:
                out.append((offset, target))
    return out


def add_events(stats, events, rows):
    stats['reads'] += 1
    stats['count_events'] += len(events)
    if not events:
        return
    stats['matched_reads'] += 1
    distinct = {i for _, i in events}
    genes = {rows[i]['gene'] for i in distinct}
    stats['extra_events'] += len(events) - 1
    stats['repeated_same_target_events'] += len(events) - len(distinct)
    stats['additional_distinct_target_events'] += len(distinct) - 1
    if len(events) > 1:
        stats['multiply_counted_reads'] += 1
    if len(distinct) > 1:
        stats['multiple_distinct_guides'] += 1
        stats['within_gene_multiple' if len(genes) == 1 else 'cross_gene_multiple'] += 1


def check_conservation(stats):
    if stats['count_events'] != stats['matched_reads'] + stats['extra_events']:
        raise AssertionError('Count-event conservation failed')
    if stats['extra_events'] != stats['repeated_same_target_events'] + stats['additional_distinct_target_events']:
        raise AssertionError('Multiplicity decomposition failed')


def read_counts(path, rows, sample):
    with Path(path).open(newline='') as f:
        table = csv.DictReader(f, delimiter='\t')
        fields = table.fieldnames
        if not fields or sample not in fields or len(fields) < 3:
            raise ValueError('Invalid count table columns')
        ident, gene = fields[:2]
        lookup = {}
        for r in table:
            if r[ident] in lookup or not r[sample].isdigit():
                raise ValueError('Duplicate target row or non-integer count')
            lookup[r[ident]] = (r[gene], int(r[sample]))
    if set(lookup) != {r['id'] for r in rows}:
        raise ValueError('Count/reference target identities differ')
    if any(lookup[r['id']][0] != r['gene'] for r in rows):
        raise ValueError('Count/reference gene identities differ')
    return [lookup[r['id']][1] for r in rows]


def run_upstream(binary, path, lib, sample, output, exact_only=False):
    command = [str(binary), 'count', '--input', str(path), '--samples', sample, '--library', str(lib),
               '--output', str(output), '--offset-sample-size', '100000', '--offset-min-fraction', '0.0025']
    if exact_only:
        command.append('--exact-match')
    with Path(str(output) + '.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    return command


def synthetic_controls(out, rows, index, context, binary, lib, accession):
    from dotmatch.sensitivity import run_sensitivity
    experiments = [('balanced', [(i, 1) for i in range(len(rows))], context['left'])]
    witness_id = 'chr11:69998591-69998614'
    witness = next((i for i, r in enumerate(rows) if r['id'] == witness_id), None)
    if accession == 'ERR376998':
        if witness is None or not context['left']:
            raise ValueError('Prespecified witness is absent')
        changed = next(b for b in 'ACGT' if b != context['left'][-1])
        experiments += [('witness_control', [(witness, 128)], context['left']),
                        ('witness_depleted', [(witness, 16)], context['left']),
                        ('witness_flank_changed', [(witness, 128)], context['left'][:-1] + changed)]
    results = []
    for name, origins, left in experiments:
        folder = out / 'synthetic' / name
        folder.mkdir(parents=True, exist_ok=False)
        path = folder / 'synthetic.fastq.gz'
        truth = [0] * len(rows)
        with path.open('xb') as raw, gzip.GzipFile(fileobj=raw, filename='', mode='wb', mtime=0) as f:
            ordinal = 0
            for i, n in origins:
                truth[i] += n
                seq = left + rows[i]['sequence'] + context['right']
                for _ in range(n):
                    ordinal += 1
                    f.write(f'@synthetic_{ordinal}\n{seq}\n+\n{"I" * len(seq)}\n'.encode('ascii'))
        tsv(folder / 'truth.tsv', ['id', 'gene', 'known_source_reads'],
            ((r['id'], r['gene'], truth[i]) for i, r in enumerate(rows)))
        native = run_sensitivity(targets=lib, reads=path, target_start=len(left), target_length=index.length,
                                 sample_label=name, out_dir=folder / 'native')
        if read_counts(folder / 'native/exact.counts.tsv', rows, name) != truth:
            raise AssertionError('Fixed exact native counts differ from known synthetic origins')
        for exact_only in (False, True):
            policy = 'exact' if exact_only else 'one_mismatch'
            prefix = folder / ('upstream_' + policy)
            command = run_upstream(binary, path, lib, name, prefix, exact_only)
            offsets = discover(path, index, len(left), exact_only=exact_only)['selected']
            independent, stats = [0] * len(rows), Counter()
            for _, _, seq in fastq(path):
                events = counted_events(index, seq, offsets, exact_only)
                add_events(stats, events, rows)
                for _, i in events:
                    independent[i] += 1
            if independent != read_counts(str(prefix) + '.counts.txt', rows, name):
                raise AssertionError('Synthetic upstream/reference mismatch')
            check_conservation(stats)
            results.append({'experiment': name, 'data_kind': 'synthetic_known_origin_no_errors', 'policy': policy,
                            'source_reads': sum(truth), 'count_events': sum(independent), 'statistics': stats,
                            'changed_guide_rows': sum(a != b for a, b in zip(truth, independent)),
                            'zero_origin_guides_receiving_counts': sum(t == 0 and c > 0 for t, c in zip(truth, independent)),
                            'counts_to_zero_origin_guides': sum(c for t, c in zip(truth, independent) if t == 0),
                            'all_guide_reconciliation': True, 'fixed_exact_native_matches_truth': True,
                            'selected_offsets': offsets, 'command': command, 'input_sha256': sha(path)})
    return results


def audit(args):
    import numpy as np
    output = Path(args.out)
    if output.exists():
        raise FileExistsError(output)
    stage = output.with_name(output.name + '.pending')
    stage.mkdir(parents=True, exist_ok=False)
    lib = one(args.raw, 'library.tsv')
    path = one(args.raw, args.accession + '.fastq.gz')
    provenance = json.loads(one(args.evidence, args.accession + '.provenance.json').read_text())
    library_provenance = json.loads(one(args.evidence, 'library-provenance.json').read_text())
    native_paths = [p for p in Path(args.evidence).rglob('summary.json')
                    if p.parent.name == 'full' and p.parent.parent.name == args.accession]
    if len(native_paths) != 1:
        raise ValueError('Expected exactly one native full-run summary')
    native_path = native_paths[0]
    native = json.loads(native_path.read_text())
    if native['completion'] != 'complete' or native['parameters']['target_start'] != args.start:
        raise ValueError('Native completion or extraction mismatch')
    n = int(provenance['observed_records'])
    if n != int(provenance['read_count']) or n != native['read_count']:
        raise ValueError('ENA/native record totals differ')
    if path.stat().st_size != int(provenance['fastq_bytes']):
        raise ValueError('Archive byte count mismatch')
    if sha(path, 'md5') != provenance['fastq_md5'] or sha(path) != provenance['sha256']:
        raise ValueError('Archive digest mismatch')
    if sha(lib) != library_provenance['derived_sha256']:
        raise ValueError('Derived library digest mismatch')
    binary = Path(args.binary).resolve()
    if sha(binary) != UPSTREAM_SHA:
        raise ValueError('Unpinned upstream executable')
    rows = library(lib)
    index = Index([r['sequence'] for r in rows])
    if native['parameters']['target_length'] != index.length:
        raise ValueError('Native target length mismatch')
    code_identity = {'replay_source_sha256': sha(__file__), 'source_commit': os.getenv('GITHUB_SHA'),
                     'workflow_run': os.getenv('GITHUB_RUN_ID'), 'engine_commit': ENGINE_COMMIT,
                     'upstream_binary_sha256': sha(binary), 'python': sys.version, 'platform': platform.platform()}
    dump(stage / 'execution-start.json', code_identity)
    context = discover(path, index, args.start)
    if not context['selected'] or len(context['left']) != args.start:
        raise ValueError('Cannot construct discovery-only context model')
    model_rows, predictions = [], []
    for i, r in enumerate(rows):
        events = counted_events(index, context['left'] + r['sequence'] + context['right'], context['selected'])
        predictions.append(len(events) > 1)
        model_rows.append((r['id'], r['gene'], len(events), json.dumps([(o, rows[j]['id']) for o, j in events])))
    tsv(stage / 'template-model.tsv', ['source_guide', 'gene', 'predicted_count_events', 'events'], model_rows)
    dump(stage / 'discovery-model.json', context)
    model_hash = sha(stage / 'template-model.tsv')  # Written and closed before evaluation.
    command = run_upstream(binary, path, lib, args.accession, stage / 'upstream')
    counts = {mode: [0] * len(rows) for mode in POLICIES}
    states = {mode: Counter({'unique': 0, 'none': 0, 'ambiguous': 0, 'invalid': 0}) for mode in POLICIES}
    up_counts, stats, heldout, confusion = [0] * len(rows), Counter(), Counter(), Counter()
    event_sets, examples = Counter(), []
    selected_ordinals = set(random.Random(20260906).sample(range(1, n + 1), min(200, n)))
    reference = np.frombuffer(''.join(index.sequences).encode('ascii'), dtype=np.uint8).reshape(len(rows), index.length)
    oracle_checked, evaluated, excluded = 0, 0, 0
    observed_records = changed_reads = 0
    for ordinal, read_id, seq in fastq(path):
        observed_records = ordinal
        window = seq[args.start:args.start + index.length] if len(seq) >= args.start + index.length else None
        calls = index.calls(window)
        changed_reads += len(set(calls)) > 1
        for policy, call in zip(POLICIES, calls):
            states[policy]['unique' if call >= 0 else STATE[call]] += 1
            if call >= 0:
                counts[policy][call] += 1
        if ordinal in selected_ordinals:
            if window is None:
                expected = (-3, -3, -3)
            else:
                distances = np.count_nonzero(reference != np.frombuffer(window.encode('ascii'), dtype=np.uint8), axis=1)
                xs = [np.flatnonzero(distances == 0), np.flatnonzero(distances <= 1)]
                xs.append(np.flatnonzero(distances == distances.min()) if distances.min() <= 1 else np.array([], dtype=int))
                expected = tuple(int(x[0]) if len(x) == 1 else -2 if len(x) else -1 for x in xs)
            if expected != calls:
                raise AssertionError(f'All-target oracle disagreement at {ordinal}')
            oracle_checked += 1
        events = counted_events(index, seq, context['selected'])
        add_events(stats, events, rows)
        for _, i in events:
            up_counts[i] += 1
        if len(events) > 1:
            key = tuple(sorted(set(i for _, i in events)))
            event_sets[key] += 1
            if len(examples) < 20:
                examples.append({'record_ordinal': ordinal, 'read_id': read_id,
                                 'events': [{'offset': o, 'target_id': rows[i]['id'], 'gene': rows[i]['gene']} for o, i in events]})
        if ordinal > context['records']:
            add_events(heldout, events, rows)
            if calls[0] >= 0:
                evaluated += 1
                predicted, observed = predictions[calls[0]], len(events) > 1
                confusion[('TP' if observed else 'FP') if predicted else ('FN' if observed else 'TN')] += 1
            else:
                excluded += 1
    if observed_records != n or sha(path) != provenance['sha256']:
        raise AssertionError('Full FASTQ changed or record count differs')
    check_conservation(stats)
    check_conservation(heldout)
    for policy in POLICIES:
        if counts[policy] != read_counts(native_path.parent / (policy + '.counts.tsv'), rows, args.accession):
            raise AssertionError('Native count disagreement: ' + policy)
        if dict(states[policy]) != native['outcomes'][policy]:
            raise AssertionError('Native state disagreement: ' + policy)
    if changed_reads != native['changed_reads']:
        raise AssertionError('Native changed-read count differs')
    if up_counts != read_counts(stage / 'upstream.counts.txt', rows, args.accession):
        raise AssertionError('Upstream per-guide count disagreement')
    if sha(stage / 'template-model.tsv') != model_hash:
        raise AssertionError('Model changed during evaluation')
    tsv(stage / 'all-guide-counts.tsv', ['id', 'gene', *POLICIES, 'upstream_count_events'],
        ((r['id'], r['gene'], *(counts[p][i] for p in POLICIES), up_counts[i]) for i, r in enumerate(rows)))
    tsv(stage / 'multiple-guide-read-classes.tsv', ['target_ids', 'genes', 'read_count'],
        ((json.dumps([rows[i]['id'] for i in key]), json.dumps(sorted({rows[i]['gene'] for i in key})), count)
         for key, count in sorted(event_sets.items())))
    dump(stage / 'examples.json', {'selection': 'first_20_multiple_count_records_not_a_representative_sample', 'examples': examples})
    controls = synthetic_controls(stage, rows, index, context, binary, lib, args.accession)
    result = {'completion': 'complete', 'accession': args.accession, 'scope': 'complete_public_archive',
              'target_count': len(rows), 'start': args.start, 'length': index.length, 'input': provenance,
              'library_sha256': sha(lib), 'source': code_identity, 'upstream_command': command,
              'statistics': stats, 'heldout_statistics': heldout, 'native_outcomes': states,
              'native_changed_reads': changed_reads, 'selected_offsets': context['selected'],
              'all_guide_native_comparisons': len(rows) * 3, 'all_guide_upstream_comparisons': len(rows),
              'all_target_oracle_windows': oracle_checked, 'assignment_disagreements': 0,
              'template_prediction': {'model_sha256': model_hash, 'predicted_multiple_templates': sum(predictions),
                                      'evaluation_reads': evaluated, 'excluded_evaluation_reads': excluded,
                                      'confusion': confusion, 'truth_accuracy_claim': False},
              'synthetic_controls': controls, 'biological_inference': False,
              'files': {str(p.relative_to(stage)): sha(p) for p in sorted(stage.rglob('*')) if p.is_file()}}
    dump(stage / 'completion.json', result)
    stage.replace(output)
    print(json.dumps({k: result[k] for k in ('accession', 'statistics', 'template_prediction', 'assignment_disagreements')}, sort_keys=True))


class ScientificTests(unittest.TestCase):
    def test_complete_candidate_and_policy_oracle(self):
        for seed in range(20):
            rng = random.Random(seed)
            seqs = [''.join(rng.choice('ACGT') for _ in range(4)) for _ in range(8)]
            if seed % 2 == 0:
                seqs[-1] = seqs[0]
            index = Index(seqs)
            for letters in itertools.product('ACGTN', repeat=4):
                window = ''.join(letters)
                expected = tuple((i, sum(a != b for a, b in zip(window, s))) for i, s in enumerate(seqs)
                                 if sum(a != b for a, b in zip(window, s)) <= 1)
                self.assertEqual(index.candidates(window), expected)
                self.assertEqual(index.calls(window), oracle(seqs, window))
            self.assertEqual(index.calls(None), (-3, -3, -3))

    def test_duplicate_and_n_semantics(self):
        index = Index(['AAAA', 'AAAA', 'AAAC'])
        self.assertEqual(index.upstream('AAAA'), 1)
        self.assertEqual(index.calls('AAAA'), (-2, -2, -2))
        self.assertIsNone(index.upstream('AAAN'))
        self.assertEqual(Index(['AAAA']).calls('AAAN'), (-1, 0, 0))

    def test_conservation_can_hide_multiplicity(self):
        rows = [{'gene': 'G'}, {'gene': 'G'}]
        stats = Counter()
        add_events(stats, [(0, 0), (1, 1)], rows)
        add_events(stats, [], rows)
        add_events(stats, [], rows)
        check_conservation(stats)
        self.assertLess(stats['count_events'], stats['reads'])
        self.assertEqual(stats['extra_events'], 1)
        self.assertEqual(stats['within_gene_multiple'], 1)

    def test_real_sequence_shift_witness(self):
        a, b = 'AAGTATGCGTCAATGATGT', 'GAAGTATGCGTCAATGATG'
        self.assertEqual(sum(x != y for x, y in zip(a, b)), 17)
        index = Index([a, b])
        read = 'CTTGTGGAAAGGACGAAACACCG' + a + 'GTTTTAGA'
        self.assertEqual(counted_events(index, read, [22, 23], True), [(22, 1), (23, 0)])
        changed = read[:22] + 'A' + read[23:]
        self.assertEqual(counted_events(index, changed, [22, 23], True), [(23, 0)])
        self.assertEqual(counted_events(index, changed, [22, 23]), [(22, 1), (23, 0)])

    def test_flank_tie(self):
        self.assertEqual(mode(Counter(A=2, C=2)), 'N')
        self.assertEqual(mode(Counter(A=2, C=1)), 'A')

    def test_truncated_fastq(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'bad.gz'
            with gzip.open(path, 'wt') as f:
                f.write('@r\nAAAA\n+\n')
            with self.assertRaises(ValueError):
                list(fastq(path))


def main():
    if len(sys.argv) == 2 and sys.argv[1] == '--test':
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(ScientificTests)
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        if not result.wasSuccessful():
            raise SystemExit(1)
        print('Independent grid: 12,500 windows x three policies; all candidates also exhaustively verified.')
        return
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw', required=True)
    p.add_argument('--evidence', required=True)
    p.add_argument('--accession', choices=['ERR376998', 'ERR376999', 'SRR8297997'], required=True)
    p.add_argument('--start', type=int, required=True)
    p.add_argument('--binary', required=True)
    p.add_argument('--out', required=True)
    args = p.parse_args()
    if args.start < 0:
        p.error('start must be non-negative')
    audit(args)

if __name__ == '__main__':
    main()
