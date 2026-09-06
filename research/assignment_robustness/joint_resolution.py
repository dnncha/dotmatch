#!/usr/bin/env python3
"""AR003: joint target/position explanations, not independent offset votes.

Research-only extension. Counts at guide and gene resolution are overlapping
views and must never be added together. Bounds are conditional, not confidence
intervals. See ar003/PROTOCOL.md. No production package behaviour is changed.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import random
import resource
import sqlite3
import sys
import time
import traceback
import unittest
import replay as prior

POLICIES = ('exact', 'radius_k1', 'best_k1')
STATES = ('unique', 'ambiguous', 'none', 'invalid')
DNA = frozenset('ACGTRYSWKMBDHVN')


class Decoder:
    def __init__(self, rows, offsets, cache_size=65536):
        if not rows or any(not r.get(k) for r in rows for k in ('id', 'sequence', 'gene')):
            raise ValueError('Complete reference annotations required')
        if len({r['id'] for r in rows}) != len(rows):
            raise ValueError('Duplicate target ID')
        if not offsets or any(type(o) is not int or o < 0 for o in offsets) or len(set(offsets)) != len(offsets):
            raise ValueError('Distinct nonnegative integer offsets required')
        if type(cache_size) is not int or cache_size < 0:
            raise ValueError('Invalid cache bound')
        length = len(rows[0]['sequence'])
        if length < 2 or any(len(r['sequence']) != length or not set(r['sequence']) <= DNA for r in rows):
            raise ValueError('Uniform literal DNA targets of length >=2 required')
        self.rows = rows
        self.genes = tuple(r['gene'] for r in rows)
        self.offsets = tuple(sorted(offsets))
        self.index = prior.Index([r['sequence'] for r in rows])
        self.decode = lru_cache(maxsize=cache_size)(self._decode)
        self.fixed_calls = lru_cache(maxsize=131072)(self.index.calls)

    def resolution(self, hits):
        return (tuple(sorted({i for _, i, _ in hits})),
                tuple(sorted({self.genes[i] for _, i, _ in hits})),
                tuple(sorted({o for o, _, _ in hits})))

    def _decode(self, sequence):
        if not set(sequence) <= DNA or len(sequence) < self.offsets[-1] + self.index.length:
            return None, (None, None, None)
        hits = tuple((o, i, d) for o in self.offsets
                     for i, d in self.index.candidates(sequence[o:o+self.index.length]))
        exact = tuple(h for h in hits if h[2] == 0)
        best_distance = min((d for _, _, d in hits), default=2)
        best = tuple(h for h in hits if h[2] == best_distance)
        return hits, tuple(self.resolution(h) for h in (exact, hits, best))


def exhaustive(decoder, sequence):
    """All targets, all positions: independent of seed/index candidates."""
    if not set(sequence) <= DNA or len(sequence) < decoder.offsets[-1] + decoder.index.length:
        return None
    hits = []
    for offset in decoder.offsets:
        window = sequence[offset:offset+decoder.index.length]
        for i, target in enumerate(decoder.index.sequences):
            d = sum(a != b for a, b in zip(window, target))
            if d <= 1:
                hits.append((offset, i, d))
    return tuple(hits)


def state(values):
    return 'invalid' if values is None else 'none' if not values else 'unique' if len(values) == 1 else 'ambiguous'


class Classes:
    """Bounded write buffer; complete sufficient statistics persist on disk."""
    def __init__(self, path, limit=65536):
        self.conn = sqlite3.connect(path)
        self.conn.execute('CREATE TABLE classes(policy INTEGER, targets TEXT, n INTEGER NOT NULL CHECK(n>0), PRIMARY KEY(policy,targets)) WITHOUT ROWID')
        self.pending = Counter()
        self.limit = limit

    def add(self, policy, targets):
        self.pending[(policy, ','.join(map(str, targets))) ] += 1
        if len(self.pending) >= self.limit:
            self.flush()

    def flush(self):
        with self.conn:
            self.conn.executemany('INSERT INTO classes VALUES(?,?,?) ON CONFLICT(policy,targets) DO UPDATE SET n=n+excluded.n',
                                  ((p, t, n) for (p, t), n in self.pending.items()))
        self.pending.clear()

    def rows(self):
        self.flush()
        yield from self.conn.execute('SELECT policy,targets,n FROM classes ORDER BY policy,targets')

    def close(self):
        self.conn.close()


def verify_exhaustive_sample(decoder, sequence, hits, matrix, native_windows):
    import numpy as np
    if hits is None:
        if set(sequence) <= DNA and len(sequence) >= decoder.offsets[-1] + decoder.index.length:
            raise AssertionError('Unexpected invalid record')
        return
    expected = []
    for offset in decoder.offsets:
        window = sequence[offset:offset+decoder.index.length]
        distances = np.count_nonzero(matrix != np.frombuffer(window.encode('ascii'), dtype=np.uint8), axis=1)
        expected.extend((offset, int(i), int(distances[i])) for i in np.flatnonzero(distances <= 1))
        native_windows.add(window)
    if tuple(expected) != hits:
        raise AssertionError('All-target/all-position oracle disagreement')


def verify_native(decoder, windows):
    import dotmatch
    windows = sorted(windows)
    with dotmatch.Matcher(decoder.index.sequences) as matcher:
        outputs = matcher.assign_hamming(windows, k=1, policy='best')
    if len(outputs) != len(windows):
        raise AssertionError('Native output length mismatch')
    for window, output in zip(windows, outputs):
        hits = decoder.index.candidates(window)
        closest = min((d for _, d in hits), default=2)
        best = tuple(i for i, d in hits if d == closest)
        if dotmatch.status_name(output.status) != state(best) or output.match_count != len(hits):
            raise AssertionError('Native candidate-cardinality/status mismatch')
        if best and output.best_distance != closest:
            raise AssertionError('Native best-distance mismatch')
        if len(best) == 1 and output.target_index != best[0]:
            raise AssertionError('Native unique identity mismatch')
    return len(windows)


def controls(decoder, folder):
    """Archived balanced error-free constructs, not biological truth labels."""
    path = folder / 'synthetic/balanced/synthetic.fastq.gz'
    truth_path = folder / 'synthetic/balanced/truth.tsv'
    with truth_path.open(newline='') as f:
        truth = list(csv.DictReader(f, delimiter='\t'))
    if [(r['id'], r['gene']) for r in truth] != [(r['id'], r['gene']) for r in decoder.rows]:
        raise ValueError('Known-origin control reference order/identity mismatch')
    if any(r['known_source_reads'] != '1' for r in truth):
        raise ValueError('Expected one balanced source record per target')
    metrics = [Counter() for _ in POLICIES]
    records = 0
    for ordinal, _, sequence in prior.fastq(path):
        if ordinal > len(truth):
            raise ValueError('Too many balanced control records')
        source = ordinal - 1
        _, calls = decoder.decode(sequence)
        records += 1
        for p, call in enumerate(calls):
            metrics[p]['records'] += 1
            if call is None:
                metrics[p]['invalid'] += 1
                continue
            ids, genes, _ = call
            metrics[p]['origin_in_candidate_set'] += source in ids
            metrics[p]['guide_' + state(ids)] += 1
            metrics[p]['gene_' + state(genes)] += 1
            if len(ids) == 1:
                metrics[p]['guide_correct' if ids[0] == source else 'guide_incorrect'] += 1
            if len(genes) == 1:
                metrics[p]['gene_correct' if genes[0] == truth[source]['gene'] else 'gene_incorrect'] += 1
    if records != len(truth):
        raise ValueError('Incomplete balanced control')
    return {'kind': 'previously_archived_balanced_error_free_constructs', 'records': records,
            'fastq_sha256': prior.sha(path), 'truth_sha256': prior.sha(truth_path),
            'policies': {p: dict(metrics[i]) for i, p in enumerate(POLICIES)},
            'biological_samples': False, 'real_error_rate_or_accuracy_estimate': False}


def audit(args):
    import numpy as np
    import dotmatch
    output = Path(args.out)
    if output.exists():
        raise FileExistsError(output)
    stage = output.with_name(output.name + '.pending')
    stage.mkdir(parents=True, exist_ok=False)
    store = None
    started = time.monotonic()
    try:
        lib = prior.one(args.raw, 'library.tsv')
        reads = prior.one(args.raw, args.accession + '.fastq.gz')
        previous = prior.one(args.replay, 'completion.json')
        old = json.loads(previous.read_text())
        if old.get('completion') != 'complete' or old.get('accession') != args.accession:
            raise ValueError('Wrong or incomplete prior evidence')
        # Verify every earlier scientific artifact, not only the headline JSON.
        for name, expected in old['files'].items():
            p = previous.parent / name
            if p.resolve().is_relative_to(previous.parent.resolve()) is False or prior.sha(p) != expected:
                raise ValueError('Prior evidence hash mismatch: ' + name)
        provenance = old['input']
        for actual, expected, label in ((prior.sha(lib), old['library_sha256'], 'library'),
                                       (prior.sha(reads), provenance['sha256'], 'reads'),
                                       (prior.sha(reads, 'md5'), provenance['fastq_md5'], 'archive MD5')):
            if actual != expected:
                raise ValueError(label + ' digest mismatch')
        if reads.stat().st_size != int(provenance['fastq_bytes']):
            raise ValueError('Archive byte count mismatch')
        offsets = tuple(int(x) for x in args.offsets.split(','))
        if tuple(old['selected_offsets']) != tuple(sorted(offsets)):
            raise ValueError('Offsets differ from the already frozen discovery model')
        rows = prior.library(lib)
        decoder = Decoder(rows, offsets)
        if args.start != old['start'] or decoder.index.length != old['length']:
            raise ValueError('Fixed-window comparator differs from prior protocol')
        expected_records = int(provenance['read_count'])
        if expected_records <= 0:
            raise ValueError('Positive expected archive length required')
        selected = set(random.Random(20260906).sample(range(1, expected_records + 1), min(200, expected_records)))
        matrix = np.array([list(r['sequence'].encode('ascii')) for r in rows], dtype=np.uint8)
        native_windows = set()
        exhaustive_records = 0
        fixed_counts = [[0] * len(rows) for _ in POLICIES]
        fixed_states = [Counter() for _ in POLICIES]
        qc = [[Counter() for _ in range(3)] for _ in POLICIES]
        store = Classes(stage / 'candidate-classes.sqlite')
        examples = []
        records = 0
        prior.dump(stage / 'execution.json', {'accession': args.accession, 'input': provenance,
                   'library_sha256': old['library_sha256'], 'offsets': offsets,
                   'protocol_commit': 'de0002a37259cf27000015f43dd335343cd01234',
                   'source_sha256': prior.sha(Path(__file__)), 'index_source_sha256': prior.sha(Path(prior.__file__))})
        for ordinal, read_id, sequence in prior.fastq(reads):
            records += 1
            hits, calls = decoder.decode(sequence)
            if ordinal in selected:
                verify_exhaustive_sample(decoder, sequence, hits, matrix, native_windows)
                exhaustive_records += 1
            w = sequence[args.start:args.start+decoder.index.length] if len(sequence) >= args.start+decoder.index.length else None
            for p, call in enumerate(decoder.fixed_calls(w)):
                if call >= 0:
                    fixed_counts[p][call] += 1
                    fixed_states[p]['unique'] += 1
                else:
                    fixed_states[p][prior.STATE[call]] += 1
            for p, call in enumerate(calls):
                for r in range(3):
                    qc[p][r][state(None if call is None else call[r])] += 1
                store.add(p, () if call is None else call[0])
            if len(examples) < 50 and calls[1] is not None and len(calls[1][0]) > 1:
                examples.append((ordinal, read_id, json.dumps(hits), json.dumps(calls[1][1])))
            if records % 1000000 == 0:
                print(args.accession, records, 'records processed', flush=True)
        if records != expected_records or exhaustive_records != len(selected):
            raise AssertionError('Archive record count or exhaustive-sample coverage mismatch')
        native_checks = verify_native(decoder, native_windows)
        with (previous.parent / 'all-guide-counts.tsv').open(newline='') as f:
            table = list(csv.DictReader(f, delimiter='\t'))
        if [(r['id'], r['gene']) for r in table] != [(r['id'], r['gene']) for r in rows]:
            raise AssertionError('Native comparison reference identity/order mismatch')
        for p, policy in enumerate(POLICIES):
            if fixed_counts[p] != [int(r[policy]) for r in table]:
                raise AssertionError('Complete fixed-window/native count mismatch: ' + policy)
            if any(fixed_states[p][s] != old['native_outcomes'][policy].get(s, 0) for s in STATES):
                raise AssertionError('Complete fixed-window/native state mismatch: ' + policy)
        genes = sorted(set(decoder.genes))
        guide_counts = [Counter() for _ in POLICIES]
        bounds = [{g: [0, 0] for g in genes} for _ in POLICIES]
        class_records = Counter()
        with (stage / 'candidate-classes.tsv').open('w', newline='') as f:
            writer = csv.writer(f, delimiter='\t', lineterminator='\n')
            writer.writerow(('policy', 'target_indices', 'records'))
            for p, text, n in store.rows():
                writer.writerow((POLICIES[p], text, n))
                class_records[p] += n
                ids = tuple(int(i) for i in text.split(',') if i)
                gset = set(decoder.genes[i] for i in ids)
                if len(ids) == 1:
                    guide_counts[p][ids[0]] += n
                for gene in gset:
                    bounds[p][gene][1] += n
                    if len(gset) == 1:
                        bounds[p][gene][0] += n
        for p in range(3):
            if class_records[p] != records or any(sum(c.values()) != records for c in qc[p]):
                raise AssertionError('Complete read-state/class budget failed')
            if sum(guide_counts[p].values()) != qc[p][0]['unique']:
                raise AssertionError('Guide-unique count budget failed')
            if sum(b[0] for b in bounds[p].values()) != qc[p][1]['unique']:
                raise AssertionError('Gene-unique/lower-bound budget failed')
            if qc[p][1]['unique'] < qc[p][0]['unique']:
                raise AssertionError('Gene resolution lost a guide-unique record')
        for p in (0, 1):
            if any(guide_counts[p][i] > guide_counts[2][i] for i in range(len(rows))):
                raise AssertionError('Policy containment failed at guide level')
            if any(bounds[p][g][0] > bounds[2][g][0] for g in genes):
                raise AssertionError('Policy containment failed at gene level')
        prior.tsv(stage / 'guide-counts.tsv', ['id', 'gene', *POLICIES],
                  ((r['id'], r['gene'], *(guide_counts[p][i] for p in range(3))) for i, r in enumerate(rows)))
        prior.tsv(stage / 'gene-counts-and-bounds.tsv', ['gene', *(p+'_'+k for p in POLICIES for k in ('lower', 'upper'))],
                  ((g, *(v for p in range(3) for v in bounds[p][g])) for g in genes))
        prior.tsv(stage / 'qc.tsv', ['policy', 'resolution', *STATES],
                  ((policy, res, *(qc[p][r][s] for s in STATES)) for p, policy in enumerate(POLICIES) for r, res in enumerate(('guide', 'gene', 'position'))))
        prior.tsv(stage / 'ambiguity-examples.tsv', ['record_ordinal', 'read_id', 'candidate_hits', 'radius_genes'], examples)
        control_results = controls(decoder, previous.parent)
        prior.dump(stage / 'known-origin-controls.json', control_results)
        if prior.sha(reads) != provenance['sha256'] or prior.sha(lib) != old['library_sha256']:
            raise AssertionError('Input changed during execution')
        store.close()
        store = None
        summary = {'schema': 'dotmatch.research.ar003.v1', 'completion': 'complete', 'accession': args.accession,
                   'scope': 'full_archive', 'records': records, 'target_count': len(rows), 'gene_annotations': len(genes),
                   'offsets': offsets, 'fixed_start': args.start, 'target_length': decoder.index.length,
                   'input': provenance, 'library_sha256': old['library_sha256'],
                   'qc': {policy: {res: {s: qc[p][r][s] for s in STATES} for r, res in enumerate(('guide', 'gene', 'position'))} for p, policy in enumerate(POLICIES)},
                   'additional_gene_resolved_records': {policy: qc[p][1]['unique']-qc[p][0]['unique'] for p, policy in enumerate(POLICIES)},
                   'validation': {'full_fixed_window_native_count_cells': 3*len(rows), 'full_fixed_window_native_states_reconciled': True,
                                  'exhaustive_full_library_full_position_records': exhaustive_records, 'sampling_seed': 20260906,
                                  'native_selected_window_checks': native_checks, 'assignment_disagreements': 0},
                   'environment': {'python': sys.version, 'platform': platform.platform(), 'dotmatch': dotmatch.__version__,
                                   'wall_seconds': time.monotonic()-started, 'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                                   'source_sha256': prior.sha(Path(__file__)), 'index_source_sha256': prior.sha(Path(prior.__file__))},
                   'limitations': ['No assay-flank validation, indels, base-quality weighting or orientation inference.',
                                   'Candidate bounds assume true origin is included; they are not confidence intervals.',
                                   'Gene and guide counts overlap and cannot be added.', 'Sequencing reads are not deduplicated molecules or biological replicates.',
                                   'No biological hit calling or accuracy superiority claim.'],
                   'files': {str(p.relative_to(stage)): prior.sha(p) for p in sorted(stage.rglob('*')) if p.is_file() and p.suffix != '.sqlite'}}
        prior.dump(stage / 'completion.json', summary)
        stage.rename(output)
        print(json.dumps({'accession': args.accession, 'completion': 'complete', 'records': records}), flush=True)
    except Exception:
        prior.dump(stage / 'FAILED.json', {'completion': 'failed', 'traceback': traceback.format_exc()})
        raise
    finally:
        if store is not None:
            store.close()


class Tests(unittest.TestCase):
    def test_exhaustive_all_positions(self):
        rng = random.Random(20260906)
        for n in range(12):
            seqs = [''.join(rng.choices('ACGTN', k=4)) for _ in range(10)]
            if n % 2 == 0:
                seqs[1] = seqs[0]
            rows = [dict(id=str(i), sequence=s, gene=str(i//2)) for i, s in enumerate(seqs)]
            decoder = Decoder(rows, (0, 1, 2), cache_size=32)
            for _ in range(200):
                read = ''.join(rng.choices('ACGTN', k=6))
                hits, calls = decoder.decode(read)
                brute = exhaustive(decoder, read)
                self.assertEqual(hits, brute)
                self.assertEqual(calls[0], decoder.resolution(tuple(h for h in brute if h[2] == 0)))
                self.assertEqual(calls[1], decoder.resolution(brute))
                distance = min((h[2] for h in brute), default=2)
                self.assertEqual(calls[2], decoder.resolution(tuple(h for h in brute if h[2] == distance)))

    def test_gene_resolution_and_global_best(self):
        rows = [dict(id='a', sequence='ACGT', gene='g'), dict(id='b', sequence='CGTA', gene='g')]
        d = Decoder(rows, (0, 1))
        _, calls = d.decode('ACGTA')
        self.assertEqual(calls[0], ((0, 1), ('g',), (0, 1)))
        self.assertEqual(calls[2], calls[0])
        rows[1] = dict(id='b', sequence='CGTA', gene='other')
        self.assertEqual(Decoder(rows, (0, 1)).decode('ACGTA')[1][0][1], ('g', 'other'))

    def test_repeated_positions_not_repeated_targets(self):
        d = Decoder([dict(id='a', sequence='AAAA', gene='g')], (0, 1))
        self.assertEqual(d.decode('AAAAA')[1][0], ((0,), ('g',), (0, 1)))

    def test_invalid_and_literals(self):
        d = Decoder([dict(id='a', sequence='ACGT', gene='g')], (0, 1))
        self.assertEqual(d.decode('ACGT')[0], None)
        self.assertEqual(d.decode('ACGT?')[0], None)
        self.assertEqual(d.decode('NCGTA')[1][0][0], ())
        self.assertEqual(d.decode('NCGTA')[1][1][0], (0,))

    def test_duplicate_identity_and_offset_rejection(self):
        row = dict(id='a', sequence='ACGT', gene='g')
        for offsets in [(), (0, 0), (-1,), (True,)]:
            with self.assertRaises(ValueError):
                Decoder([row], offsets)
        with self.assertRaises(ValueError):
            Decoder([row, row], (0,))
        self.assertEqual(Decoder([row, dict(row, id='b')], (0,)).decode('ACGT')[1][0][0], (0, 1))

    def test_class_conservation_and_flush(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            c = Classes(Path(tmp)/'classes.sqlite', limit=1)
            try:
                for _ in range(3):
                    c.add(0, (1, 2))
                c.add(0, ())
                self.assertEqual(list(c.rows()), [(0, '', 1), (0, '1,2', 3)])
            finally:
                c.close()

    def test_order_and_cache_invariance(self):
        rows = [dict(id='a', sequence='ACGT', gene='g'), dict(id='b', sequence='CGTA', gene='h')]
        for read in ('ACGTA', 'TCGTA', 'NNNNN', 'ACGT'):
            a = Decoder(rows, (0, 1), cache_size=0).decode(read)
            b = Decoder(rows, (1, 0), cache_size=2).decode(read)
            self.assertEqual(a, b)
            reverse = list(reversed(rows))
            c = Decoder(reverse, (0, 1)).decode(read)
            for left, right in zip(a[1], c[1]):
                if left is None:
                    self.assertIsNone(right)
                else:
                    self.assertEqual({rows[i]['id'] for i in left[0]}, {reverse[i]['id'] for i in right[0]})
                    self.assertEqual(left[1:], right[1:])

    def test_model_violation_not_truth_bound(self):
        rows = [dict(id='a', sequence='AAAA', gene='a'), dict(id='b', sequence='AACC', gene='b')]
        _, calls = Decoder(rows, (0,)).decode('AACC')
        # A two-substitution A-origin read can be a unique B call. Do not call
        # a narrow candidate interval an unconditional biological guarantee.
        self.assertEqual(calls[1][0], (1,))
        self.assertNotIn(0, calls[1][0])


def main():
    if '--test' in sys.argv:
        unittest.main(argv=[sys.argv[0]], verbosity=2)
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', required=True)
    parser.add_argument('--replay', required=True)
    parser.add_argument('--accession', required=True, choices=('ERR376998', 'ERR376999', 'SRR8297997'))
    parser.add_argument('--offsets', required=True)
    parser.add_argument('--start', type=int, required=True)
    parser.add_argument('--out', required=True)
    audit(parser.parse_args())


if __name__ == '__main__':
    main()
