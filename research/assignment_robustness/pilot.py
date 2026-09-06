#!/usr/bin/env python3
"""AR-001 exploratory pilot. Standard library only; DotMatch is an external CLI.

No synthetic fallback, no gene-level inference, no production code modifications.
Run `python pilot.py test` before `python pilot.py run --out NEW_DIRECTORY`.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import itertools
import json
import platform
import random
import subprocess
import sys
import time
import unittest
import urllib.parse
import urllib.request
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
POLICIES = ('exact', 'radius_k1', 'best_k1')
STATES = ('unique', 'ambiguous', 'none', 'invalid')
PAIRS = tuple(itertools.combinations(range(3), 2))
BASE = '11d159fa1648365f2a4e96917b483c33aa5d9fe7'
SOURCES = {
    'yusa': {
        'length': 19,
        'library_url': 'https://sourceforge.net/projects/mageck/files/libraries/yusa_library.csv.zip/download',
        'member': 'yusa_library.csv',
        'samples': {'plasmid': ['ERR376998'], 'ESC1': ['ERR376999']},
    },
    'brunello': {
        'length': 20,
        'library_url': 'https://sourceforge.net/projects/mageck/files/libraries/broadgpp-brunello-library-corrected.txt.zip/download',
        'member': 'broadgpp-brunello-library-corrected.txt',
        'samples': {'plasmid': ['SRR8297997'], 'RepA': ['SRR8297837', 'SRR8297836'],
                    'RepB': ['SRR8297839', 'SRR8297838'], 'RepC': ['SRR8297841', 'SRR8297840']},
    },
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def save_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + '\n', encoding='utf-8')


def table(path, header, rows):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)


def fetch_bytes(url, limit=32 * 1024 * 1024):
    request = urllib.request.Request(url, headers={'User-Agent': 'DotMatch-AR001-public-research/1.0'})
    with urllib.request.urlopen(request, timeout=120) as handle:
        data = handle.read(limit + 1)
    require(len(data) <= limit, 'download exceeds explicit byte limit')
    return data


def metadata(accession):
    fields = 'run_accession,sample_accession,sample_alias,experiment_title,read_count,fastq_ftp,fastq_bytes,fastq_md5'
    url = 'https://www.ebi.ac.uk/ena/portal/api/filereport?' + urllib.parse.urlencode({
        'accession': accession, 'result': 'read_run', 'fields': fields, 'format': 'tsv'})
    raw = fetch_bytes(url, 1024 * 1024)
    rows = list(csv.DictReader(io.StringIO(raw.decode('utf-8')), delimiter='\t'))
    require(len(rows) == 1 and rows[0]['run_accession'] == accession, 'ENA accession mismatch')
    row = rows[0]
    paths = row['fastq_ftp'].split(';')
    require(len(paths) == 1 and paths[0], 'pilot requires one single-end FASTQ per run')
    row['url'] = 'https://' + paths[0].removeprefix('ftp://')
    require(urllib.parse.urlsplit(row['url']).hostname == 'ftp.sra.ebi.ac.uk', 'unexpected ENA file host')
    row['metadata_url'] = url
    row['metadata_sha256'] = hashlib.sha256(raw).hexdigest()
    row['archive_md5_locally_verified'] = False
    return row


def records(handle):
    """Strict four-line FASTQ reader; rejects partial records and invalid quality."""
    ordinal = 0
    while True:
        header = handle.readline()
        if not header:
            return
        ordinal += 1
        sequence, plus, quality = (handle.readline() for _ in range(3))
        require(sequence and plus and quality, f'truncated FASTQ record {ordinal}')
        lines = [line.rstrip(b'\r\n').decode('ascii', errors='strict') for line in (header, sequence, plus, quality)]
        h, s, p, q = lines
        require(h.startswith('@') and len(h) > 1 and p.startswith('+'), f'invalid FASTQ framing {ordinal}')
        require(s and len(s) == len(q), f'sequence/quality length mismatch {ordinal}')
        require(all(33 <= ord(c) <= 126 for c in q), f'invalid quality {ordinal}')
        require(all(c in 'ACGTRYSWKMBDHVN' for c in s), f'unsupported sequence alphabet {ordinal}')
        yield lines


def write_fastq(path, recs):
    with Path(path).open('wb') as raw:
        with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0, filename='') as handle:
            for rec in recs:
                handle.write(('\n'.join(rec) + '\n').encode('ascii'))


def fetch_prefix(accessions, n):
    recs, provenance = [], []
    for accession in accessions:
        meta = metadata(accession)
        meta['contributed_prefix_records'] = 0
        if len(recs) < n:
            with urllib.request.urlopen(meta['url'], timeout=120) as remote:
                with gzip.GzipFile(fileobj=remote, mode='rb') as stream:
                    for rec in itertools.islice(records(stream), n - len(recs)):
                        recs.append(rec)
                        meta['contributed_prefix_records'] += 1
        provenance.append(meta)
    require(len(recs) == n, f'expected {n} prefix records, obtained {len(recs)}')
    return recs, provenance


def parse_library(text, length):
    require(text.strip(), 'empty reference library')
    delimiter = '\t' if '\t' in text.splitlines()[0] else ','
    rows = list(csv.reader(io.StringIO(text), delimiter=delimiter))
    require(rows and all(len(row) == 3 for row in rows if row), 'library must have three columns')
    rows = [row for row in rows if row]
    first = rows[0]
    seq_names = {'seq', 'sequence', 'grna.sequence', 'sgrna_sequence', 'guide_seq', 'sgrna sequence'}
    if first[1].strip().lower() in seq_names:
        rows = rows[1:]
    require(rows, 'empty reference library')
    seen = set()
    result = []
    for ident, sequence, gene in rows:
        ident, sequence, gene = ident.strip(), sequence.strip().upper(), gene.strip()
        require(ident and ident not in seen, f'duplicate or empty guide ID: {ident!r}')
        require(len(sequence) == length and set(sequence) <= set('ACGT'), 'reference length/alphabet mismatch')
        require(gene, f'empty gene annotation for {ident}')
        seen.add(ident)
        result.append((ident, sequence, gene))
    return result


def fetch_library(name, work):
    config = SOURCES[name]
    archive = fetch_bytes(config['library_url'])
    with zipfile.ZipFile(io.BytesIO(archive)) as zf:
        members = [x for x in zf.infolist() if x.filename.split('/')[-1] == config['member']]
        require(len(members) == 1 and members[0].file_size < 32 * 1024 * 1024, 'invalid library archive')
        raw = zf.read(members[0])
    rows = parse_library(raw.decode('utf-8-sig'), config['length'])
    path = work / f'{name}.targets.tsv'
    table(path, ['sgRNA', 'sequence', 'Gene'], rows)
    return rows, path, {'url': config['library_url'], 'member': config['member'], 'rows': len(rows),
        'source_sha256': hashlib.sha256(raw).hexdigest(), 'archive_sha256': hashlib.sha256(archive).hexdigest(),
        'normalized_sha256': digest(path)}


class Reference:
    """Independent candidate enumeration: no DotMatch imports or index code."""
    def __init__(self, library):
        self.library = library
        self.length = len(library[0][1])
        self.lookup = defaultdict(list)
        for i, (_ident, sequence, _gene) in enumerate(library):
            self.lookup[sequence].append(i)

    def candidates(self, window):
        if window is None:
            return [], []
        require(len(window) == self.length, 'wrong query length')
        exact = list(self.lookup.get(window, []))
        near = []
        for pos, old in enumerate(window):
            for base in 'ACGT':
                if base != old:
                    near.extend(self.lookup.get(window[:pos] + base + window[pos + 1:], []))
        return exact, near

    @staticmethod
    def unique(indices):
        return ('unique', indices[0]) if len(indices) == 1 else ('ambiguous', -1) if indices else ('none', -1)

    def calls(self, window):
        if window is None:
            return (('invalid', -1),) * 3
        exact, near = self.candidates(window)
        return self.unique(exact), self.unique(exact + near), self.unique(exact if exact else near)

    def brute(self, window):
        if window is None:
            return (('invalid', -1),) * 3
        exact, near = [], []
        for i, (_ident, sequence, _gene) in enumerate(self.library):
            d = 0
            for a, b in zip(sequence, window):
                d += a != b
                if d > 1:
                    break
            if d == 0:
                exact.append(i)
            elif d == 1:
                near.append(i)
        return self.unique(exact), self.unique(exact + near), self.unique(exact if exact else near)


def geometry(ref, out, name):
    edges, near_ids, duplicate_ids = [], set(), set()
    duplicate_groups = 0
    for sequence, indices in ref.lookup.items():
        if len(indices) > 1:
            duplicate_groups += 1
            duplicate_ids.update(indices)
        _exact, near = ref.candidates(sequence)
        for left in indices:
            for right in near:
                if left < right:
                    near_ids.update((left, right))
                    a, b = ref.library[left], ref.library[right]
                    edges.append((a[0], b[0], a[2], b[2], int(a[2] == b[2])))
    table(out / f'{name}.distance_one_edges.tsv', ['guide_a', 'guide_b', 'gene_a', 'gene_b', 'same_gene'], edges)
    return {'guide_rows': len(ref.library), 'distinct_sequences': len(ref.lookup),
        'duplicate_sequence_groups': duplicate_groups, 'guides_in_duplicate_groups': len(duplicate_ids),
        'distance_one_pairs': len(edges), 'guides_with_distance_one_neighbour': len(near_ids),
        'cross_gene_distance_one_pairs': sum(1 for x in edges if not x[-1])}


def offsets(recs, ref):
    distribution = Counter({i: 0 for i in range(41)})
    for _h, seq, _p, _q in recs:
        for start in range(41):
            if seq[start:start + ref.length] in ref.lookup:
                distribution[start] += 1
    return distribution


def evaluate(label, recs, start, ref, library_path, out, work, command):
    sample_out = out / label
    sample_out.mkdir()
    reads_path = work / f'{label}.evaluation.fastq.gz'
    write_fastq(reads_path, recs)
    cli = [command, 'sensitivity', '--targets', str(library_path), '--reads', str(reads_path),
           '--target-start', str(start), '--target-length', str(ref.length), '--sample-label', label,
           '--write-read-changes', '--out-dir', str(sample_out / 'dotmatch')]
    started = time.perf_counter()
    proc = subprocess.run(cli, capture_output=True, text=True, timeout=600)
    (sample_out / 'command.stdout.txt').write_text(proc.stdout)
    (sample_out / 'command.stderr.txt').write_text(proc.stderr)
    save_json(sample_out / 'command.json', {'argv': cli, 'returncode': proc.returncode,
              'seconds': time.perf_counter() - started})
    require(proc.returncode == 0, f'DotMatch failed for {label}: {proc.stderr[-2000:]}')
    summary = json.loads((sample_out / 'dotmatch' / 'summary.json').read_text())
    require(summary['completion'] == 'complete', 'incomplete DotMatch output')
    require(summary['read_count'] == len(recs), 'DotMatch read denominator differs')
    require(summary['inputs']['reads']['sha256'] == digest(reads_path), 'DotMatch read digest differs')
    require(summary['inputs']['targets']['sha256'] == digest(library_path), 'DotMatch reference digest differs')
    for filename, info in summary['artifacts'].items():
        require(digest(sample_out / 'dotmatch' / filename) == info['sha256'], f'DotMatch output digest differs: {filename}')
    counts = [Counter() for _ in POLICIES]
    states = [Counter({state: 0 for state in STATES}) for _ in POLICIES]
    changes = []
    pair_changes = Counter()
    mechanisms = Counter()
    transitions = [Counter() for _ in PAIRS]
    cache, brute_cases = {}, defaultdict(list)
    for ordinal, (_h, seq, _p, _q) in enumerate(recs, 1):
        window = seq[start:start + ref.length] if len(seq) >= start + ref.length else None
        if window not in cache:
            calls = ref.calls(window)
            cache[window] = calls
            stratum = '|'.join(call[0] for call in calls)
            if len(brute_cases[stratum]) < 4:
                brute_cases[stratum].append(window)
            if window and 'N' in window and len(brute_cases['contains_N']) < 4:
                brute_cases['contains_N'].append(window)
        calls = cache[window]
        for i, (status, target) in enumerate(calls):
            states[i][status] += 1
            if status == 'unique':
                counts[i][target] += 1
        for j, (a, b) in enumerate(PAIRS):
            transitions[j][(calls[a][0], calls[b][0])] += 1
            if calls[a] != calls[b]:
                pair_changes[f'{POLICIES[a]}__{POLICIES[b]}'] += 1
        if len(set(calls)) > 1:
            values = [str(ordinal)]
            for status, target in calls:
                values.extend((status, ref.library[target][0] if status == 'unique' else ''))
            changes.append(values)
            if calls[0][0] == 'none' and calls[2][0] == 'unique':
                mechanisms['nonexact_unique_nearest'] += 1
            if calls[0][0] == 'unique' and calls[1][0] == 'ambiguous':
                mechanisms['exact_hit_with_distance_one_alternative'] += 1
            if calls[0][0] == 'none' and calls[2][0] == 'ambiguous':
                mechanisms['nonexact_nearest_tie'] += 1
    checked = set()
    for cases in brute_cases.values():
        for window in cases:
            if window not in checked:
                require(ref.brute(window) == cache[window], 'independent oracle failed exhaustive all-target check')
                checked.add(window)
    for i, policy in enumerate(POLICIES):
        require(sum(states[i].values()) == len(recs), 'oracle read conservation failed')
        require(sum(counts[i].values()) == states[i]['unique'], 'oracle count conservation failed')
        require(all(summary['outcomes'][policy].get(s, 0) == states[i][s] for s in STATES), 'read-state totals differ')
        with (sample_out / 'dotmatch' / f'{policy}.counts.tsv').open() as handle:
            reader = csv.reader(handle, delimiter='\t')
            require(next(reader) == ['sgRNA', 'Gene', label], 'sample axis differs')
            observed = list(reader)
        require(len(observed) == len(ref.library), 'guide axis length differs')
        for j, row in enumerate(observed):
            require(row == [ref.library[j][0], ref.library[j][2], str(counts[i][j])], f'guide count disagreement: {policy}/{j}')
    with (sample_out / 'dotmatch' / 'read_changes.tsv').open() as handle:
        reader = csv.reader(handle, delimiter='\t')
        next(reader)
        observed_changes = [[row[0], *row[2:]] for row in reader]
    require(observed_changes == changes, 'changed-read calls or ordinals disagree')
    require(summary['changed_reads'] == len(changes), 'changed-read total differs')
    with (sample_out / 'dotmatch' / 'transitions.tsv').open() as handle:
        observed_transitions = list(csv.DictReader(handle, delimiter='\t'))
    require(len(observed_transitions) == 48, 'transition table incomplete')
    for row in observed_transitions:
        pair = (POLICIES.index(row['from_policy']), POLICIES.index(row['to_policy']))
        expected = transitions[PAIRS.index(pair)][(row['from_status'], row['to_status'])]
        require(int(row['reads']) == expected, 'transition count disagrees')
    result = {'sample': label, 'scope': 'discovery-excluded archival prefix; not whole-screen prevalence',
        'records': len(recs), 'target_start': start, 'target_length': ref.length,
        'read_sha256': digest(reads_path), 'states': dict(zip(POLICIES, states)),
        'changed_records': len(changes), 'changed_fraction': len(changes) / len(recs),
        'pair_changed_records': {f'{POLICIES[a]}__{POLICIES[b]}': pair_changes[f'{POLICIES[a]}__{POLICIES[b]}'] for a, b in PAIRS},
        'mechanisms': dict(mechanisms), 'unique_windows': len(cache),
        'all_reference_bruteforce_windows': len(checked),
        'independent_count_state_transition_changed_call_agreement': True,
        'dotmatch': summary}
    save_json(sample_out / 'validation.json', result)
    print(json.dumps({k: result[k] for k in ('sample', 'records', 'target_start', 'changed_records', 'states', 'mechanisms')}), flush=True)
    return result


def run(args):
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    work = out / 'work'
    work.mkdir()
    protocol_hash = digest(HERE / 'PROTOCOL.md')
    summary = {'schema': 'dotmatch.research.AR001.v1', 'status': 'running', 'baseline': BASE,
               'harness_sha256': digest(__file__), 'protocol_sha256': protocol_hash,
               'python': sys.version, 'platform': platform.platform(), 'libraries': {}, 'samples': [],
               'read_evaluation_records': 100000, 'discovery_records': 2000,
               'limitations': ['Exploratory selected pilot, not a registered confirmatory study.',
                   'Nonrandom prefixes; no full-archive checksum validation.', 'No gene-level inference.',
                   'Brunello mixed offsets may invalidate fixed-window use as a complete assay workflow.',
                   'No biological accuracy or competitor superiority claim.']}
    summary['git_commit'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    summary['dotmatch_version'] = subprocess.check_output([args.dotmatch, '--version'], text=True).strip()
    save_json(out / 'run.json', summary)
    try:
        for name, config in SOURCES.items():
            library, lib_path, provenance = fetch_library(name, work)
            ref = Reference(library)
            summary['libraries'][name] = {'provenance': provenance, 'geometry': geometry(ref, out, name)}
            print(json.dumps({'library': name, **summary['libraries'][name]}), flush=True)
            for sample, accessions in config['samples'].items():
                label = f'{name}.{sample}'
                recs, metas = fetch_prefix(accessions, 102000)
                discovery = recs[:2000]
                dist = offsets(discovery, ref)
                start = 23 if name == 'yusa' else min(dist, key=lambda k: (-dist[k], k))
                require(dist[start] > 0, 'no discovery exact hits at selected offset')
                write_fastq(work / f'{label}.discovery.fastq.gz', discovery)
                write_fastq(work / f'{label}.prefix.fastq.gz', recs)
                source = {'sample': label, 'sources': metas, 'offset_discovery_exact_counts': dict(dist),
                          'selected_start': start, 'selection': 'documented' if name == 'yusa' else 'modal_discovery_exact_offset',
                          'discovery_sha256': digest(work / f'{label}.discovery.fastq.gz'),
                          'whole_prefix_sha256': digest(work / f'{label}.prefix.fastq.gz'),
                          'evaluation_first_source_ordinal': 2001, 'evaluation_records': 100000}
                save_json(out / f'{label}.source.json', source)
                result = evaluate(label, recs[2000:], start, ref, lib_path, out, work, args.dotmatch)
                summary['samples'].append({k: v for k, v in result.items() if k != 'dotmatch'})
                save_json(out / 'run.json', summary)
        require(digest(HERE / 'PROTOCOL.md') == protocol_hash, 'protocol changed during run')
        summary['status'] = 'complete'
    except Exception as exc:
        summary['status'] = 'failed'
        summary['error'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        save_json(out / 'run.json', summary)
        files = {str(p.relative_to(out)): {'bytes': p.stat().st_size, 'sha256': digest(p)}
                 for p in sorted(out.rglob('*')) if p.is_file() and p.name != 'MANIFEST.json' and 'work' not in p.relative_to(out).parts}
        save_json(out / 'MANIFEST.json', {'status': summary['status'], 'files': files})


class Tests(unittest.TestCase):
    def test_policy_semantics(self):
        ref = Reference([('a', 'AAAA', 'g1'), ('b', 'AAAC', 'g2'), ('c', 'CCCC', 'g3'), ('d', 'GGGG', 'g4'), ('e', 'GGGG', 'g5')])
        self.assertEqual(ref.calls('AAAA'), (('unique', 0), ('ambiguous', -1), ('unique', 0)))
        self.assertEqual(ref.calls('AAAG'), (('none', -1), ('ambiguous', -1), ('ambiguous', -1)))
        self.assertEqual(ref.calls('CCCT'), (('none', -1), ('unique', 2), ('unique', 2)))
        self.assertEqual(ref.calls('GGGG'), (('ambiguous', -1),) * 3)
        self.assertEqual(ref.calls('TTTT'), (('none', -1),) * 3)
        self.assertEqual(ref.calls(None), (('invalid', -1),) * 3)
        self.assertEqual(ref.calls('CCCN'), (('none', -1), ('unique', 2), ('unique', 2)))
        self.assertEqual(ref.calls('CCNN'), (('none', -1),) * 3)

    def test_exhaustive_all_queries(self):
        rng = random.Random(20260906)
        sequences = [''.join(x) for x in itertools.product('ACGT', repeat=4)]
        for trial in range(5):
            selected = rng.sample(sequences, 20) + ['AAAA', 'AAAA']
            ref = Reference([(str(i), s, str(i)) for i, s in enumerate(selected)])
            for bases in itertools.product('ACGTN', repeat=4):
                q = ''.join(bases)
                self.assertEqual(ref.calls(q), ref.brute(q), (trial, q))

    def test_library(self):
        self.assertEqual(len(parse_library('id,gRNA.sequence,Gene\na,AAAA,g\nb,AAAA,h\n', 4)), 2)
        for text in ('', 'a\tAAAA\tg\na\tCCCC\th\n', 'a\tAAAN\tg\n', 'a\tAAA\tg\n', 'a\tAAAA\t\n'):
            with self.assertRaises(ValueError):
                parse_library(text, 4)

    def test_brunello_source_header(self):
        rows = parse_library('sgRNAID\tSeq\tgene\ng1\tACGTACGTACGTACGTACGT\tGENE1\n', 20)
        self.assertEqual(rows, [('g1', 'ACGTACGTACGTACGTACGT', 'GENE1')])

    def test_fastq(self):
        self.assertEqual(len(list(records(io.BytesIO(b'@a\nACGT\n+\nIIII\n')))), 1)
        for raw in (b'@a\nACGT\n+\n', b'@a\nACGT\n+\nIII\n', b'a\nACGT\n+\nIIII\n', b'@a\nACGT\n-\nIIII\n'):
            with self.assertRaises(ValueError):
                list(records(io.BytesIO(raw)))

    def test_geometry(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            ref = Reference([('a', 'AAAA', 'g'), ('b', 'AAAC', 'h'), ('c', 'AAAA', 'g')])
            g = geometry(ref, Path(tmp), 'test')
            self.assertEqual(g['distance_one_pairs'], 2)
            self.assertEqual(g['guides_in_duplicate_groups'], 2)
            self.assertEqual(g['cross_gene_distance_one_pairs'], 2)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        unittest.main(argv=[sys.argv[0]], verbosity=2)
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('action', choices=['run'])
        parser.add_argument('--out', required=True)
        parser.add_argument('--dotmatch', default='dotmatch')
        run(parser.parse_args())
