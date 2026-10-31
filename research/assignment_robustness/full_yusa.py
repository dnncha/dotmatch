#!/usr/bin/env python3
"""AR-001-F: full checksummed Yusa archives, no biological hit calling.

Independent fast decoding and event reconstruction are checked against DotMatch
and guide-counter. Read records, offset events and biological molecules are not
interchangeable units. See PROTOCOL.md and AMENDMENTS.md.
"""
from __future__ import annotations
import argparse
import csv
import functools
import gzip
import hashlib
import itertools
import json
import platform
import random
import shutil
import subprocess
import sys
import time
import unittest
import urllib.request
from collections import Counter
from pathlib import Path
import pilot
import forensics

NONE, AMBIGUOUS, INVALID = -1, -2, -3
STATE = {NONE: 'none', AMBIGUOUS: 'ambiguous', INVALID: 'invalid'}


def download_full(accession: str, work: Path) -> tuple[Path, dict]:
    meta = pilot.metadata(accession)
    expected_bytes, expected_reads = int(meta['fastq_bytes']), int(meta['read_count'])
    pilot.require(0 < expected_bytes < 2_000_000_000 and expected_reads > 0, 'archive outside bounded study scope')
    path = work / f'{accession}.fastq.gz'
    tmp = path.with_suffix('.partial')
    h, md5, size = hashlib.sha256(), hashlib.md5(), 0
    started = time.perf_counter()
    try:
        with urllib.request.urlopen(meta['url'], timeout=120) as remote, tmp.open('xb') as out:
            while True:
                chunk = remote.read(1024 * 1024)
                if not chunk:
                    break
                size += len(chunk)
                pilot.require(size <= expected_bytes, 'archive larger than ENA metadata')
                h.update(chunk)
                md5.update(chunk)
                out.write(chunk)
        pilot.require(size == expected_bytes, 'archive byte count mismatch')
        pilot.require(md5.hexdigest() == meta['fastq_md5'].lower(), 'full archive MD5 mismatch')
        tmp.rename(path)
    finally:
        tmp.unlink(missing_ok=True)
    meta.update(local_bytes=size, local_sha256=h.hexdigest(), local_md5=md5.hexdigest(),
                archive_md5_locally_verified=True, download_seconds=time.perf_counter() - started)
    print(json.dumps({'full_download_verified': accession, 'bytes': size, 'expected_reads': expected_reads,
                      'sha256': h.hexdigest(), 'md5': md5.hexdigest()}), flush=True)
    return path, meta


class Decoder:
    """Separate codeword lookup, checked against query enumeration and brute force."""
    def __init__(self, library):
        self.ref = pilot.Reference(library)
        self.lookup = forensics.matching_event_lookup(library)
        self.exact = {row[1]: i for i, row in enumerate(library)}
        self.with_neighbour = {i for i, row in enumerate(library) if self.ref.candidates(row[1])[1]}

    @functools.lru_cache(maxsize=131072)
    def calls(self, window: str | None) -> tuple[int, int, int]:
        if window is None:
            return (INVALID,) * 3
        pilot.require(len(window) == self.ref.length, 'query length mismatch')
        exact = self.exact.get(window)
        if exact is not None:
            return exact, AMBIGUOUS if exact in self.with_neighbour else exact, exact
        value = self.lookup.get(window)
        if value is not None:
            nearest = value if value >= 0 else AMBIGUOUS
            return NONE, nearest, nearest
        if set(window) - set('ACGT'):
            return tuple(index if state == 'unique' else {'none': NONE, 'ambiguous': AMBIGUOUS, 'invalid': INVALID}[state]
                         for state, index in self.ref.calls(window))
        return (NONE,) * 3

    def check(self, window):
        observed = tuple(('unique', c) if c >= 0 else (STATE[c], -1) for c in self.calls(window))
        pilot.require(observed == self.ref.calls(window), 'fast decoder differs from independent query enumeration')


def calibration(path: Path, decoder: Decoder) -> tuple[list[int], Counter]:
    n, found = 0, Counter()
    with gzip.open(path, 'rb') as handle:
        for _h, sequence, _p, _q in itertools.islice(pilot.records(handle), 100000):
            n += 1
            for offset in range(len(sequence) - decoder.ref.length + 1):
                if decoder.lookup.get(sequence[offset:offset + decoder.ref.length], -1) >= 0:
                    found[offset] += 1
    pilot.require(n == 100000 and sum(found.values()) > 0, 'insufficient comparator calibration records')
    total = sum(found.values())
    selected = sorted(offset for offset, count in found.items() if count * 400 >= total)
    return selected, found


def independent_sample(path: Path, label: str, decoder: Decoder, metadata: dict, out: Path) -> dict:
    selected, offsets = calibration(path, decoder)
    pilot.table(out / f'{label}.calibration_offsets.tsv', ['offset', 'matching_events', 'selected'],
                [(k, v, int(k in selected)) for k, v in sorted(offsets.items())])
    n_guides = len(decoder.ref.library)
    counts = [[0] * n_guides for _ in range(3)]
    states = [Counter({s: 0 for s in pilot.STATES}) for _ in range(3)]
    event_counts = [0] * n_guides
    multiplicity, mechanisms, pairs = Counter(), Counter(), Counter()
    traces, brute_checked = [], set()
    same_guide = different_guide = cross_gene = same_gene = read_any = changed = n = 0
    fixed_only = events_only = both = neither = 0
    started = time.perf_counter()
    with gzip.open(path, 'rb') as handle:
        for n, (_h, sequence, _p, _q) in enumerate(pilot.records(handle), 1):
            window = sequence[23:42] if len(sequence) >= 42 else None
            calls = decoder.calls(window)
            if n <= 1000 or n % 100000 == 0:
                decoder.check(window)
            stratum = tuple('unique' if c >= 0 else STATE[c] for c in calls)
            if stratum not in brute_checked:
                expected = decoder.ref.brute(window)
                observed = tuple(('unique', c) if c >= 0 else (STATE[c], -1) for c in calls)
                pilot.require(expected == observed, 'all-reference brute-force disagreement')
                brute_checked.add(stratum)
            for policy, call in enumerate(calls):
                if call >= 0:
                    counts[policy][call] += 1
                    states[policy]['unique'] += 1
                else:
                    states[policy][STATE[call]] += 1
            if not calls[0] == calls[1] == calls[2]:
                changed += 1
            for a, b in pilot.PAIRS:
                if calls[a] != calls[b]:
                    pairs[f'{pilot.POLICIES[a]}__{pilot.POLICIES[b]}'] += 1
            if calls[0] == NONE and calls[2] >= 0:
                mechanisms['nonexact_unique_nearest'] += 1
            if calls[0] >= 0 and calls[1] == AMBIGUOUS:
                mechanisms['exact_hit_with_distance_one_alternative'] += 1
            if calls[0] == NONE and calls[2] == AMBIGUOUS:
                mechanisms['nonexact_nearest_tie'] += 1
            events = []
            for offset in selected:
                if offset + 19 <= len(sequence):
                    target = decoder.lookup.get(sequence[offset:offset + 19], -1)
                    if target >= 0:
                        event_counts[target] += 1
                        events.append((offset, target))
            multiplicity[len(events)] += 1
            if events:
                read_any += 1
                if calls[2] >= 0:
                    both += 1
                else:
                    events_only += 1
            elif calls[2] >= 0:
                fixed_only += 1
            else:
                neither += 1
            if len(events) > 1:
                ids = {target for _, target in events}
                genes = {decoder.ref.library[target][2] for target in ids}
                same_guide += len(ids) == 1
                different_guide += len(ids) > 1
                cross_gene += len(genes) > 1
                same_gene += len(genes) == 1
                if len(traces) < 1000:
                    traces.append([n, ';'.join(f'{offset}:{decoder.ref.library[target][0]}' for offset, target in events),
                                   len(ids), len(genes)])
            if n % 2000000 == 0:
                print(json.dumps({'sample': label, 'independently_checked_records': n}), flush=True)
    pilot.require(n == int(metadata['read_count']), 'full parsed read denominator differs from ENA')
    for i in range(3):
        pilot.require(sum(states[i].values()) == n and sum(counts[i]) == states[i]['unique'], 'full read conservation failed')
    pilot.require(sum(event_counts) == sum(k * v for k, v in multiplicity.items()), 'full event conservation failed')
    pilot.require(fixed_only + both == states[2]['unique'] and events_only + both == read_any, 'membership reconciliation failed')
    for policy, vector in zip(pilot.POLICIES, counts):
        pilot.table(out / f'{label}.oracle.{policy}.tsv', ['sgRNA', 'Gene', label],
                    [(row[0], row[2], vector[i]) for i, row in enumerate(decoder.ref.library)])
    pilot.table(out / f'{label}.oracle.events.tsv', ['sgRNA', 'Gene', label],
                [(row[0], row[2], event_counts[i]) for i, row in enumerate(decoder.ref.library)])
    pilot.table(out / f'{label}.first1000_multi_offset_records.tsv', ['record_ordinal', 'offset_and_guide', 'distinct_guides', 'distinct_gene_annotations'], traces)
    result = {'sample': label, 'records': n, 'states': dict(zip(pilot.POLICIES, states)), 'changed_records': changed,
        'changed_fraction': changed / n, 'pair_changed_records': dict(pairs), 'mechanisms': dict(mechanisms),
        'brute_force_all_reference_strata': len(brute_checked), 'selected_offsets': selected,
        'event_multiplicity': dict(multiplicity), 'matching_events': sum(event_counts),
        'reads_with_any_event': read_any, 'excess_events_above_matched_reads': sum(event_counts) - read_any,
        'multi_offset_same_guide_reads': same_guide, 'multi_offset_different_guide_reads': different_guide,
        'multi_offset_same_gene_reads': same_gene, 'multi_offset_cross_gene_reads': cross_gene,
        'membership_vs_fixed_best': {'both': both, 'event_only': events_only, 'fixed_only': fixed_only, 'neither': neither},
        'trace_selection': 'first at most 1000 multi-event records; not a random sample',
        'seconds': time.perf_counter() - started}
    pilot.save_json(out / f'{label}.independent.json', result)
    print(json.dumps(result), flush=True)
    return result


def compare_matrices(left: Path, right: Path, sample: str) -> None:
    ls, a = forensics.matrix(left)
    rs, b = forensics.matrix(right)
    pilot.require(sample in ls and sample in rs and a.keys() == b.keys(), 'guide or sample axes differ')
    bad = [(key, a[key], b[key]) for key in a if a[key][0] != b[key][0] or a[key][1][sample] != b[key][1][sample]]
    pilot.require(not bad, f'{left.name} vs {right.name}: {len(bad)} guide discrepancies; first: {bad[:2]}')


def run(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=False)
    work = out / 'work'
    work.mkdir()
    summary = {'schema': 'dotmatch.research.AR001F.v1', 'status': 'running', 'baseline': pilot.BASE,
        'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=pilot.ROOT, text=True).strip(),
        'python': sys.version, 'platform': platform.platform(),
        'code_sha256': {p.name: pilot.digest(p) for p in [Path(__file__), Path(pilot.__file__), Path(forensics.__file__)]},
        'protocol_sha256': pilot.digest(pilot.HERE / 'PROTOCOL.md'),
        'amendments_sha256': pilot.digest(pilot.HERE / 'AMENDMENTS.md'),
        'samples': {}, 'results': {}, 'limitations': ['One plasmid and one cellular library, not replicated biological contrasts.',
            'No gene-level significance, false-assignment rate or biological accuracy claim.',
            'Reference sequence and gene annotations do not establish molecular origin.']}
    pilot.save_json(out / 'run.json', summary)
    try:
        library, lib_path, provenance = pilot.fetch_library('yusa', work)
        summary['library'] = provenance
        paths = []
        labels = list(pilot.SOURCES['yusa']['samples'])
        for label, accessions in pilot.SOURCES['yusa']['samples'].items():
            pilot.require(len(accessions) == 1, 'full Yusa requires one source per sample')
            path, metadata = download_full(accessions[0], work)
            paths.append(path)
            summary['samples'][label] = metadata
            pilot.save_json(out / 'run.json', summary)
        installed = subprocess.check_output(['cargo', 'install', '--list'], text=True)
        pilot.require('guide-counter v0.1.3:' in installed.splitlines(), 'wrong comparator version')
        summary['guide_counter_version_evidence'] = installed
        summary['guide_counter_binary_sha256'] = pilot.digest(shutil.which('guide-counter'))
        summary['dotmatch_version'] = subprocess.check_output(['dotmatch', '--version'], text=True).strip()
        summary['crate_hashes'] = {}
        for crate in (Path.home() / '.cargo/registry/src').glob('*/guide-counter-0.1.3'):
            for rel in ('Cargo.lock', '.cargo_vcs_info.json', 'src/commands/count.rs', 'src/guide.rs'):
                if (crate / rel).exists():
                    summary['crate_hashes'][rel] = pilot.digest(crate / rel)
        gc = ['guide-counter', 'count', '--input', *map(str, paths), '--samples', *labels,
            '--library', str(lib_path), '--output', str(out / 'guide_counter'),
            '--offset-sample-size', '100000', '--offset-min-fraction', '0.0025']
        forensics.command(gc, out / 'guide_counter')
        for label, path in zip(labels, paths):
            argv = ['dotmatch', 'sensitivity', '--targets', str(lib_path), '--reads', str(path),
                '--target-start', '23', '--target-length', '19', '--sample-label', label,
                '--out-dir', str(out / f'{label}.dotmatch')]
            forensics.command(argv, out / f'{label}.dotmatch_command')
            native = json.loads((out / f'{label}.dotmatch/summary.json').read_text())
            pilot.require(native['completion'] == 'complete' and native['read_count'] == int(summary['samples'][label]['read_count']), 'native completion/denominator mismatch')
            pilot.require(native['inputs']['reads']['sha256'] == summary['samples'][label]['local_sha256'], 'native input digest mismatch')
            pilot.require(native['inputs']['targets']['sha256'] == provenance['normalized_sha256'], 'native library digest mismatch')
            for filename, info in native['artifacts'].items():
                pilot.require(pilot.digest(out / f'{label}.dotmatch' / filename) == info['sha256'], 'native artifact hash mismatch')
        decoder = Decoder(library)
        for label, path in zip(labels, paths):
            result = independent_sample(path, label, decoder, summary['samples'][label], out)
            summary['results'][label] = result
            native = json.loads((out / f'{label}.dotmatch/summary.json').read_text())
            for policy in pilot.POLICIES:
                compare_matrices(out / f'{label}.oracle.{policy}.tsv', out / f'{label}.dotmatch/{policy}.counts.tsv', label)
                pilot.require(all(native['outcomes'][policy].get(s, 0) == result['states'][policy][s] for s in pilot.STATES), 'full native state disagreement')
            pilot.require(native['changed_reads'] == result['changed_records'], 'full changed-read disagreement')
            compare_matrices(out / f'{label}.oracle.events.tsv', out / 'guide_counter.counts.txt', label)
            result['independent_all_guide_count_agreement'] = True
            pilot.save_json(out / 'run.json', summary)
        summary['status'] = 'complete'
    except Exception as exc:
        summary['status'] = 'failed'
        summary['error'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        pilot.save_json(out / 'run.json', summary)
        files = {str(p.relative_to(out)): {'sha256': pilot.digest(p), 'bytes': p.stat().st_size}
            for p in sorted(out.rglob('*')) if p.is_file() and 'work' not in p.relative_to(out).parts and p.name != 'MANIFEST.json'}
        pilot.save_json(out / 'MANIFEST.json', {'status': summary['status'], 'files': files})


class Tests(unittest.TestCase):
    def test_fast_decoder_exhaustive(self):
        rng = random.Random(20260906)
        strings = [''.join(s) for s in itertools.product('ACGT', repeat=4)]
        for _ in range(5):
            library = [(str(i), s, str(i)) for i, s in enumerate(rng.sample(strings, 20))]
            decoder = Decoder(library)
            for bases in itertools.product('ACGTN', repeat=4):
                window = ''.join(bases)
                decoder.check(window)
                observed = tuple(('unique', c) if c >= 0 else (STATE[c], -1) for c in decoder.calls(window))
                self.assertEqual(observed, decoder.ref.brute(window))
            decoder.check(None)

    def test_duplicate_reference_is_rejected(self):
        with self.assertRaises(ValueError):
            Decoder([('a', 'AAAA', 'g'), ('b', 'AAAA', 'h')])

    def test_literal_noncanonical(self):
        decoder = Decoder([('a', 'AAAA', 'g'), ('b', 'CCCC', 'h')])
        for window in ('AAAN', 'AAAR', 'AANN', 'TTTT', None):
            decoder.check(window)
        self.assertEqual(decoder.calls('AAAR'), (NONE, 0, 0))
        self.assertEqual(decoder.calls('AANN'), (NONE, NONE, NONE))

    def test_named_sample_axes_and_annotation(self):
        import tempfile
        with tempfile.TemporaryDirectory() as directory:
            a, b = Path(directory) / 'a.tsv', Path(directory) / 'b.tsv'
            pilot.table(a, ['sgRNA', 'Gene', 'x', 'y'], [('a', 'G', 1, 2)])
            pilot.table(b, ['sgRNA', 'Gene', 'y', 'x'], [('a', 'G', 2, 1)])
            compare_matrices(a, b, 'x')
            pilot.table(b, ['sgRNA', 'Gene', 'y', 'x'], [('a', 'H', 2, 1)])
            with self.assertRaises(ValueError):
                compare_matrices(a, b, 'x')


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        unittest.main(argv=[sys.argv[0]], verbosity=2)
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--out', required=True)
        run(Path(parser.parse_args().out).resolve())
