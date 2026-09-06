#!/usr/bin/env python3
"""AR-001 source diagnostics and historical Yusa count reconstruction.

The reconstruction counts matching offset events, then separately counts reads.
It deliberately does not interpret an event count as a molecule count.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import itertools
import json
import shutil
import subprocess
import sys
import time
import zipfile
from collections import Counter
from pathlib import Path
import pilot


def inspect(out):
    out.mkdir(parents=True, exist_ok=False)
    report = {}
    for name, config in pilot.SOURCES.items():
        data = pilot.fetch_bytes(config['library_url'])
        with zipfile.ZipFile(io.BytesIO(data)) as zf:
            member = next(x for x in zf.infolist() if x.filename.split('/')[-1] == config['member'])
            pilot.require(member.file_size < 32 * 1024 * 1024, 'oversized source')
            raw = zf.read(member)
        text = raw.decode('utf-8-sig')
        delimiter = '\t' if '\t' in text.splitlines()[0] else ','
        rows = list(csv.reader(io.StringIO(text), delimiter=delimiter))
        body = rows[1:]
        errors = [(i + 2, row) for i, row in enumerate(body) if len(row) != 3 or len(row[1]) != config['length'] or set(row[1].upper()) - set('ACGT')]
        result = {'url': config['library_url'], 'source_sha256': hashlib.sha256(raw).hexdigest(),
            'archive_sha256': hashlib.sha256(data).hexdigest(), 'header': rows[0], 'data_rows': len(body),
            'raw_length_distribution': dict(Counter(len(row[1]) for row in body if len(row) >= 2)),
            'rows_failing_expected_length_or_alphabet': len(errors), 'first_failing_rows': errors[:50]}
        report[name] = result
        pilot.save_json(out / 'reference-inspection.json', report)
        print(json.dumps({'reference_inspection': name, **result}), flush=True)
    return report


def command(argv, prefix):
    started = time.perf_counter()
    run = subprocess.run(argv, text=True, capture_output=True, timeout=600)
    Path(str(prefix) + '.stdout.txt').write_text(run.stdout)
    Path(str(prefix) + '.stderr.txt').write_text(run.stderr)
    pilot.save_json(Path(str(prefix) + '.command.json'), {'argv': argv, 'returncode': run.returncode, 'seconds': time.perf_counter() - started})
    pilot.require(run.returncode == 0, f'command failed: {argv[0]}: {run.stderr[-3000:]}')


def matrix(path):
    with Path(path).open() as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        pilot.require(len(header) >= 3 and len(set(header)) == len(header), 'invalid matrix header')
        result = {}
        for row in reader:
            pilot.require(len(row) == len(header) and row[0] not in result, 'duplicate/malformed matrix row')
            pilot.require(all(x.isdecimal() for x in row[2:]), 'non-integer matrix count')
            result[row[0]] = (row[1], dict(zip(header[2:], map(int, row[2:]))))
    return header[2:], result


def matching_event_lookup(library):
    """Enumerate ACGT one-substitution codewords; colliding mutants abstain.

    Exact sequences replace their mutant entries. For this historical lane,
    duplicate reference sequences are rejected so no library-order tie is hidden.
    """
    lookup = {}
    pilot.require(len({row[1] for row in library}) == len(library), 'duplicate sequence in historical lane')
    for index, (_ident, sequence, _gene) in enumerate(library):
        for pos, old in enumerate(sequence):
            for base in 'ACGT':
                if base == old:
                    continue
                mutant = sequence[:pos] + base + sequence[pos + 1:]
                previous = lookup.get(mutant)
                if previous is None:
                    lookup[mutant] = index
                elif previous != index:
                    lookup[mutant] = -1
    for index, (_ident, sequence, _gene) in enumerate(library):
        lookup[sequence] = index
    return lookup


def reconstruct(recs, library, lookup, sample, out):
    length = len(library[0][1])
    all_offsets = Counter()
    for _h, seq, _p, _q in recs:
        for offset in range(len(seq) - length + 1):
            if lookup.get(seq[offset:offset + length], -1) >= 0:
                all_offsets[offset] += 1
    total = sum(all_offsets.values())
    pilot.require(total > 0, 'no matching offset events')
    # Integer comparison is exactly n/total >= 0.0025 = 1/400.
    selected = sorted(offset for offset, n in all_offsets.items() if n * 400 >= total)
    counts, once_consensus, multiplicity = Counter(), Counter(), Counter()
    read_any = same = different = 0
    multi = []
    for ordinal, (_h, seq, _p, _q) in enumerate(recs, 1):
        events = [(offset, lookup.get(seq[offset:offset + length], -1)) for offset in selected if offset + length <= len(seq)]
        events = [(offset, guide) for offset, guide in events if guide >= 0]
        multiplicity[len(events)] += 1
        if events:
            read_any += 1
            ids = {i for _, i in events}
            if len(ids) == 1:
                once_consensus[next(iter(ids))] += 1
            for _, guide in events:
                counts[guide] += 1
            if len(events) > 1:
                same += len(ids) == 1
                different += len(ids) > 1
                multi.append([ordinal, ';'.join(f'{offset}:{library[guide][0]}' for offset, guide in events), len(ids)])
    pilot.table(out / f'{sample}.multi_offset_records.tsv', ['record_ordinal', 'offset_and_guide', 'distinct_guides'], multi)
    pilot.table(out / f'{sample}.offsets.tsv', ['offset', 'matching_events', 'selected'], [(k, v, int(k in selected)) for k, v in sorted(all_offsets.items())])
    pilot.require(sum(counts.values()) == sum(k * v for k, v in multiplicity.items()), 'event conservation failed')
    result = {'records': len(recs), 'selected_offsets': selected, 'all_offset_discovery_events': total,
        'matching_events': sum(counts.values()), 'reads_with_at_least_one_event': read_any,
        'excess_events_above_matched_reads': sum(counts.values()) - read_any,
        'events_above_total_records': max(0, sum(counts.values()) - len(recs)),
        'multi_offset_same_guide_reads': same, 'multi_offset_different_guide_reads': different,
        'event_multiplicity': dict(multiplicity), 'single_guide_consensus_reads': sum(once_consensus.values())}
    print(json.dumps({'reconstruction': sample, **result}), flush=True)
    return counts, result


def yusa(out):
    out.mkdir(parents=True, exist_ok=False)
    work = out / 'work'
    work.mkdir()
    summary = {'status': 'running', 'scope': 'historical first 100000 records per Yusa sample, includes discovery reads',
               'harness_sha256': pilot.digest(__file__), 'protocol_sha256': pilot.digest(pilot.HERE / 'PROTOCOL.md'),
               'samples': {}, 'results': {}, 'claims': 'Counting events and read records are different units; not biological accuracy.'}
    try:
        library, lib_path, provenance = pilot.fetch_library('yusa', work)
        summary['library'] = provenance
        paths, recsets = [], []
        for label, accessions in pilot.SOURCES['yusa']['samples'].items():
            recs, metas = pilot.fetch_prefix(accessions, 100000)
            path = work / f'{label}.fastq.gz'
            pilot.write_fastq(path, recs)
            paths.append(path)
            recsets.append(recs)
            summary['samples'][label] = {'sources': metas, 'records': len(recs), 'sha256': pilot.digest(path)}
        labels = list(pilot.SOURCES['yusa']['samples'])
        summary['guide_counter_version'] = subprocess.check_output(['guide-counter', '--version'], text=True).strip()
        summary['guide_counter_binary_sha256'] = pilot.digest(shutil.which('guide-counter'))
        summary['dotmatch_version'] = subprocess.check_output(['dotmatch', '--version'], text=True).strip()
        crate_paths = list((Path.home() / '.cargo/registry/src').glob('*/guide-counter-0.1.3'))
        summary['crate_source_files'] = {}
        for crate in crate_paths:
            for rel in ('Cargo.lock', '.cargo_vcs_info.json', 'src/commands/count.rs', 'src/guide.rs'):
                path = crate / rel
                if path.exists():
                    summary['crate_source_files'][rel] = {'sha256': pilot.digest(path)}
                    if rel == '.cargo_vcs_info.json':
                        summary['crate_source_files'][rel]['content'] = json.loads(path.read_text())
        for exact in (False, True):
            label = 'guide_counter_exact' if exact else 'guide_counter_hamming'
            argv = ['guide-counter', 'count', '--input', *map(str, paths), '--samples', *labels,
                    '--library', str(lib_path), '--output', str(out / label), '--offset-sample-size', '100000', '--offset-min-fraction', '0.0025']
            if exact:
                argv.append('--exact-match')
            command(argv, out / label)
        common = ['dotmatch', 'count', '--targets', str(lib_path), '--sample-label', ','.join(labels),
                  '--target-start', '23', '--target-length', '19', '--metric', 'hamming', '--format', 'mageck']
        for path in paths:
            common += ['--reads', str(path)]
        for policy in ('exact', 'hamming'):
            argv = common + ['--k', '0' if policy == 'exact' else '1', '--out', str(out / f'dotmatch_{policy}.counts.tsv'),
                             '--summary', str(out / f'dotmatch_{policy}.summary.json')]
            if policy == 'hamming':
                argv += ['--ambiguity-policy', 'best', '--auto-offset', '5', '--auto-offset-sample', '100000', '--offset-min-fraction', '0.0025']
            command(argv, out / f'dotmatch_{policy}')
        gc_labels, gc = matrix(out / 'guide_counter_hamming.counts.txt')
        dm_labels, dm = matrix(out / 'dotmatch_hamming.counts.tsv')
        pilot.require(set(gc_labels) == set(dm_labels) == set(labels), 'sample axes differ')
        pilot.require(gc.keys() == dm.keys() == {row[0] for row in library}, 'guide axes differ')
        lookup = matching_event_lookup(library)
        ref = pilot.Reference(library)
        # Check lookup semantics on deterministic real windows against a separate query-neighbour oracle.
        checked = 0
        for recs in recsets:
            for _h, seq, _p, _q in recs[:1000]:
                window = seq[23:42]
                if len(window) == 19 and set(window) <= set('ACGT'):
                    call = ref.calls(window)[2]
                    pilot.require(lookup.get(window, -1) == (call[1] if call[0] == 'unique' else -1), 'lookup/oracle discrepancy')
                    checked += 1
        summary['lookup_oracle_checks'] = checked
        for label, recs in zip(labels, recsets):
            predicted, result = reconstruct(recs, library, lookup, label, out)
            mismatches = []
            for i, (ident, _seq, gene) in enumerate(library):
                pilot.require(gc[ident][0] == dm[ident][0] == gene, 'gene annotation differs')
                if gc[ident][1][label] != predicted[i]:
                    mismatches.append((ident, gc[ident][1][label], predicted[i]))
            result['reconstruction_count_disagreements'] = len(mismatches)
            result['dotmatch_total'] = sum(row[1][label] for row in dm.values())
            result['guide_counter_total'] = sum(row[1][label] for row in gc.values())
            result['dotmatch_vs_guide_counter_differing_guides'] = sum(dm[x][1][label] != gc[x][1][label] for x in dm)
            summary['results'][label] = result
            pilot.table(out / f'{label}.reconstruction_disagreements.tsv', ['guide', 'guide_counter', 'reconstructed'], mismatches)
            pilot.require(not mismatches, 'independent reconstruction does not match guide-counter')
        dm_total = {x: sum(dm[x][1].values()) for x in dm}
        gc_total = {x: sum(gc[x][1].values()) for x in gc}
        observed = {'dotmatch': sum(dm_total.values()), 'guide_counter': sum(gc_total.values()),
                    'differing_guides': sum(dm_total[x] != gc_total[x] for x in dm)}
        summary['historical_aggregate'] = observed
        summary['historical_report_exactly_reproduced'] = observed == {'dotmatch': 184167, 'guide_counter': 208700, 'differing_guides': 13537}
        summary['status'] = 'complete'
    except Exception as exc:
        summary['status'] = 'failed'
        summary['error'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        pilot.save_json(out / 'reconstruction.json', summary)
        print(json.dumps(summary), flush=True)
        files = {str(p.relative_to(out)): {'bytes': p.stat().st_size, 'sha256': pilot.digest(p)}
                 for p in sorted(out.rglob('*')) if p.is_file() and 'work' not in p.relative_to(out).parts and p.name != 'MANIFEST.json'}
        pilot.save_json(out / 'MANIFEST.json', {'status': summary['status'], 'files': files})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['inspect', 'yusa'])
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    (inspect if args.action == 'inspect' else yusa)(Path(args.out).resolve())
