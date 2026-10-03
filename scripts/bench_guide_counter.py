#!/usr/bin/env python3
"""Paired complete-command comparison with unmodified fulcrumgenomics/guide-counter."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SUFFIXES = ('counts.txt', 'extended-counts.txt', 'stats.txt')
STATS_FIELDS = {'file', 'label', 'total_guides', 'total_reads', 'mapped_reads',
                'frac_mapped', 'mean_reads_per_guide', 'mean_reads_essential',
                'mean_reads_nonessential', 'mean_reads_control', 'mean_reads_other', 'zero_read_guides'}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def library_rows(path):
    with Path(path).open() as stream:
        header = stream.readline()
        delimiter = '\t' if header.count('\t') >= 2 else ','
        stream.seek(0)
        rows = list(csv.reader(stream, delimiter=delimiter))[1:]
    if not rows or len({len(row[1]) for row in rows}) != 1 or len({row[1] for row in rows}) != len(rows):
        raise ValueError('Benchmark requires distinct guide sequences of one length')
    return rows


def write_reads(path, targets, n_reads, seed):
    """70% exact, 20% one substitution, 5% two substitutions, 5% random decoys.

    Guides start at offsets 23/24/25 in 50-base reads; every byte is uppercase ACGT.
    The generation label describes construction, not guaranteed source identity.
    """
    rng = random.Random(seed)
    length = len(targets[0])
    with path.open('wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0, compresslevel=6) as stream:
        for i in range(n_reads):
            seq = bytearray(rng.choices(b'ACGT', k=50))
            guide = bytearray(targets[rng.randrange(len(targets))])
            category = i % 20
            errors = 0 if category < 14 else (1 if category < 18 else 2)
            if category == 19:
                guide = bytearray(rng.choices(b'ACGT', k=length))
            else:
                for position in rng.sample(range(length), errors):
                    guide[position] = rng.choice([base for base in b'ACGT' if base != guide[position]])
            offset = 23 + i % 3
            seq[offset:offset + length] = guide
            stream.write(b'@controlled_' + str(i).encode() + b'\n' + seq + b'\n+\n' + b'I' * len(seq) + b'\n')


def prepare(folder, library, sizes, seed, experimental):
    targets = [row[1].encode('ascii') for row in library_rows(library)]
    cases = []
    for n_reads, samples in [(size, 1) for size in sizes] + [(max(sizes), 4)]:
        if n_reads % samples:
            raise ValueError('read sizes must be divisible by four')
        name = f'controlled_yusa_{n_reads}_{samples}sample'
        paths = []
        for sample in range(samples):
            path = folder / f'{name}-{sample}.fastq.gz'
            write_reads(path, targets, n_reads // samples, seed + n_reads + sample + 100 * samples)
            paths.append(path)
        cases.append({'workload': name, 'scope': 'controlled_simulation', 'library': str(library.resolve()),
                      'reads': [str(p.resolve()) for p in paths], 'labels': [f'sample{i}' for i in range(samples)],
                      'n_reads': n_reads, 'n_targets': len(targets), 'guide_length': len(targets[0]),
                      'seed': seed + n_reads + 100 * samples})
    if experimental:
        sizes = []
        for path in experimental:
            with gzip.open(path, 'rb') as stream:
                sizes.append(sum(1 for _ in stream) // 4)
        cases.append({'workload': 'yusa_experimental_cached_prefix', 'scope': 'experimental_prefix',
                      'library': str(library.resolve()), 'reads': [str(p.resolve()) for p in experimental],
                      'labels': [p.name.split('.')[0] for p in experimental], 'n_reads': sum(sizes),
                      'n_targets': len(targets), 'guide_length': len(targets[0]), 'records_by_sample': sizes,
                      'source_urls': [f'https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR376/{p.name.split(".")[0]}/{p.name}'
                                      for p in experimental]})
    protocol = {'schema_version': 1, 'seed': seed, 'cases': cases,
                'controlled_read_recipe': write_reads.__doc__,
                'files': {str(path.resolve()): digest(path) for path in [library, *[Path(p) for c in cases for p in c['reads']]]}}
    (folder / 'protocol.json').write_text(json.dumps(protocol, indent=2) + '\n')
    return protocol


def timed(command, log, launcher):
    """A small fresh launcher prevents Python heap size from inflating child RSS."""
    measurement = log.with_suffix('.measurement.json')
    with log.open('wb') as stream:
        process = subprocess.run([str(launcher.resolve()), str(measurement.resolve()), *command],
                                 stdout=stream, stderr=stream, env={**os.environ, 'RUST_LOG': 'info'})
    if process.returncode:
        raise RuntimeError(f'Command failed ({process.returncode}); see {log}')
    values = json.loads(measurement.read_text())
    return values['seconds'], values['cpu_seconds'], values['peak_rss_kib'], values['exit_code']


def outputs(prefix):
    count_path = Path(str(prefix) + '.counts.txt')
    with count_path.open() as stream:
        reader = csv.reader(stream, delimiter='\t')
        header = next(reader)
        labels = header[2:]
        rows = list(reader)
    if header[:2] != ['guide', 'gene'] or not labels or any(len(row) != len(header) for row in rows):
        raise ValueError('Malformed count output axes or rows')
    guides = {row[0]: row[1] for row in rows}
    if len(guides) != len(rows) or len(set(labels)) != len(labels):
        raise ValueError('Duplicate output axes')
    matrix = {row[0]: tuple(int(value) for value in row[2:]) for row in rows}
    if any(len(values) != len(labels) or min(values) < 0 for values in matrix.values()):
        raise ValueError('Malformed count matrix')
    canonical = json.dumps({'labels': labels, 'genes': guides, 'matrix': matrix}, sort_keys=True, separators=(',', ':'))
    with Path(str(prefix) + '.extended-counts.txt').open() as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError('Duplicate or missing extended output columns')
        if set(reader.fieldnames) != {'guide', 'gene', 'guide_type', *labels}:
            raise ValueError('Extended output columns differ')
        extended_rows = list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in extended_rows):
            raise ValueError('Malformed extended output row')
        extended = {row['guide']: row for row in extended_rows}
        if len(extended) != len(extended_rows):
            raise ValueError('Duplicate extended output guides')
    if set(extended) != set(matrix):
        raise ValueError('Extended output axes differ')
    for guide, values in matrix.items():
        if tuple(int(extended[guide][label]) for label in labels) != values or extended[guide]['gene'] != guides[guide]:
            raise ValueError('Extended count matrix differs')
    with Path(str(prefix) + '.stats.txt').open() as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError('Duplicate or missing stats columns')
        if set(reader.fieldnames) != STATS_FIELDS:
            raise ValueError('Statistics columns differ')
        stats = list(reader)
    for row in stats:
        if None in row or any(value is None for value in row.values()):
            raise ValueError('Malformed statistics row')
        if any(not math.isfinite(float(value)) for key, value in row.items() if key not in ('file', 'label')):
            raise ValueError('Nonfinite output statistic')
    if [row['label'] for row in stats] != labels:
        raise ValueError('Stats sample axis differs')
    totals = [sum(values[j] for values in matrix.values()) for j in range(len(labels))]
    for j, row in enumerate(stats):
        if int(row['mapped_reads']) != totals[j] or int(row['total_guides']) != len(matrix):
            raise ValueError('Stats/count disagreement')
    return {'labels': labels, 'genes': guides, 'matrix': matrix, 'stats': stats, 'extended': extended,
            'counts_sha256': hashlib.sha256(canonical.encode()).hexdigest(), 'mapped_by_sample': totals}


def compare(a, b):
    if a['labels'] != b['labels'] or a['genes'] != b['genes']:
        raise ValueError('Output library/sample identity differs')
    if len(a['stats']) != len(a['labels']) or len(b['stats']) != len(b['labels']):
        raise ValueError('Statistics sample coverage differs')
    count_errors = sum(x != y for guide in a['matrix'] for x, y in zip(a['matrix'][guide], b['matrix'][guide]))
    extended_errors = sum(a['extended'][guide] != b['extended'][guide] for guide in a['extended'])
    stats_errors = 0
    for x, y in zip(a['stats'], b['stats']):
        if set(x) != set(y):
            raise ValueError('Statistics fields differ')
        for key in x:
            if key in ('file', 'label'):
                stats_errors += x[key] != y[key]
            else:
                if not math.isfinite(float(x[key])) or not math.isfinite(float(y[key])):
                    raise ValueError('Nonfinite output statistic')
                stats_errors += abs(float(x[key]) - float(y[key])) > 1e-6
    return count_errors, extended_errors, stats_errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--guide-counter', required=True, type=Path)
    parser.add_argument('--guide-counter-source', required=True, type=Path)
    parser.add_argument('--dotmatch', default=ROOT / 'dotmatch', type=Path)
    parser.add_argument('--rustc', type=Path, help='compiler used for the competitor build, for provenance')
    parser.add_argument('--library', default=ROOT / 'benchmarks/real/data/yusa_library.csv', type=Path)
    parser.add_argument('--sizes', default='100000,1000000')
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--seed', type=int, default=20261002)
    parser.add_argument('--experimental-fastq', action='append', type=Path, default=[])
    parser.add_argument('--work-dir', default=ROOT / 'benchmarks/work/guide_counter', type=Path)
    parser.add_argument('--out', default=ROOT / 'benchmarks/raw/guide_counter_comparison.csv', type=Path)
    parser.add_argument('--reuse-inputs', action='store_true')
    args = parser.parse_args()
    sizes = [int(value) for value in args.sizes.split(',')]
    if min(sizes) < 1 or args.repeats < 2:
        parser.error('positive sizes and at least two repeats required')
    if args.dotmatch.resolve() != (ROOT / 'dotmatch').resolve():
        parser.error('this source-audited comparison requires the in-tree DotMatch binary')
    subprocess.run(['make', 'dotmatch', 'CC=cc'], cwd=ROOT, check=True)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    launcher = args.work_dir / 'measure-command'
    subprocess.run(['cc', '-O2', '-std=c11', '-Wall', '-Wextra', str(ROOT / 'benchmarks/measure_command.c'),
                    '-o', str(launcher)], check=True)
    if args.reuse_inputs:
        protocol = json.loads((args.work_dir / 'protocol.json').read_text())
        for path, sha in protocol['files'].items():
            if digest(path) != sha:
                raise ValueError(f'Input changed: {path}')
    else:
        protocol = prepare(args.work_dir, args.library, sizes, args.seed, args.experimental_fastq)
    source = args.guide_counter_source.resolve()
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=source, text=True).strip()
    if subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=source, text=True).strip():
        raise ValueError('Competitor source must be unmodified')
    metadata = {'schema_version': 1, 'command': [sys.executable, *sys.argv], 'platform': platform.platform(),
                'python': platform.python_version(), 'cpu_count': os.cpu_count(),
                'dotmatch_version': subprocess.check_output([str(args.dotmatch.resolve()), '--version'], text=True).strip(),
                'cc_version': subprocess.check_output(['cc', '--version'], text=True).splitlines()[0],
                'rustc_version': subprocess.check_output([str(args.rustc.resolve()), '--version'], text=True).strip() if args.rustc else 'not supplied',
                'affinity': sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None,
                'guide_counter_commit': commit, 'guide_counter_version': '0.1.3',
                'guide_counter_build': 'cargo build --release --locked; upstream fat LTO and one codegen unit',
                'dotmatch_build': 'make dotmatch CC=cc; default -O3 and -mavx2 on x86-64',
                'thread_scope': 'both commands single-threaded; fresh process per timing',
                'timing_scope': 'includes index/library construction, offset detection, gzip parsing, counting and all three output files; excludes generation and checks',
                'offset_sample_size': 100000, 'offset_min_fraction': .0025,
                'counting_unit': 'one unique best-distance guide assignment per selected offset; multiple assignments per read permitted',
                'protocol': protocol,
                'files': {str(path.resolve()): digest(path) for path in [args.dotmatch, args.guide_counter,
                           ROOT / 'src/qda.c', ROOT / 'src/qdalign.c', Path(__file__), source / 'Cargo.lock', source / 'Cargo.toml']}}
    metadata['measurement'] = {'method': 'fresh native launcher with monotonic clock and per-child wait4 rusage',
                               'source_sha256': digest(ROOT / 'benchmarks/measure_command.c'),
                               'binary_sha256': digest(launcher)}
    cpuinfo = Path('/proc/cpuinfo')
    metadata['cpu_model'] = next((line.split(':', 1)[1].strip() for line in cpuinfo.read_text().splitlines()
                                  if line.startswith('model name')), 'unknown') if cpuinfo.exists() else platform.processor()
    rows = []
    for case in protocol['cases']:
        print(f'{case["workload"]}: {case["n_reads"]} reads / {case["n_targets"]} guides', flush=True)
        for k in (0, 1):
            expected = None
            for repeat in range(args.repeats):
                order = ('guide_counter', 'dotmatch') if repeat % 2 == 0 else ('dotmatch', 'guide_counter')
                paired = {}
                for order_index, tool in enumerate(order):
                    prefix = args.work_dir / f'{case["workload"]}-k{k}-r{repeat}-{tool}'
                    command = [str(args.guide_counter.resolve()), 'count'] if tool == 'guide_counter' else [str(args.dotmatch.resolve()), 'guide-counter', 'count']
                    command += ['--input', *case['reads'], '--samples', *case['labels'], '--library', case['library'],
                                '--output', str(prefix.resolve()), '--offset-sample-size', '100000', '--offset-min-fraction', '.0025']
                    if k == 0:
                        command.append('--exact-match')
                    elapsed, cpu, rss, code = timed(command, prefix.with_suffix('.log'), launcher)
                    result = outputs(prefix)
                    if len(result['matrix']) != case['n_targets'] or sum(int(s['total_reads']) for s in result['stats']) != case['n_reads']:
                        raise RuntimeError('Output coverage differs from protocol')
                    paired[tool] = result
                    if expected is None:
                        expected = result['counts_sha256']
                    if result['counts_sha256'] != expected:
                        raise RuntimeError(f'Count disagreement: {case["workload"]}/k{k}/{tool}; see {prefix}')
                    rows.append({'workload': case['workload'], 'scope': case['scope'], 'k': k, 'repeat': repeat,
                                 'order': order_index, 'tool': tool, 'n_reads': case['n_reads'], 'n_targets': case['n_targets'],
                                 'samples': len(case['reads']), 'seconds': elapsed, 'cpu_seconds': cpu,
                                 'reads_per_second': case['n_reads'] / elapsed, 'peak_rss_kib': rss,
                                 'counts_sha256': result['counts_sha256'], 'count_mismatches': 0,
                                 'extended_mismatches': 0, 'stats_mismatches': 0,
                                 'mapped_by_sample': json.dumps(result['mapped_by_sample']), 'exit_code': code,
                                 'output_sha256': json.dumps({suffix: digest(Path(str(prefix) + '.' + suffix)) for suffix in SUFFIXES}),
                                 'command': json.dumps(command)})
                errors = compare(paired['guide_counter'], paired['dotmatch'])
                if any(errors):
                    raise RuntimeError(f'Output disagreement: {case["workload"]}/k{k}: {errors}')
                print(f'  k={k} repeat={repeat}: complete outputs agree', flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata['csv_sha256'] = digest(args.out)
    args.out.with_suffix('.json').write_text(json.dumps(metadata, indent=2) + '\n')


if __name__ == '__main__':
    main()
