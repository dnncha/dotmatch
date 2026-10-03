#!/usr/bin/env python3
"""Audit and regenerate the paired guide-counter comparison and figures."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'benchmarks/raw'
REPORT = ROOT / 'docs/benchmarks/guide_counter/README.md'
FIGURES = ROOT / 'benchmarks/figures'
STEM = 'guide_counter_comparison'
TOOLS = ('dotmatch', 'guide_counter')
CALLOUTS = (
    'README.md', 'app/guide-counter-benchmark.tsx', 'docs/index.md',
    'docs/crisprworks.md', 'docs/tutorials/crispr-count-first-run.md',
    'docs/benchmarks/README.md', 'docs/usability-comparison.md',
    'docs/scientific-claims.md', 'scripts/generate_benchmark_report.py',
    'llms.txt', 'docs/llms.txt', 'public/llms.txt',
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_validation_record(recorded, rows):
    expected_figures = {f'guide_counter_{kind}.{extension}'
                        for kind in ('throughput', 'memory') for extension in ('svg', 'png', 'pdf')}
    if (recorded.get('schema_version') != 1 or recorded.get('commands') != len(rows)
            or any(recorded.get(field) != 0 for field in
                   ('count_mismatches', 'extended_mismatches', 'stats_mismatches'))):
        raise AssertionError('invalid validation summary')
    if not isinstance(recorded.get('figures'), dict) or set(recorded['figures']) != expected_figures:
        raise AssertionError('incomplete figure validation')
    for name, sha in recorded['figures'].items():
        if not isinstance(sha, str) or len(sha) != 64 or any(c not in '0123456789abcdef' for c in sha):
            raise AssertionError(f'invalid figure digest: {name}')


def check_callout_headlines(stats):
    mismatch = [s for s in stats if s['k'] == 1 and s['case']['scope'] == 'controlled_simulation']
    intervals = [f'{min(s[field] for s in mismatch):.1f}–{max(s[field] for s in mismatch):.1f}×'
                 for field in ('speedup', 'memory_ratio')]
    for relative in CALLOUTS:
        text = (ROOT / relative).read_text()
        if any(interval not in text for interval in intervals):
            raise AssertionError(f'benchmark headline differs from raw evidence: {relative}')


def option(command, key):
    i = command.index(key) + 1
    values = []
    while i < len(command) and not command[i].startswith('--'):
        values.append(command[i])
        i += 1
    return values


def load(check_source=False):
    path = RAW / f'{STEM}.csv'
    metadata = json.loads(path.with_suffix('.json').read_text())
    if not (metadata['csv_sha256'] == digest(path)):
        raise AssertionError('CSV hash mismatch')
    if not (len(metadata['guide_counter_commit']) == 40):
        raise AssertionError('missing competitor source pin')
    if not (metadata['guide_counter_version'] == '0.1.3'):
        raise AssertionError('unexpected comparator version')
    if not (metadata['offset_sample_size'] == 100000 and metadata['offset_min_fraction'] == .0025):
        raise AssertionError('evidence validation failed')
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    cases = {case['workload']: case for case in metadata['protocol']['cases']}
    if not (len(cases) == len(metadata['protocol']['cases'])):
        raise AssertionError('duplicate workloads')
    if not ({case['n_reads'] for case in cases.values() if case['scope'] == 'controlled_simulation'} == {100000, 1000000}):
        raise AssertionError('missing read scales')
    if not (any(len(case['reads']) == 4 for case in cases.values())):
        raise AssertionError('missing multi-sample workload')
    groups = {}
    for row in rows:
        if not (row['workload'] in cases and row['tool'] in TOOLS):
            raise AssertionError('unexpected workload/tool')
        case = cases[row['workload']]
        if not (row['scope'] == case['scope']):
            raise AssertionError('input scope differs')
        if not (int(row['n_reads']) == case['n_reads'] and int(row['n_targets']) == case['n_targets'] == 87437):
            raise AssertionError('evidence validation failed')
        if not (int(row['samples']) == len(case['reads'])):
            raise AssertionError('sample coverage differs')
        if not (int(row['exit_code']) == 0):
            raise AssertionError('command failed')
        if not (all(int(row[field]) == 0 for field in ('count_mismatches', 'extended_mismatches', 'stats_mismatches'))):
            raise AssertionError('output disagreement')
        if not (all(math.isfinite(float(row[field])) and float(row[field]) > 0 for field in
                   ('seconds', 'cpu_seconds', 'peak_rss_kib', 'reads_per_second'))):
            raise AssertionError('invalid measurements')
        if not (math.isclose(float(row['reads_per_second']), int(row['n_reads']) / float(row['seconds']), rel_tol=1e-12)):
            raise AssertionError('evidence validation failed')
        if not (len(row['counts_sha256']) == 64):
            raise AssertionError('missing count digest')
        if not (set(json.loads(row['output_sha256'])) == {'counts.txt', 'extended-counts.txt', 'stats.txt'}):
            raise AssertionError('evidence validation failed')
        if not (len(json.loads(row['mapped_by_sample'])) == len(case['reads'])):
            raise AssertionError('incomplete count totals')
        command = json.loads(row['command'])
        if not (command[1:3] == ['guide-counter', 'count'] if row['tool'] == 'dotmatch' else command[1] == 'count'):
            raise AssertionError('evidence validation failed')
        if not (option(command, '--input') == case['reads'] and option(command, '--samples') == case['labels']):
            raise AssertionError('evidence validation failed')
        if not (option(command, '--library') == [case['library']]):
            raise AssertionError('different guide library')
        if not (option(command, '--offset-sample-size') == ['100000']):
            raise AssertionError('evidence validation failed')
        if not (option(command, '--offset-min-fraction') == ['.0025']):
            raise AssertionError('evidence validation failed')
        if not (('--exact-match' in command) == (int(row['k']) == 0)):
            raise AssertionError('metric differs')
        key = row['workload'], int(row['k'])
        slot = int(row['repeat']), row['tool']
        if not (slot not in groups.setdefault(key, {})):
            raise AssertionError('duplicate timing')
        groups[key][slot] = row
    expected = {(case, k) for case in cases for k in (0, 1)}
    if not (set(groups) == expected):
        raise AssertionError('missing metric/workload')
    stats = []
    for (name, k), paired in groups.items():
        repeats = sorted({repeat for repeat, _ in paired})
        if not (repeats == list(range(5))):
            raise AssertionError('five repeats required')
        if not (set(paired) == {(r, tool) for r in repeats for tool in TOOLS}):
            raise AssertionError('missing paired run')
        if not (len({row['counts_sha256'] for row in paired.values()}) == 1):
            raise AssertionError('count digests differ')
        if not (len({row['mapped_by_sample'] for row in paired.values()}) == 1):
            raise AssertionError('count totals differ')
        speedups, memory_ratios = [], []
        for repeat in repeats:
            dot, guide = (paired[repeat, tool] for tool in TOOLS)
            if not (int(guide['order']) == repeat % 2 and int(dot['order']) == 1 - repeat % 2):
                raise AssertionError('nonalternating order')
            speedups.append(float(guide['seconds']) / float(dot['seconds']))
            memory_ratios.append(float(guide['peak_rss_kib']) / float(dot['peak_rss_kib']))
        result = {'workload': name, 'k': k, 'case': cases[name], 'speedup': statistics.median(speedups),
                  'speedup_min': min(speedups), 'speedup_max': max(speedups),
                  'memory_ratio': statistics.median(memory_ratios)}
        for tool in TOOLS:
            entries = [paired[r, tool] for r in repeats]
            result[tool] = {field: statistics.median(float(row[field]) for row in entries) for field in
                            ('seconds', 'reads_per_second', 'peak_rss_kib')}
            for field in ('reads_per_second', 'peak_rss_kib'):
                result[tool][field + '_min'] = min(float(row[field]) for row in entries)
                result[tool][field + '_max'] = max(float(row[field]) for row in entries)
        stats.append(result)
    if check_source:
        for relative in ('src/qda.c', 'src/qdalign.c', 'scripts/bench_guide_counter.py'):
            recorded = [sha for name, sha in metadata['files'].items() if name.endswith('/' + relative)]
            if not (recorded == [digest(ROOT / relative)]):
                raise AssertionError(f'source changed: {relative}')
        if not (metadata['measurement']['source_sha256'] == digest(ROOT / 'benchmarks/measure_command.c')):
            raise AssertionError('evidence validation failed')
    return metadata, rows, stats


def label(case):
    reads = f'{case["n_reads"] // 1000}k' if case['n_reads'] < 1000000 else '1M'
    samples = len(case['reads'])
    return f'{reads} / {samples} sample' + ('s' if samples > 1 else '')


def render(metadata, rows, stats):
    controlled = [s for s in stats if s['case']['scope'] == 'controlled_simulation']
    mismatch = [s for s in controlled if s['k'] == 1]
    lower, upper = min(s['speedup'] for s in mismatch), max(s['speedup'] for s in mismatch)
    memory_low, memory_high = min(s['memory_ratio'] for s in mismatch), max(s['memory_ratio'] for s in mismatch)
    lines = [
        '# DotMatch versus guide-counter', '',
        f'**DotMatch is {lower:.1f}–{upper:.1f}× faster with {memory_low:.1f}–{memory_high:.1f}× lower peak memory in the tested one-mismatch workflows.**', '',
        'These are complete single-thread counting commands against the same 87,437-guide Yusa library, with identical full guide/sample count matrices. The performance workloads use controlled simulated FASTQs with 100,000 or 1,000,000 reads. They measure counting performance under the recorded conditions; they do not establish genome-wide biological accuracy or performance on full experimental screens.', '',
        'The benchmarked source includes the compatibility fixes and optimizations described here. Measurements apply to the recorded source and binary hashes; they do not benchmark a published wheel.', '',
        '## Complete-command throughput', '',
        '![Complete counting throughput for DotMatch and guide-counter, with five-run ranges](../../../benchmarks/figures/guide_counter_throughput.svg)', '',
        'Bars show median throughput; whiskers show the observed minimum and maximum across five fresh processes. Exact and one-mismatch modes are reported separately. All timings include library/index construction, offset detection, gzip parsing, counting, and the count, extended-count and statistics files.', '',
        '| Reads / samples | Mode | DotMatch seconds | guide-counter seconds | Paired speedup, median [min–max] |',
        '| --- | --- | --- | --- | --- |',
    ]
    for s in controlled:
        lines.append(f'| {label(s["case"])} | {"Exact" if s["k"] == 0 else "One mismatch"} | {s["dotmatch"]["seconds"]:.3f} | {s["guide_counter"]["seconds"]:.3f} | {s["speedup"]:.2f}× [{s["speedup_min"]:.2f}–{s["speedup_max"]:.2f}] |')
    lines += ['', 'Speedups are medians of within-repeat guide-counter/DotMatch runtime ratios. Guide-counter is faster in the recorded exact-mode million-read comparisons; the one-mismatch headline does not apply to exact mode.', '',
              '## Peak memory', '',
              '![Peak resident memory for the same complete commands](../../../benchmarks/figures/guide_counter_memory.svg)', '',
              'Memory is per-command peak resident set size from a fresh native launcher and `wait4`, converted from KiB to MiB. The launcher prevents the benchmark Python process’s validation heap from inflating child memory. Bars show medians and whiskers show the five-run range.', '',
              '| Reads / samples | Mode | DotMatch maximum MiB | guide-counter maximum MiB | Paired median memory ratio |',
              '| --- | --- | --- | --- | --- |']
    for s in controlled:
        lines.append(f'| {label(s["case"])} | {"Exact" if s["k"] == 0 else "One mismatch"} | {s["dotmatch"]["peak_rss_kib_max"] / 1024:.1f} | {s["guide_counter"]["peak_rss_kib_max"] / 1024:.1f} | {s["memory_ratio"]:.2f}× |')
    lines += ['', '## What changed', '',
              'The compatibility entrypoint now matches guide-counter’s counting rules for supported ACGT libraries of 1–32 bases:', '',
              '- Offset detection uses the same exact/one-mismatch lookup as counting.',
              '- Offset fractions use all matched windows as their denominator; offsets with zero matches are excluded, and no match means no fallback offset.',
              '- Each selected offset is counted independently. One read may therefore contribute more than one guide count.',
              '- Offset fractions retain full double precision through CLI forwarding; boundary regressions cover values immediately below, equal to and above 0.5.',
              '- Uppercase ACGT read windows are required. Windows containing `N`, IUPAC bytes or lowercase bases are excluded in this entrypoint.',
              '- Exact hits take precedence; equal-distance ties are excluded.', '',
              'The implementation uses two packed Hamming seeds to find possible distance-one targets, verifies their full packed distance, and avoids materializing every possible mutated guide. Rolling offset encoding and table-based base decoding reduce repeated work. The library is parsed once, and all three outputs are written directly from the count matrix.', '',
              'Ordinary `dotmatch count` retains its own documented policies, including literal unknown-byte matching. For guide-counter comparisons use the compatibility entrypoint and explicit shared offset parameters.', '',
              '## Count agreement and experimental smoke check', '',
              f'All {len(rows)} recorded commands completed successfully. Every guide/sample count agrees across tools and repeats, as do guide annotations, extended counts and numerical statistics. Differences in floating-point text formatting are allowed; integer counts are exact.', '',
              'Independent byte-oracle regressions cover offset detection, threshold boundaries, exact priority, ties, multiple windows, unknown/lowercase rejection, empty inputs, and guide lengths 1, 19, 20 and 32. Duplicate guide sequences are rejected.', '']
    experimental = [s for s in stats if s['case']['scope'] == 'experimental_prefix']
    for s in experimental:
        if s['k'] == 1:
            counts = s['case']['records_by_sample']
            lines.append(f'The cached ERR376998/ERR376999 experimental prefixes contain {counts[0]} and {counts[1]} reads. Their full count matrices agree in both modes. This small check verifies compatibility; its timings are excluded from the headline and figures because startup dominates. It is not a full experimental-read performance benchmark.')
    lines += ['', 'Earlier guide-counter-style comparisons in [the historical public CRISPR report](../crispr_comparison/README.md) used another counting path and reported count differences. This corrected compatibility benchmark is separate; those earlier results have not been relabeled as agreeing.', '',
              '## Protocol and provenance', '',
              '- Controlled recipe: 50-base uppercase ACGT reads; uniformly sampled library guides at offsets 23, 24 or 25; 70% exact, 20% one substitution, 5% two substitutions, 5% random decoys. Construction labels are not guaranteed biological source identities.',
              '- Read sets: 100,000 in one FASTQ; 1,000,000 in one FASTQ; 1,000,000 across four FASTQs. Generation seeds and input hashes are frozen in the JSON protocol.',
              '- Offset sample: first 100,000 reads per FASTQ, or the whole file when shorter; minimum matched-window fraction 0.0025.',
              '- Five paired repeats alternate tool order. Both tools run one thread. Inputs are reused from the local filesystem; cache flushing and cold-storage performance are not measured.',
              '- Fresh processes include startup and indexing. Generation, output checks and plotting are outside the timed region.',
              f'- Competitor: unmodified [guide-counter 0.1.3 source](https://github.com/fulcrumgenomics/guide-counter/tree/{metadata["guide_counter_commit"]}), commit `{metadata["guide_counter_commit"]}`; locked release build with upstream fat LTO and one codegen unit.',
              f'- Toolchain: `{metadata["rustc_version"]}`; `{metadata["cc_version"]}`. DotMatch uses its default `-O3` build and `-mavx2` on this x86-64 host.',
              f'- Host: `{metadata["platform"]}`; CPU `{metadata["cpu_model"]}`. Container scheduling and shared-host load remain sources of timing variation.', '',
              '## Reproduce', '',
              'Build the pinned competitor without modifying its source. Use Rust 1.95.0 for the recorded build:', '',
              '```bash', 'git clone https://github.com/fulcrumgenomics/guide-counter.git benchmarks/work/guide-counter-source',
              f'git -C benchmarks/work/guide-counter-source checkout {metadata["guide_counter_commit"]}',
              'cargo build --release --locked --manifest-path benchmarks/work/guide-counter-source/Cargo.toml',
              'python3 scripts/fetch_mageck_demo.py --out benchmarks/real/data --subsample 25',
              'make dotmatch CC=cc',
              'python3 scripts/bench_guide_counter.py \\',
              '  --guide-counter benchmarks/work/guide-counter-source/target/release/guide-counter \\',
              '  --guide-counter-source benchmarks/work/guide-counter-source \\',
              '  --rustc "$(command -v rustc)" \\',
              '  --sizes 100000,1000000 --repeats 5 \\',
              '  --experimental-fastq benchmarks/real/data/ERR376998.fastq.gz \\',
              '  --experimental-fastq benchmarks/real/data/ERR376999.fastq.gz',
              'python3 scripts/report_guide_counter.py',
              'python3 scripts/report_guide_counter.py --check --check-source',
              'make test cli-test guide-counter-gate', '```', '',
              'For an existing protocol, pass `--reuse-inputs` with the same `--work-dir`; all input hashes are checked before reuse. Do not run another benchmark or CPU-intensive test suite concurrently with timing.', '',
              'Raw [timings and agreement](https://github.com/dnncha/dotmatch/blob/main/benchmarks/raw/guide_counter_comparison.csv), [source/input protocol](https://github.com/dnncha/dotmatch/blob/main/benchmarks/raw/guide_counter_comparison.json), and [audit/figure hashes](https://github.com/dnncha/dotmatch/blob/main/benchmarks/raw/guide_counter_comparison_validation.json) are retained. SVG, PNG and PDF versions of both figures are available under `benchmarks/figures/`.', '',
              '## Limits and next validation', '',
              'This comparison covers one public library, a fixed controlled error mixture, two read scales and one/four samples on one host. More reads amortize guide-counter’s index construction, so the ratio can change substantially for full screens. Results do not cover indels, arbitrary regex compatibility, guides longer than 32 bases, full experimental FASTQs or biological gene-hit accuracy. A lower runtime or more assigned windows is not evidence of better biological inference.', '',
              'Larger experimental Yusa and Sanson/Brunello FASTQs remain the next benchmark. Downloads from ENA were blocked by this cloud environment’s host policy during this run; the required host additions were saved for review. The available 25-read prefixes have not been expanded or presented as full screens.', '']
    return '\n'.join(lines)


def figures(stats):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'svg.hashsalt': 'dotmatch-guide-counter',
                               'axes.spines.top': False, 'axes.spines.right': False})
    controlled = [s for s in stats if s['case']['scope'] == 'controlled_simulation']
    names = list(dict.fromkeys(s['workload'] for s in controlled))
    written = []
    for metric, stem, title, scale, ylabel in (
        ('reads_per_second', 'guide_counter_throughput', 'Complete guide counting', 1e6, 'Million reads / second'),
        ('peak_rss_kib', 'guide_counter_memory', 'Peak resident memory', 1024, 'Peak RSS (MiB)'),
    ):
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
        fig.suptitle(title, fontsize=17, fontweight='bold', y=.98)
        for k, ax in enumerate(axes):
            subset = [next(s for s in controlled if s['workload'] == name and s['k'] == k) for name in names]
            for i, (tool, color, display) in enumerate((('guide_counter', '#6b7280', 'guide-counter 0.1.3'),
                                                       ('dotmatch', '#0072b2', 'DotMatch'))):
                values = [s[tool][metric] / scale for s in subset]
                lower = [v - s[tool][metric + '_min'] / scale for v, s in zip(values, subset)]
                upper = [s[tool][metric + '_max'] / scale - v for v, s in zip(values, subset)]
                bars = ax.bar([j + (i - .5) * .36 for j in range(len(subset))], values, .34, color=color,
                              label=display, yerr=[lower, upper], capsize=3,
                              error_kw={'elinewidth': 1, 'ecolor': '#374151'})
                ax.bar_label(bars, labels=[f'{v:.2f}' if metric == 'reads_per_second' else f'{v:.0f}' for v in values],
                             padding=5, fontsize=10)
            ax.set_title('Exact matching' if k == 0 else 'One mismatch, no indels', fontsize=13)
            ax.set_xticks(range(len(subset)), [label(s['case']).replace(' / ', '\n') for s in subset])
            ax.set_ylabel(ylabel)
            ax.set_axisbelow(True)
            ax.grid(axis='y', alpha=.2)
            ax.set_ylim(0, ax.get_ylim()[1] * 1.18)
        axes[0].legend(frameon=False, loc='upper left', fontsize=10)
        fig.text(.5, .015, 'Same 87,437-guide Yusa library · controlled simulated FASTQs · one thread · five alternating repeats\n'
                 'Identical full count matrices · median bars; observed min–max whiskers · experimental smoke timings excluded',
                 ha='center', fontsize=10, color='#374151')
        fig.tight_layout(rect=(0, .09, 1, .92))
        for extension in ('svg', 'png', 'pdf'):
            path = FIGURES / f'{stem}.{extension}'
            metadata = {'Date': None} if extension == 'svg' else ({'CreationDate': None, 'ModDate': None} if extension == 'pdf' else {})
            fig.savefig(path, dpi=200, metadata=metadata)
            written.append(path)
        plt.close(fig)
    return {path.name: digest(path) for path in written}, matplotlib.__version__


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--check-source', action='store_true')
    args = parser.parse_args()
    metadata, rows, stats = load(args.check_source)
    content = render(metadata, rows, stats)
    validation = RAW / f'{STEM}_validation.json'
    if args.check:
        check_callout_headlines(stats)
        if not (REPORT.read_text() == content):
            raise AssertionError('report differs from raw evidence')
        recorded = json.loads(validation.read_text())
        check_validation_record(recorded, rows)
        if not (recorded['report_sha256'] == digest(REPORT)):
            raise AssertionError('evidence validation failed')
        if not (recorded['csv_sha256'] == metadata['csv_sha256']):
            raise AssertionError('evidence validation failed')
        if not (recorded['report_script_sha256'] == digest(__file__)):
            raise AssertionError('report generator changed')
        for name, sha in recorded['figures'].items():
            if not (digest(FIGURES / name) == sha):
                raise AssertionError(f'figure changed: {name}')
    else:
        REPORT.parent.mkdir(parents=True, exist_ok=True)
        REPORT.write_text(content)
        hashes, version = figures(stats)
        validation.write_text(json.dumps({'schema_version': 1, 'commands': len(rows), 'count_mismatches': 0,
                                          'extended_mismatches': 0, 'stats_mismatches': 0,
                                          'csv_sha256': metadata['csv_sha256'], 'report_sha256': digest(REPORT),
                                          'report_script_sha256': digest(__file__), 'matplotlib_version': version,
                                          'figures': hashes}, indent=2) + '\n')
    print(f'guide-counter evidence: {len(rows)} complete commands; paired counts and outputs agree')


if __name__ == '__main__':
    main()
