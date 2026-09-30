#!/usr/bin/env python3
"""AR-001-B historical Brunello prefix audit; see BRUNELLO_FORENSICS.md."""
from __future__ import annotations
import argparse
import csv
import json
import shutil
import subprocess
from pathlib import Path
import pilot
import forensics


def run(out: Path):
    out.mkdir(parents=True, exist_ok=False)
    work = out / 'work'
    work.mkdir()
    summary = {'schema': 'dotmatch.research.AR001B.v1', 'status': 'running',
        'scope': 'first 100000 records per configured sample, calibration included; nonrandom prefixes',
        'baseline': pilot.BASE, 'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=pilot.ROOT, text=True).strip(),
        'protocol_sha256': pilot.digest(pilot.HERE / 'BRUNELLO_FORENSICS.md'),
        'harness_sha256': pilot.digest(__file__), 'samples': {}, 'results': {},
        'limitations': ['No gene-level inference.', 'No whole-archive prevalence or biological accuracy claim.',
                        'Archived FASTQ MD5 values are metadata only, not locally verified for prefixes.']}
    pilot.save_json(out / 'run.json', summary)
    try:
        library, lib_path, provenance = pilot.fetch_library('brunello', work)
        summary['library'] = provenance
        installed = subprocess.check_output(['cargo', 'install', '--list'], text=True)
        pilot.require('guide-counter v0.1.3:' in installed.splitlines(), 'comparator version mismatch')
        summary['guide_counter_version_evidence'] = installed
        summary['guide_counter_binary_sha256'] = pilot.digest(shutil.which('guide-counter'))
        summary['dotmatch_version'] = subprocess.check_output(['dotmatch', '--version'], text=True).strip()
        labels, paths, recsets = [], [], []
        for label, accessions in pilot.SOURCES['brunello']['samples'].items():
            records, metadata = pilot.fetch_prefix(accessions, 100000)
            path = work / f'{label}.fastq.gz'
            pilot.write_fastq(path, records)
            labels.append(label)
            paths.append(path)
            recsets.append(records)
            summary['samples'][label] = {'sources': metadata, 'records': len(records), 'sha256': pilot.digest(path)}
            pilot.save_json(out / 'run.json', summary)
        gc = ['guide-counter', 'count', '--input', *map(str, paths), '--samples', *labels,
              '--library', str(lib_path), '--output', str(out / 'guide_counter'),
              '--offset-sample-size', '100000', '--offset-min-fraction', '0.0025']
        forensics.command(gc, out / 'guide_counter')
        dm = ['dotmatch', 'count', '--targets', str(lib_path), '--sample-label', ','.join(labels),
            '--target-start', '20', '--target-length', '20', '--metric', 'hamming', '--k', '1',
            '--ambiguity-policy', 'best', '--auto-offset', '20', '--auto-offset-sample', '100000',
            '--offset-mode', 'multi', '--offset-min-fraction', '0.0025', '--format', 'mageck',
            '--out', str(out / 'dotmatch.counts.tsv'), '--summary', str(out / 'dotmatch.summary.json')]
        for path in paths:
            dm += ['--reads', str(path)]
        forensics.command(dm, out / 'dotmatch')
        gc_labels, gc_matrix = forensics.matrix(out / 'guide_counter.counts.txt')
        dm_labels, dm_matrix = forensics.matrix(out / 'dotmatch.counts.tsv')
        pilot.require(set(gc_labels) == set(dm_labels) == set(labels), 'sample axes differ')
        pilot.require(gc_matrix.keys() == dm_matrix.keys() == {row[0] for row in library}, 'guide axes differ')
        lookup = forensics.matching_event_lookup(library)
        ref = pilot.Reference(library)
        oracle_checks = 0
        genes = {row[0]: row[2] for row in library}
        for label, records in zip(labels, recsets):
            for _h, sequence, _p, _q in records[:1000]:
                window = sequence[23:43]
                if len(window) == 20 and set(window) <= set('ACGT'):
                    call = ref.calls(window)[2]
                    pilot.require(lookup.get(window, -1) == (call[1] if call[0] == 'unique' else -1), 'lookup and independent oracle disagree')
                    oracle_checks += 1
            predicted, result = forensics.reconstruct(records, library, lookup, label, out)
            mismatches, differences = [], []
            for i, (ident, _sequence, gene) in enumerate(library):
                pilot.require(gc_matrix[ident][0] == dm_matrix[ident][0] == gene, 'annotation differs')
                a, b = gc_matrix[ident][1][label], dm_matrix[ident][1][label]
                if a != predicted[i]:
                    mismatches.append((ident, a, predicted[i]))
                if a != b:
                    differences.append((ident, gene, a, b, b - a))
            pilot.table(out / f'{label}.reconstruction_disagreements.tsv', ['guide', 'guide_counter', 'reconstructed'], mismatches)
            pilot.table(out / f'{label}.tool_count_differences.tsv', ['guide', 'gene', 'guide_counter', 'dotmatch', 'dotmatch_minus_guide_counter'], differences)
            result['reconstruction_count_disagreements'] = len(mismatches)
            result['dotmatch_total'] = sum(row[1][label] for row in dm_matrix.values())
            result['guide_counter_total'] = sum(row[1][label] for row in gc_matrix.values())
            result['differing_guides'] = len(differences)
            same_gene = cross_gene = 0
            with (out / f'{label}.multi_offset_records.tsv').open() as handle:
                for row in csv.DictReader(handle, delimiter='\t'):
                    ids = [event.split(':', 1)[1] for event in row['offset_and_guide'].split(';')]
                    annotations = {genes[ident] for ident in ids}
                    same_gene += len(annotations) == 1
                    cross_gene += len(annotations) > 1
            result['multi_offset_same_gene_reads'] = same_gene
            result['multi_offset_cross_gene_reads'] = cross_gene
            summary['results'][label] = result
            pilot.save_json(out / 'run.json', summary)
            pilot.require(not mismatches, 'comparator event reconstruction failed')
        a = {ident: sum(row[1].values()) for ident, row in dm_matrix.items()}
        b = {ident: sum(row[1].values()) for ident, row in gc_matrix.items()}
        observed = {'dotmatch': sum(a.values()), 'guide_counter': sum(b.values()),
                    'differing_guides': sum(a[ident] != b[ident] for ident in a)}
        summary['historical_aggregate'] = observed
        summary['historical_report_exactly_reproduced'] = observed == {'dotmatch': 349184, 'guide_counter': 350374, 'differing_guides': 255}
        summary['lookup_oracle_checks'] = oracle_checks
        summary['status'] = 'complete'
    except Exception as exc:
        summary['status'] = 'failed'
        summary['error'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        pilot.save_json(out / 'run.json', summary)
        files = {str(p.relative_to(out)): {'bytes': p.stat().st_size, 'sha256': pilot.digest(p)}
                 for p in sorted(out.rglob('*')) if p.is_file() and 'work' not in p.relative_to(out).parts and p.name != 'MANIFEST.json'}
        pilot.save_json(out / 'MANIFEST.json', {'status': summary['status'], 'files': files})
        print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True)
    run(Path(parser.parse_args().out).resolve())
