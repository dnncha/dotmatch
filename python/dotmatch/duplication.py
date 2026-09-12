"""Read-only exact sequence duplication QC, with disk-backed counting."""
from __future__ import annotations

import argparse
import gzip
import json
import sqlite3
import sys
import tempfile
from contextlib import closing
from itertools import zip_longest
from pathlib import Path


def _records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='ascii') as handle:
        index = 0
        while True:
            header = handle.readline()
            if not header:
                return
            index += 1
            seq, plus, qual = (handle.readline().rstrip('\r\n') for _ in range(3))
            if (not header.startswith('@') or not header[1:].strip()
                    or not plus.startswith('+') or not seq or len(seq) != len(qual)
                    or any(c not in 'ACGTNacgtn' for c in seq)
                    or any(not 33 <= ord(c) <= 126 for c in qual)):
                raise ValueError(f'{path}: malformed four-line FASTQ record {index}')
            yield header[1:].split()[0], seq.upper()


def _pair_id(identifier):
    return identifier[:-2] if identifier.endswith(('/1', '/2')) else identifier


def audit_duplication(reads, *, reads2=None, cache_mb=64, temp_dir=None, max_records=None):
    """Count full-sequence identities, retaining abundance; never infer molecules.

    cache_mb configures SQLite's page cache, not a hard process RSS limit.
    max_records audits only the input prefix; it is not a random sample.
    """
    if not isinstance(cache_mb, int) or cache_mb < 1:
        raise ValueError('cache_mb must be a positive integer')
    if max_records is not None and (not isinstance(max_records, int) or max_records < 1):
        raise ValueError('max_records must be a positive integer')
    first = _records(reads)
    second = _records(reads2) if reads2 is not None else None
    total = ambiguous = 0
    # Full strings, not hashes: no hash collision can merge distinct sequences.
    with tempfile.TemporaryDirectory(prefix='dotmatch-duplication-', dir=temp_dir) as tmp:
        with closing(sqlite3.connect(str(Path(tmp) / 'counts.sqlite'))) as db:
            db.execute(f'PRAGMA cache_size=-{cache_mb * 1024}')
            db.execute('PRAGMA temp_store=FILE')
            db.execute('PRAGMA mmap_size=0')
            db.execute('CREATE TABLE counts (r1 TEXT, r2 TEXT, n INTEGER NOT NULL, PRIMARY KEY(r1,r2)) WITHOUT ROWID')
            records = zip_longest(first, second) if second is not None else ((r, ('', '')) for r in first)
            try:
                for left, right in records:
                    if left is None or right is None:
                        raise ValueError('paired FASTQ files have different record counts')
                    if second is not None and _pair_id(left[0]) != _pair_id(right[0]):
                        raise ValueError(f'paired FASTQ identifiers differ at pair {total + 1}')
                    total += 1
                    ambiguous += int('N' in left[1] or 'N' in right[1])
                    db.execute('INSERT INTO counts VALUES (?, ?, 1) ON CONFLICT(r1,r2) DO UPDATE SET n=n+1', (left[1], right[1]))
                    if total % 10000 == 0:
                        db.commit()
                    if max_records is not None and total >= max_records:
                        break
                db.commit()
                unique, singleton, largest = db.execute('SELECT COUNT(*), COALESCE(SUM(n=1),0), COALESCE(MAX(n),0) FROM counts').fetchone()
            finally:
                first.close()
                if second is not None:
                    second.close()
    excess = total - unique
    return {
        'schema_version': 'dotmatch.duplication.v1',
        'status': 'ok',
        'inputs': {'reads': str(reads), 'reads2': str(reads2) if reads2 is not None else None},
        'scope': 'prefix' if max_records is not None else 'full',
        'max_records': max_records,
        'unit': 'read_pairs' if reads2 is not None else 'reads',
        'identity': 'full_sequence_exact_case_normalized_N_literal_orientation_preserved',
        'records': total, 'distinct_sequences': unique,
        'singleton_sequences': singleton, 'largest_sequence_group': largest,
        'repeated_records': excess,
        'repeated_fraction': excess / total if total else None,
        'records_with_N': ambiguous,
        'backend': 'sqlite_full_sequence_keys', 'sqlite_cache_mb': cache_mb,
        'input_modified': False, 'pcr_duplicate_fraction': None,
        'recommended_action': 'preserve_read_counts',
        'interpretation': 'Sequence repetition is not proof of PCR duplication. Known-target abundance can produce identical reads. For molecule counts, use an assay-specific UMI workflow with the appropriate cell/target or alignment context and UMI error handling.',
        'limitations': ['No PCR duplicate or molecule inference', 'No approximate matching or reverse-complement collapsing', 'Prefix audits are not random samples', 'SQLite cache size is not a hard process memory limit; temporary disk use grows with distinct sequence data'],
    }


def command_duplication(argv):
    parser = argparse.ArgumentParser(prog='dotmatch duplication', description='Audit exact sequence repetition without removing reads. JSON is available for agents.')
    parser.add_argument('--reads', required=True, help='Four-line FASTQ or FASTQ.gz')
    parser.add_argument('--reads2', help='Synchronized mate FASTQ or FASTQ.gz')
    parser.add_argument('--cache-mb', type=int, default=64, help='SQLite page cache in MiB (not a hard RAM limit; default: 64)')
    parser.add_argument('--temp-dir', help='Directory for automatically cleaned disk-backed counts')
    parser.add_argument('--max-records', type=int, help='Audit only the first N reads/pairs, explicitly labelled prefix')
    parser.add_argument('--json', action='store_true', help='Emit one versioned JSON object, including errors')
    args = parser.parse_args(argv)
    try:
        report = audit_duplication(args.reads, reads2=args.reads2, cache_mb=args.cache_mb, temp_dir=args.temp_dir, max_records=args.max_records)
    except (OSError, ValueError, sqlite3.Error, EOFError) as exc:
        if args.json:
            print(json.dumps({'schema_version': 'dotmatch.duplication.v1', 'status': 'error', 'message': str(exc), 'input_modified': False}))
        else:
            print(f'dotmatch duplication: {exc}', file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(report, sort_keys=True))
    else:
        fraction = report['repeated_fraction']
        rate = f'{fraction:.1%}' if fraction is not None else 'not available (empty input)'
        print(f"{report['records']:,} {report['unit']} audited ({report['scope']}); {report['distinct_sequences']:,} distinct sequences; repetition {rate}.")
        print(report['interpretation'])
    return 0
