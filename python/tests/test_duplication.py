import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest
from dotmatch.duplication import audit_duplication


def fastq(path, sequences, side=None):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt') as f:
        for i, seq in enumerate(sequences):
            suffix = f'/{side}' if side else ''
            f.write(f'@r{i}{suffix}\n{seq}\n+\n{"I" * len(seq)}\n')
    return path


def test_exact_counts_and_preserved_input(tmp_path):
    p = fastq(tmp_path / 'r.gz', ['ACGT', 'acgt', 'ACG', 'ACGN', 'ACGN', 'TGCA'])
    before = p.read_bytes()
    r = audit_duplication(p, cache_mb=1, temp_dir=tmp_path)
    assert (r['records'], r['distinct_sequences'], r['repeated_records']) == (6, 4, 2)
    assert r['singleton_sequences'] == 2
    assert r['largest_sequence_group'] == 2
    assert r['records_with_N'] == 2
    assert r['pcr_duplicate_fraction'] is None
    assert r['recommended_action'] == 'preserve_read_counts'
    assert p.read_bytes() == before
    assert not list(tmp_path.glob('dotmatch-duplication-*'))


def test_both_mates_define_identity(tmp_path):
    a = fastq(tmp_path / 'a.fq', ['AAA', 'AAA', 'AAA'], 1)
    b = fastq(tmp_path / 'b.fq', ['CCC', 'GGG', 'CCC'], 2)
    r = audit_duplication(a, reads2=b)
    assert r['unit'] == 'read_pairs'
    assert r['distinct_sequences'] == 2
    assert r['repeated_records'] == 1


@pytest.mark.parametrize('mode', ['length', 'id', 'malformed'])
def test_fail_closed_and_clean(tmp_path, mode):
    a = fastq(tmp_path / 'a.fq', ['AAA', 'CCC'], 1)
    b = fastq(tmp_path / 'b.fq', ['AAA'] if mode == 'length' else ['AAA', 'CCC'], 2)
    if mode == 'id':
        b.write_text(b.read_text().replace('@r1', '@other'))
    if mode == 'malformed':
        b.write_text('@r0/2\nAAA\n+\nII\n')
    with pytest.raises(ValueError):
        audit_duplication(a, reads2=b, temp_dir=tmp_path)
    assert not list(tmp_path.glob('dotmatch-duplication-*'))


def test_empty_and_prefix(tmp_path):
    p = fastq(tmp_path / 'r.fq', [])
    assert audit_duplication(p)['repeated_fraction'] is None
    fastq(p, ['AAA', 'AAA', 'CCC'])
    r = audit_duplication(p, max_records=2)
    assert r['scope'] == 'prefix' and r['records'] == 2
    assert r['repeated_fraction'] == 0.5


@pytest.mark.parametrize('kwargs', [{'cache_mb': 0}, {'max_records': 0}, {'cache_mb': -1}])
def test_invalid_options(tmp_path, kwargs):
    with pytest.raises(ValueError):
        audit_duplication(tmp_path / 'missing', **kwargs)


def test_cli_routes_without_native_delegation(tmp_path):
    p = fastq(tmp_path / 'r.fq', ['AAA', 'AAA'])
    command = [sys.executable, '-m', 'dotmatch.cli', 'duplication', '--reads', str(p), '--json']
    run = subprocess.run(command, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    assert json.loads(run.stdout)['repeated_records'] == 1
    p.write_text('@bad\nAAA\n+\nI\n')
    run = subprocess.run(command, capture_output=True, text=True)
    assert run.returncode == 2
    assert json.loads(run.stdout)['status'] == 'error'
    assert run.stderr == ''


def test_order_invariant_and_quality_independent(tmp_path):
    sequences = ['AAA', 'CCC', 'AAA', 'NNN', 'TTT', 'NNN']
    a = fastq(tmp_path / 'a.fq', sequences)
    b = fastq(tmp_path / 'b.fq', list(reversed(sequences)))
    b.write_text(b.read_text().replace('III', '!!!'))
    x, y = audit_duplication(a), audit_duplication(b)
    for key in ['records', 'distinct_sequences', 'repeated_records', 'singleton_sequences', 'largest_sequence_group']:
        assert x[key] == y[key]
