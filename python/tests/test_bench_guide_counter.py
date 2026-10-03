"""Benchmark validator rejects malformed complete-command artifacts."""

import csv
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('counter_benchmark', ROOT / 'scripts/bench_guide_counter.py')
BENCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BENCH)


def artifacts(tmp_path):
    prefix = tmp_path / 'run'
    Path(str(prefix) + '.counts.txt').write_text('guide\tgene\tsample\ng\tG\t1\n')
    Path(str(prefix) + '.extended-counts.txt').write_text('guide\tgene\tguide_type\tsample\ng\tG\tother\t1\n')
    stats = {field: '0' for field in sorted(BENCH.STATS_FIELDS)}
    stats.update(file='reads.fastq', label='sample', total_guides='1', total_reads='1', mapped_reads='1')
    with Path(str(prefix) + '.stats.txt').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(stats), delimiter='\t')
        writer.writeheader()
        writer.writerow(stats)
    return prefix


@pytest.mark.parametrize('value', ['nan', 'inf', '-inf'])
def test_nonfinite_statistics_cannot_disappear_in_comparison(tmp_path, value):
    prefix = artifacts(tmp_path)
    before, after = BENCH.outputs(prefix), BENCH.outputs(prefix)
    after['stats'][0]['mean_reads_per_guide'] = value
    with pytest.raises(ValueError, match='Nonfinite'):
        BENCH.compare(before, after)
    path = Path(str(prefix) + '.stats.txt')
    path.write_text(path.read_text().replace('\t0', f'\t{value}', 1))
    with pytest.raises(ValueError, match='Nonfinite'):
        BENCH.outputs(prefix)


def test_duplicate_extended_rows_are_rejected(tmp_path):
    prefix = artifacts(tmp_path)
    path = Path(str(prefix) + '.extended-counts.txt')
    with path.open('a') as stream:
        stream.write('g\tG\tother\t1\n')
    with pytest.raises(ValueError, match='Duplicate extended'):
        BENCH.outputs(prefix)


def test_missing_statistic_columns_and_samples_are_rejected(tmp_path):
    prefix = artifacts(tmp_path)
    before, after = BENCH.outputs(prefix), BENCH.outputs(prefix)
    after['stats'][0].pop('mean_reads_per_guide')
    with pytest.raises(ValueError, match='fields differ'):
        BENCH.compare(before, after)
    after['stats'] = []
    with pytest.raises(ValueError, match='sample coverage'):
        BENCH.compare(before, after)


def test_valid_artifacts_compare_without_differences(tmp_path):
    prefix = artifacts(tmp_path)
    result = BENCH.outputs(prefix)
    assert BENCH.compare(result, result) == (0, 0, 0)
