"""Counting failures preserve research inputs and previously published outputs."""
import os
import gzip
import json
from pathlib import Path

import pytest

from dotmatch.entrypoint import main


def count_inputs(tmp_path):
    library = tmp_path / 'targets.tsv'
    library.write_text('target_id\tsequence\ng1\tACGT\n')
    reads = tmp_path / 'reads.fastq'
    reads.write_text('@r1\nACGT\n+\nIIII\n')
    return library, reads


def count_args(library, reads, output):
    return ['count', '--targets', str(library), '--reads', str(reads),
            '--target-start', '0', '--target-length', '4', '--out', str(output)]


def test_failed_python_count_preserves_all_previous_outputs(tmp_path, monkeypatch):
    monkeypatch.setenv('DOTMATCH_PYTHON_NO_DELEGATE', '1')
    library, reads = count_inputs(tmp_path)
    reads.write_text(reads.read_text() + '@truncated\nACGT\n+\n')
    outputs = {name: tmp_path / name for name in ('counts.tsv', 'assignments.tsv', 'summary.json')}
    for path in outputs.values():
        path.write_text('previous ' + path.name)
    before = {path: path.read_bytes() for path in (library, reads, *outputs.values())}
    assert main(count_args(library, reads, outputs['counts.tsv']) + [
        '--assignments', str(outputs['assignments.tsv']), '--summary', str(outputs['summary.json'])]) == 2
    assert {path: path.read_bytes() for path in before} == before
    assert set(tmp_path.iterdir()) == set(before)


def test_failed_assignment_iterator_preserves_existing_file(tmp_path):
    from dotmatch import stream_assign, write_assignments_tsv

    library, reads = count_inputs(tmp_path)
    reads.write_text(reads.read_text() + '@truncated\nACGT\n+\n')
    output = tmp_path / 'assignments.tsv'
    output.write_text('previous assignments')
    with pytest.raises(ValueError):
        write_assignments_tsv(stream_assign(reads, ['ACGT'], target_start=0, target_length=4), output)
    assert output.read_text() == 'previous assignments'
    assert sorted(path.name for path in tmp_path.iterdir()) == ['assignments.tsv', 'reads.fastq', 'targets.tsv']


def test_destination_changed_during_publication_is_preserved(tmp_path, monkeypatch):
    monkeypatch.setenv('DOTMATCH_PYTHON_NO_DELEGATE', '1')
    library, reads = count_inputs(tmp_path)
    counts = tmp_path / 'counts.tsv'
    assignments = tmp_path / 'assignments.tsv'
    counts.write_bytes(b'previous counts')
    assignments.write_bytes(b'previous assignments')
    replace = os.replace

    def change_next_destination(source, destination):
        replace(source, destination)
        if Path(destination) == counts:
            assignments.write_bytes(b'external assignments')

    monkeypatch.setattr(os, 'replace', change_next_destination)
    assert main(count_args(library, reads, counts) + ['--assignments', str(assignments)]) == 2
    assert counts.read_bytes() == b'previous counts'
    assert assignments.read_bytes() == b'external assignments'
    assert set(tmp_path.iterdir()) == {library, reads, counts, assignments}


@pytest.mark.parametrize('role', ['out', 'assignments', 'summary'])
@pytest.mark.parametrize('source_role', ['targets', 'reads'])
@pytest.mark.parametrize('alias_kind', ['same', 'hardlink', 'symlink'])
def test_python_count_refuses_input_aliases(tmp_path, monkeypatch, role, source_role, alias_kind):
    monkeypatch.setenv('DOTMATCH_PYTHON_NO_DELEGATE', '1')
    library, reads = count_inputs(tmp_path)
    outputs = {name: tmp_path / (name + '.txt') for name in ('out', 'assignments', 'summary')}
    for path in outputs.values():
        path.write_bytes(b'previous result')
    source = library if source_role == 'targets' else reads
    alias = tmp_path / 'alias'
    if alias_kind == 'same':
        alias = source
    elif alias_kind == 'hardlink':
        os.link(source, alias)
    else:
        alias.symlink_to(source)
    outputs[role] = alias
    before = {path: path.read_bytes() for path in tmp_path.iterdir()}
    assert main(count_args(library, reads, outputs['out']) + [
        '--assignments', str(outputs['assignments']), '--summary', str(outputs['summary'])]) == 2
    assert {path: path.read_bytes() for path in tmp_path.iterdir()} == before
    if alias_kind == 'symlink':
        assert alias.is_symlink()


@pytest.mark.parametrize('linked', [False, True])
def test_python_count_refuses_output_collisions(tmp_path, monkeypatch, linked):
    monkeypatch.setenv('DOTMATCH_PYTHON_NO_DELEGATE', '1')
    library, reads = count_inputs(tmp_path)
    counts = tmp_path / 'counts.tsv'
    counts.write_bytes(b'previous counts')
    assignments = counts
    if linked:
        assignments = tmp_path / 'assignments.tsv'
        os.link(counts, assignments)
    assert main(count_args(library, reads, counts) + ['--assignments', str(assignments)]) == 2
    assert counts.read_bytes() == assignments.read_bytes() == b'previous counts'


def test_failed_fresh_count_leaves_no_artifacts(tmp_path, monkeypatch):
    monkeypatch.setenv('DOTMATCH_PYTHON_NO_DELEGATE', '1')
    library, reads = count_inputs(tmp_path)
    reads.write_text(reads.read_text() + '@truncated\nACGT\n+\n')
    assert main(count_args(library, reads, tmp_path / 'counts.tsv') + [
        '--assignments', str(tmp_path / 'assignments.tsv'), '--summary', str(tmp_path / 'summary.json')]) == 2
    assert set(tmp_path.iterdir()) == {library, reads}


def test_successful_gzip_reruns_publish_all_outputs(tmp_path, monkeypatch):
    monkeypatch.setenv('DOTMATCH_PYTHON_NO_DELEGATE', '1')
    library, reads = count_inputs(tmp_path)
    counts, assignments, summary = [tmp_path / name for name in ('counts.tsv.gz', 'assignments.tsv.gz', 'summary.json.gz')]
    args = count_args(library, reads, counts) + ['--assignments', str(assignments), '--summary', str(summary)]
    for _ in range(2):
        assert main(args) == 0
        with gzip.open(counts, 'rt') as handle:
            assert handle.read().splitlines()[1].startswith('g1\tACGT\t\t1\t')
        with gzip.open(assignments, 'rt') as handle:
            assert handle.read().splitlines()[1].startswith('r1\tACGT\tg1\t')
        with gzip.open(summary, 'rt') as handle:
            assert json.load(handle)['assigned_unique'] == 1
    assert set(tmp_path.iterdir()) == {library, reads, counts, assignments, summary}


def test_rollback_preserves_foreign_replacement_and_retains_previous_file(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv('DOTMATCH_PYTHON_NO_DELEGATE', '1')
    library, reads = count_inputs(tmp_path)
    counts = tmp_path / 'counts.tsv'
    assignments = tmp_path / 'assignments.tsv'
    counts.write_bytes(b'previous counts')
    assignments.write_bytes(b'previous assignments')
    replace = os.replace

    def fail_after_foreign_replacement(source, destination):
        if Path(destination) == assignments:
            raise OSError('injected publication failure')
        replace(source, destination)
        if Path(destination) == counts:
            foreign = tmp_path / 'foreign'
            foreign.write_bytes(b'external counts')
            replace(foreign, counts)

    monkeypatch.setattr(os, 'replace', fail_after_foreign_replacement)
    assert main(count_args(library, reads, counts) + ['--assignments', str(assignments)]) == 2
    assert counts.read_bytes() == b'external counts'
    assert assignments.read_bytes() == b'previous assignments'
    backups = set(tmp_path.iterdir()) - {library, reads, counts, assignments}
    assert len(backups) == 1
    backup = backups.pop()
    assert backup.read_bytes() == b'previous counts'
    assert str(backup) in capsys.readouterr().err
