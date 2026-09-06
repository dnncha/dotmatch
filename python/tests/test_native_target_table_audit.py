"""Exercise the actual C CLI, not only the Python library reader."""
import csv
import gzip
import io
import subprocess
from pathlib import Path
import pytest
from dotmatch.target_io import read_target_table

ROOT = Path(__file__).resolve().parents[2]


def count(path, tmp_path):
    binary = ROOT / 'dotmatch'
    if not binary.is_file():
        pytest.skip('native CLI must be built (make dotmatch)')
    reads = tmp_path / 'reads.fastq'
    reads.write_text('@read\nACGTACGT\n+\nIIIIIIII\n')
    output = tmp_path / 'counts.tsv'
    return subprocess.run([str(binary), 'count', '--targets', str(path), '--reads', str(reads),
        '--sample-label', 'sample', '--target-start', '0', '--target-length', '8', '--k', '0',
        '--metric', 'hamming', '--format', 'mageck', '--out', str(output)],
        text=True, capture_output=True, timeout=15), output


@pytest.mark.parametrize('suffix', ['.csv', '.csv.gz', '.CSV.GZ', '.tsv', '.tsv.gz'])
def test_native_and_python_read_reordered_quoted_libraries(suffix, tmp_path):
    delimiter = ',' if '.csv' in suffix.lower() else '\t'
    content = io.StringIO()
    writer = csv.writer(content, delimiter=delimiter)
    writer.writerows([['gene', 'sequence', 'target_id'], ['G,1', 'acgtacgt', 'guide,1']])
    data = ('\ufeff' + content.getvalue()).encode()
    path = tmp_path / ('targets' + suffix)
    path.write_bytes(gzip.compress(data) if suffix.lower().endswith('.gz') else data)
    expected = read_target_table(path)
    result, output = count(path, tmp_path)
    assert result.returncode == 0, result.stderr
    row = list(csv.DictReader(output.open(), delimiter='\t'))[0]
    assert row['sgRNA'] == expected[0].target_id == 'guide,1'
    assert row['Gene'] == expected[0].gene == 'G,1'
    assert row['sample'] == '1'


@pytest.mark.parametrize('content', ['ACGTACGT\n', 'sequence\nACGTACGT\n'])
def test_sequence_only_ids_are_consistent(content, tmp_path):
    path = tmp_path / 'targets.tsv'
    path.write_text(content)
    result, output = count(path, tmp_path)
    assert result.returncode == 0, result.stderr
    assert list(csv.DictReader(output.open(), delimiter='\t'))[0]['sgRNA'] == read_target_table(path)[0].target_id == 'target_0'


@pytest.mark.parametrize('content', [
    'id\tsequence\ng\tACGTACGT\textra\n', 'id\tsequence\ng\n',
    'id\tsequence\ng\tACGTACGT\ng\tACGTACGT\n',
    'id\tsequence\tseq\ng\tACGTACGT\tACGTACGT\n',
    'id\tsequence\tmeta\tmeta\ng\tACGTACGT\tx\tx\n',
    'id\tsequence\n\tACGTACGT\n', 'id\tsequence\ng\tACGT ACGT\n',
    'id\tsequence\ng\tACGT\x00ACGT\n', 'id\tsequence\n"unterminated\tACGTACGT\n',
])
def test_malformed_native_libraries_fail_before_output_creation(content, tmp_path):
    path = tmp_path / 'targets.tsv'
    path.write_text(content)
    result, output = count(path, tmp_path)
    assert result.returncode != 0
    assert not output.exists()
    with pytest.raises((ValueError, csv.Error)):
        read_target_table(path)


def test_native_parser_rejects_corrupt_gzip_and_long_rows(tmp_path):
    path = tmp_path / 'targets.tsv.gz'
    path.write_bytes(gzip.compress(b'id\tsequence\ng\tACGTACGT\n')[:-6])
    result, output = count(path, tmp_path)
    assert result.returncode != 0 and not output.exists()
    path = tmp_path / 'targets.tsv'
    path.write_text('id\tsequence\ng\t' + 'A' * (1024 * 1024 + 1) + '\n')
    result, output = count(path, tmp_path)
    assert result.returncode != 0 and '1 MiB' in result.stderr and not output.exists()


def test_dynamic_target_lines_exceed_old_16k_limit_without_truncation(tmp_path):
    path = tmp_path / 'targets.tsv'
    gene = 'G' * 20000
    path.write_text(f'id\tsequence\tgene\ng\tACGTACGT\t{gene}\n')
    result, output = count(path, tmp_path)
    assert result.returncode == 0, result.stderr
    assert list(csv.DictReader(output.open(), delimiter='\t'))[0]['Gene'] == gene
