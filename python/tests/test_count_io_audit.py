"""Lossless import regressions: malformed input must not become scientific data."""
import gzip
import pytest
from dotmatch.count_io import parse_count_value, read_count_table, read_crispr_count_matrix


def test_counts_preserve_integers_above_float_precision(tmp_path):
    path = tmp_path / 'counts.tsv'
    path.write_text('sgRNA\tGene\tA\tB\n0007\tNA\t9007199254740993\t0\n')
    result = read_crispr_count_matrix(path)
    assert result['guides'] == [{'id': '0007', 'gene': 'NA', 'counts': {'A': 9007199254740993, 'B': 0}}]


@pytest.mark.parametrize('text', ['', ' ', 'NA', 'nan', 'inf', '-1', '1.5', '1e999999999', '0.00001', 'oops'])
def test_missing_invalid_fractional_or_unbounded_values_are_rejected(text):
    with pytest.raises(ValueError):
        parse_count_value(text, 'guide', 'sample')


@pytest.mark.parametrize('text,expected', [('1.0', 1), ('1e3', 1000), ('0', 0), ('9007199254740993.0', 9007199254740993)])
def test_integral_decimal_notation_is_exact(text, expected):
    assert parse_count_value(text, 'guide', 'sample') == expected


@pytest.mark.parametrize('body', ['g\tG\t1\n', 'g\tG\t1\t\n', 'g\tG\t1\t2\textra\n', 'g\tG\t1\t2\ng\tG\t3\t4\n', '\tG\t1\t2\n'])
def test_ragged_missing_and_duplicate_rows_fail(tmp_path, body):
    path = tmp_path / 'counts.tsv'
    path.write_text('sgRNA\tGene\tA\tB\n' + body)
    with pytest.raises(ValueError):
        read_crispr_count_matrix(path)


def test_selected_samples_retain_requested_order_and_gz_bom(tmp_path):
    path = tmp_path / 'counts.TSV.GZ'
    path.write_bytes(gzip.compress('\ufeffsgRNA\tGene\tA\tB\ng\tNA\t7\t9\n'.encode()))
    table = read_count_table(path, sample_cols=['B', 'A'])
    assert table.sample_names == ('B', 'A')
    assert table.counts == ((9, 7),)
    assert table.metadata['Gene'] == ('NA',)


def test_dotmatch_detailed_output_does_not_double_count_components(tmp_path):
    path = tmp_path / 'counts.tsv'
    path.write_text('target_id\ttarget_seq\tgene\tambiguous_nearby\tA_count_exact\tA_count_total\tB_count_total\ng\tACGT\tNA\t1\t7\t10\t0\n')
    table = read_count_table(path)
    assert table.sample_names == ('A', 'B')
    assert table.counts == ((10, 0),)
    assert table.metadata['target_seq'] == ('ACGT',)
    assert read_count_table(path, sample_cols=['B_count_total']).counts == ((0,),)


@pytest.mark.parametrize('selection', [[], ['missing'], ['A', 'A']])
def test_bad_sample_selection_fails(tmp_path, selection):
    path = tmp_path / 'counts.tsv'
    path.write_text('sgRNA\tGene\tA\ng\tG\t1\n')
    with pytest.raises(ValueError, match='sample_cols'):
        read_count_table(path, sample_cols=selection)
