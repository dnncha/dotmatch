"""Reject ambiguous dataframe columns and malformed assignment tables."""
import pytest
import dotmatch

pd = pytest.importorskip("pandas")


def test_guide_id_name_does_not_conflict_with_canonical_sequence_column():
    targets = pd.DataFrame({'guide': ['g'], 'sequence': ['acgt']})
    assert dotmatch.targets_from_dataframe(targets) == [('g', 'ACGT')]


@pytest.mark.parametrize('text', [
    'cell_barcode\ttarget_name\tstatus\nc\tg\t1\textra\n',
    'cell_barcode\ttarget_name\tstatus\nc\tg\n',
    'cell_barcode\ttarget_name\tstatus\tstatus\nc\tg\t1\t1\n',
])
def test_assignment_table_rejects_ragged_rows_and_duplicate_columns(tmp_path, text):
    from dotmatch.dataframes import _read_assignment_table
    path = tmp_path / 'assignments.tsv'
    path.write_text(text)
    with pytest.raises(ValueError):
        _read_assignment_table(path)


def test_assignment_table_gzip_bom_and_textual_identity(tmp_path):
    import gzip
    from dotmatch.dataframes import _read_assignment_table
    path = tmp_path / 'assignments.TSV.GZ'
    path.write_bytes(gzip.compress('\ufeffcell_barcode\ttarget_name\tstatus\n0001\tNA\t1\n'.encode()))
    frame = _read_assignment_table(path)
    assert frame.iloc[0].tolist() == ['0001', 'NA', '1']
