"""Real optional-stack integration; CI installs and exercises these dependencies."""
import pytest
import dotmatch

pd = pytest.importorskip('pandas')


def test_dataframe_named_columns_and_labels_are_honoured():
    reads = pd.DataFrame({'sequence': ['acgt', 'tTTt'], 'read_id': ['0001', 'NA']})
    targets = pd.DataFrame({'sequence': ['ACGT', 'TTTT'], 'target_id': ['g', 'h']})
    result = dotmatch.assign_dataframe(reads, targets, k=0, metric='exact')
    assert result.target_name.tolist() == ['g', 'h']
    assert result.read_id.tolist() == ['0001', 'NA']
    assert reads.sequence.tolist() == ['acgt', 'tTTt']


def test_explicit_dataframe_columns_and_integer_indexes():
    reads = pd.DataFrame({'observed': ['ACGT'], 'label': ['r']})
    targets = pd.DataFrame({'expected': ['acgt'], 'label': ['g']})
    result = dotmatch.assign_dataframe(reads, targets, k=0, metric='exact', read_seq_col='observed',
                                       read_id_col='label', target_seq_col='expected', target_id_col='label')
    assert result.target_name.tolist() == ['g']
    assert result.read_id.tolist() == ['r']


@pytest.mark.parametrize('bad', [None, float('nan'), pd.NA, 123])
def test_missing_sequences_are_not_stringified(bad):
    with pytest.raises(ValueError):
        dotmatch.targets_from_dataframe(pd.DataFrame({'id': ['g'], 'sequence': [bad]}))


def test_duplicate_target_ids_and_label_length_mismatch_fail():
    with pytest.raises(ValueError, match='duplicate'):
        dotmatch.targets_from_dataframe(pd.DataFrame({'id': ['g', 'g'], 'sequence': ['ACGT', 'TTTT']}))
    with pytest.raises(ValueError, match='read_ids'):
        dotmatch.assign_dataframe(['ACGT'], ['ACGT'], read_ids=[])
    with pytest.raises(ValueError, match='target_names'):
        dotmatch.assign_dataframe(['ACGT'], ['ACGT'], target_names=[])


def test_numpy_label_array_has_no_ambiguous_truth_value():
    import numpy as np
    result = dotmatch.assign_dataframe(pd.Series(['ACGT', 'TTTT']), ['ACGT', 'TTTT'],
                                       k=0, metric='exact', read_ids=np.array(['a', 'b']))
    assert result.read_id.tolist() == ['a', 'b']


def test_empty_assignment_dataframe_keeps_its_schema():
    result = dotmatch.assign_dataframe([], ['ACGT'], k=0, metric='exact', read_ids=[])
    assert {'status', 'status_name', 'target_name', 'read_id'} <= set(result.columns)
    assert len(result) == 0


def test_polars_does_not_require_a_pyarrow_roundtrip():
    pl = pytest.importorskip('polars')
    targets = pl.DataFrame({'sequence': ['acgt'], 'id': ['g']})
    reads = pl.DataFrame({'read_id': ['r'], 'sequence': ['acgt']})
    assert dotmatch.targets_from_dataframe(targets) == [('g', 'ACGT')]
    assert dotmatch.assign_dataframe(reads, targets, k=0, metric='exact').target_name.tolist() == ['g']


def test_counts_anndata_is_sparse_exact_and_selects_samples(tmp_path):
    pytest.importorskip('anndata')
    from scipy import sparse
    path = tmp_path / 'counts.tsv'
    path.write_text('sgRNA\tGene\tA\tB\n0001\tNA\t9007199254740993\t0\ng2\tG\t0\t3\n')
    adata = dotmatch.counts_tsv_to_anndata(path, sample_cols=['B', 'A'])
    assert sparse.isspmatrix_csr(adata.X)
    assert adata.X.dtype.name == 'int64'
    assert adata.obs_names.tolist() == ['B', 'A']
    assert adata.var_names.tolist() == ['0001', 'g2']
    assert adata.X[1, 0] == 9007199254740993
    assert adata.var.loc['0001', 'Gene'] == 'NA'
    out = tmp_path / 'roundtrip.h5ad'
    adata.write_h5ad(out)
    import anndata
    assert anndata.read_h5ad(out).X[1, 0] == 9007199254740993


def test_count_import_has_no_permissive_fallback(tmp_path):
    pytest.importorskip('anndata')
    path = tmp_path / 'bad.tsv'
    for value in ['', '-1', '0.5', 'NaN']:
        path.write_text(f'sgRNA\tGene\tA\ng\tG\t{value}\n')
        with pytest.raises(ValueError):
            dotmatch.counts_tsv_to_anndata(path)
    path.write_text(f'sgRNA\tGene\tA\ng\tG\t{2**63}\n')
    with pytest.raises(OverflowError, match='int64'):
        dotmatch.counts_tsv_to_anndata(path)


def test_assignment_matrix_keeps_zero_count_axes_and_does_not_mutate_input(tmp_path):
    pytest.importorskip('anndata')
    from scipy import sparse
    frame = pd.DataFrame({'cell_barcode': ['0001', 'NA', 'NA', 'empty'],
                          'target_name': ['', 'g', 'g', ''],
                          'status_name': ['ambiguous', 'unique', 'unique', 'none']})
    before = frame.copy(deep=True)
    adata = dotmatch.assignments_to_anndata(frame, feature_names=['unused', 'g'],
                cell_names=['empty', 'NA', '0001', 'unobserved'], include_ambiguous_per_cell=True)
    assert sparse.isspmatrix_csr(adata.X)
    assert adata.shape == (4, 2)
    assert adata.X.sum() == 2 and adata.X[1, 1] == 2
    assert adata.obs.loc['0001', 'ambiguous_count'] == 1
    assert adata.obs.loc['empty', 'n_observations'] == 1
    assert adata.uns['dotmatch']['umi_deduplicated'] is False
    pd.testing.assert_frame_equal(frame, before)
    path = tmp_path / 'assignments.tsv'
    frame.to_csv(path, sep='\t', index=False)
    from_path = dotmatch.assignments_to_anndata(path)
    assert from_path.obs_names.tolist() == ['0001', 'NA', 'empty']


def test_cells_are_never_invented_from_read_ids():
    pytest.importorskip('anndata')
    frame = pd.DataFrame({'read_id': ['madeup_cell_read'], 'target_name': ['g'], 'status': [1]})
    with pytest.raises(ValueError, match='never inferred'):
        dotmatch.assignments_to_anndata(frame)
    assert 'cell_barcode' not in frame


@pytest.mark.parametrize('changes', [{'cell_names': ['other']}, {'feature_names': ['other']},
                                      {'count_unique_only': False}, {'cell_names': ['c', 'c']}])
def test_assignment_conversion_never_silently_drops_counts(changes):
    pytest.importorskip('anndata')
    frame = pd.DataFrame({'cell_barcode': ['c'], 'target_name': ['g'], 'status': [1]})
    with pytest.raises(ValueError):
        dotmatch.assignments_to_anndata(frame, **changes)


def test_unknown_status_fails_and_all_unassigned_cells_survive():
    pytest.importorskip('anndata')
    frame = pd.DataFrame({'cell_barcode': ['c', 'd'], 'target_name': ['', ''], 'status': [2, 0]})
    assert dotmatch.assignments_to_anndata(frame).shape == (2, 0)
    frame.loc[0, 'status'] = 99
    with pytest.raises(ValueError, match='unknown assignment status'):
        dotmatch.assignments_to_anndata(frame)


def test_explicit_integer_column_positions_use_original_dataframe_order():
    reads = pd.DataFrame({'observed': ['ACGT'], 'label': ['r']})
    targets = pd.DataFrame({'expected': ['acgt'], 'label': ['g']})
    result = dotmatch.assign_dataframe(reads, targets, k=0, metric='exact',
        read_seq_col=0, read_id_col=1, target_seq_col=0, target_id_col=1)
    assert result.read_id.tolist() == ['r'] and result.target_name.tolist() == ['g']
    with pytest.raises(ValueError):
        dotmatch.targets_from_dataframe(targets, id_col=0, seq_col=0)


def test_target_path_is_supported_but_scalar_read_text_is_not(tmp_path):
    path = tmp_path / 'targets.tsv'
    path.write_text('target_id\tsequence\ng\tACGT\n')
    assert dotmatch.assign_dataframe(['acgt'], path, k=0, metric='exact').target_name.tolist() == ['g']
    with pytest.raises(TypeError):
        dotmatch.assign_dataframe('ACGT', ['ACGT'], k=0, metric='exact')


def test_large_sparse_axes_do_not_allocate_a_dense_cell_feature_rectangle():
    pytest.importorskip('anndata')
    from scipy import sparse
    size = 50000
    frame = pd.DataFrame({'cell_barcode': ['c0'], 'target_name': ['g0'], 'status': [1]})
    adata = dotmatch.assignments_to_anndata(frame,
        cell_names=[f'c{i}' for i in range(size)], feature_names=[f'g{i}' for i in range(size)])
    assert adata.shape == (size, size) and sparse.isspmatrix_csr(adata.X)
    assert adata.X.nnz == 1
    assert adata.X.data.nbytes + adata.X.indices.nbytes + adata.X.indptr.nbytes < 1000000
