from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
import subprocess
import sys
import sysconfig

import pytest
import dotmatch


def test_basic_package_import_does_not_load_optional_science_stacks():
    root = str(Path(__file__).resolve().parents[1])
    script = (
        f"import sys;sys.path[:0]={[root, sysconfig.get_path('purelib'), sysconfig.get_path('platlib')]!r};"
        f"import dotmatch;assert dotmatch.__file__.startswith({root!r});"
        "assert dotmatch.distance('ACGT','AGGT')==1;"
        "assert not {'pandas','numpy','scipy','anndata','polars'}.intersection(sys.modules),set(sys.modules)"
    )
    subprocess.run([sys.executable, '-I', '-S', '-c', script], check=True, capture_output=True, text=True)


def test_dataframe_results_do_not_assign_ambiguous_candidate_names():
    pytest.importorskip('pandas')
    results = dotmatch.assign_hamming(['AAAAAAAA', 'TTTTTTTT', 'GGGGGGGG'], ['AAAAAAAA', 'CAAAAAAA', 'TTTTTTTT'], k=1, policy='radius')
    frame = dotmatch.results_to_dataframe(results, target_names=['a', 'b', 'c'])
    assert frame['target_name'].tolist() == ['', 'c', '']
    assert frame['status_name'].tolist() == ['ambiguous', 'unique', 'none']


def test_optional_dataframe_result_api_still_works():
    pytest.importorskip('pandas')
    frame = dotmatch.assign_dataframe(['AAAA', 'CCCC'], ['AAAA', 'CCCC'], k=0, metric='exact', target_names=['a', 'c'])
    assert frame['target_name'].tolist() == ['a', 'c']


def test_csv_gzip_library_has_correct_delimiter(tmp_path):
    path = tmp_path / 'targets.csv.gz'
    with gzip.open(path, 'wt', encoding='utf-8') as handle:
        handle.write('target_id,sequence,gene\na,acgtacgt,GENE1\nb,tgcatgca,GENE2\n')
    assert dotmatch.load_targets(path) == [('a', 'ACGTACGT'), ('b', 'TGCATGCA')]


def test_missing_fastq_identifier_has_actionable_error(tmp_path):
    path = tmp_path / 'reads.fastq'
    path.write_text('@  \nACGT\n+\nIIII\n')
    with pytest.raises(ValueError, match='missing read identifier'):
        list(dotmatch.iter_fastq(path))


@pytest.mark.parametrize('gz', [False, True])
def test_decompressed_content_hash_preserves_headers_case_and_line_endings(tmp_path, gz):
    text = b'@read1 description\r\nacgt\r\n+\r\nIIII\r\n@read2\nTGCA\n+\n####\n'
    path = tmp_path / ('reads.fastq.gz' if gz else 'reads.fastq')
    path.write_bytes(gzip.compress(text) if gz else text)
    digest = hashlib.sha256()
    records = list(dotmatch.iter_fastq(path, content_digest=digest))
    assert digest.hexdigest() == hashlib.sha256(text).hexdigest()
    assert [record.seq for record in records] == ['ACGT', 'TGCA']


def test_streamed_dataframes_normalize_like_sequence_lists(tmp_path):
    pd = pytest.importorskip('pandas')
    path = tmp_path / 'reads.fastq'
    path.write_text('@r\nACGTACGT\n+\nIIIIIIII\n')
    library = pd.DataFrame({'id': ['a'], 'sequence': ['acgtacgt']})
    records = list(dotmatch.stream_assign(path, library, k=0, metric='exact'))
    assert records[0].status_name == 'unique' and records[0].target_name == 'a'


def test_quickstart_instructions_preserve_reviewed_project(tmp_path, monkeypatch, capsys):
    from argparse import Namespace
    from dotmatch import cli
    library = tmp_path / 'library.tsv'
    library.write_text('target_id\tsequence\na\tACGTACGT\n')
    reads = tmp_path / 'reads.fastq'
    reads.write_text('@r\nACGTACGT\n+\nIIIIIIII\n')
    project = tmp_path / 'project'
    monkeypatch.setattr(cli, 'scaffold_assay_project', lambda **kwargs: {'project': str(project)})
    args = Namespace(library=str(library), fastq=[str(reads)], out=str(project), threads=1, max_reads=100, max_start=20, no_run=False, accept_inference=False)
    assert cli.command_crispr_quickstart(args) == 0
    output = capsys.readouterr().out
    assert 'change status = "draft" to "ready"' in output
    assert 'rerun with --accept-inference' not in output
    assert not (tmp_path / '.project.dotmatch-inputs').exists()
