"""Public-boundary, posterior, FASTQ and matcher-lifetime regressions."""
import gzip
import hashlib
import io
import threading
from concurrent.futures import ThreadPoolExecutor

import pytest
import dotmatch
from dotmatch import core
from dotmatch.fastq_io import iter_fastq_records


@pytest.mark.parametrize('bad', [True, False, 0.5, 1.9, '1', -1, 2**31, 2**32, 2**32+1])
def test_thresholds_never_round_or_wrap(bad):
    with dotmatch.Matcher(['AAAA']) as matcher:
        for call in (
            lambda: dotmatch.distance_leq('AAAA', 'AAAT', bad),
            lambda: dotmatch.assign(['AAAT'], ['AAAA'], k=bad),
            lambda: matcher.assign(['AAAT'], k=bad),
            lambda: matcher.assign_hamming(['AAAT'], k=bad),
            lambda: matcher.assign_status_with_stats(['AAAT'], k=bad),
        ):
            with pytest.raises((ValueError, TypeError)):
                call()


def test_numpy_integer_threshold_is_supported_but_bool_is_not():
    np = pytest.importorskip('numpy')
    assert dotmatch.distance_leq('AAAA', 'AAAT', np.int64(1))
    with pytest.raises(TypeError):
        dotmatch.distance_leq('AAAA', 'AAAT', np.bool_(True))


@pytest.mark.parametrize('priors', [[float('nan'), 1], [float('inf'), 1], [0, 0], [-1, 1]])
def test_invalid_priors_fail_explicitly(priors):
    with pytest.raises(ValueError, match='priors'):
        dotmatch.assign_posterior('AAAA', ['AAAA', 'AAAT'], 'IIII', priors=priors)


def test_large_finite_priors_and_ties_have_well_defined_results():
    result = dotmatch.assign_posterior('AAAA', ['AAAA', 'AAAT'], 'IIII', priors=[1e308, 1e308])
    assert result.status == dotmatch.MATCH_UNIQUE
    assert sum(result.posteriors) == pytest.approx(1)
    for cutoff in [0, .1, .5, 1]:
        tied = dotmatch.assign_posterior('AAAA', ['AAAA', 'AAAA'], 'IIII', min_posterior=cutoff)
        assert tied.status == dotmatch.MATCH_AMBIGUOUS


@pytest.mark.parametrize('qual', [b'III\x7f', b'III\xff', b'III '])
def test_posterior_rejects_non_phred33_bytes(qual):
    with pytest.raises(ValueError, match='Phred'):
        dotmatch.assign_posterior('AAAA', ['AAAA'], qual)


def test_single_sequence_needs_explicit_batch():
    with pytest.raises(TypeError, match='wrap'):
        dotmatch.assign('ACGT', ['ACGT'])
    with pytest.raises(TypeError, match='wrap'):
        dotmatch.Matcher('ACGT')


def test_close_waits_for_native_call_before_freeing_index(monkeypatch):
    entered, release, closing, closed = (threading.Event() for _ in range(4))
    original = core._LIB.qdaln_index_assign_stats
    def paused(*args):
        entered.set()
        assert release.wait(5)
        return original(*args)
    monkeypatch.setattr(core._LIB, 'qdaln_index_assign_stats', paused)
    matcher = dotmatch.Matcher(['ACGT'])
    def close():
        closing.set()
        matcher.close()
        closed.set()
    with ThreadPoolExecutor(2) as pool:
        query = pool.submit(matcher.assign, ['ACGT'], 0)
        assert entered.wait(5)
        closing_job = pool.submit(close)
        assert closing.wait(5)
        try:
            assert not closed.wait(.05), 'index was freed while ctypes could still be using it'
        finally:
            release.set()
        assert query.result(timeout=5)[0].status == dotmatch.MATCH_UNIQUE
        closing_job.result(timeout=5)
    matcher.close()  # idempotent
    with pytest.raises(ValueError, match='closed'):
        matcher.assign_exact(['ACGT'])


@pytest.mark.parametrize('text', ['@\nACGT\n+\nIIII\n', '@r\nACGT\n+\nIII\n', '@r\nAC T\n+\nIIII\n', '@r\nACGT\n+\nIII \n', '@r\nAC\x00T\n+\nIIII\n', '@r\nACGT\n+x\nIIII\n'])
def test_python_workflows_share_fastq_failures(tmp_path, text):
    from dotmatch import cli, sensitivity
    path = tmp_path / 'bad.fastq'
    path.write_text(text)
    for read in (lambda: list(dotmatch.iter_fastq(path)), lambda: list(cli._iter_fastq(path)),
                 lambda: list(sensitivity._fastq(io.StringIO(text), path))):
        with pytest.raises(ValueError, match='record 1'):
            read()


def test_shared_fastq_digest_preserves_original_case_and_line_endings(tmp_path):
    original = b'@r description\r\nacgn\r\n+r\r\nI###\r\n'
    path = tmp_path / 'reads.FASTQ.GZ'
    path.write_bytes(gzip.compress(original))
    digest = hashlib.sha256()
    assert list(dotmatch.iter_fastq(path, content_digest=digest))[0].seq == 'ACGN'
    assert digest.hexdigest() == hashlib.sha256(original).hexdigest()


def test_empty_read_file_cannot_hide_invalid_configuration(tmp_path):
    path = tmp_path / 'empty.fastq'
    path.write_text('')
    for settings in ({'k': -1}, {'k': 1.5}, {'batch_size': 1.5}, {'policy': 'invented'}, {'target_start': True}):
        with pytest.raises((ValueError, TypeError)):
            list(dotmatch.stream_assign(path, [('g', 'ACGT')], **settings))


def test_unrelated_pyproject_does_not_override_package_version(tmp_path, monkeypatch):
    package = tmp_path / 'vendor' / 'dotmatch'
    package.mkdir(parents=True)
    (tmp_path / 'pyproject.toml').write_text('[project]\nname = "unrelated"\nversion = "99.0.0"\n')
    monkeypatch.setattr(dotmatch, '__file__', str(package / '__init__.py'))
    assert dotmatch._source_tree_version() is None


def test_native_lookup_never_searches_the_working_directory(monkeypatch, tmp_path):
    from dotmatch import core, native
    for key in ('DOTMATCH_LIB', 'QUICKDNA_LIB', 'DOTMATCH_NATIVE_CLI'):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.chdir(tmp_path)
    assert all(path.parent != tmp_path for path in core._candidate_paths())
    assert all(path.parent != tmp_path for path in native.native_cli_candidates())
    monkeypatch.setenv('DOTMATCH_NATIVE_CLI', str(tmp_path / 'missing'))
    with pytest.raises(FileNotFoundError, match='missing'):
        native.find_native_cli()
    monkeypatch.setenv('DOTMATCH_LIB', str(tmp_path / 'missing.so'))
    with pytest.raises(RuntimeError, match='missing.so'):
        core._load_lib()


def test_installed_native_lookup_does_not_walk_unrelated_ancestors(monkeypatch, tmp_path):
    from dotmatch import core, native
    for key in ('DOTMATCH_LIB', 'QUICKDNA_LIB', 'DOTMATCH_NATIVE_CLI'):
        monkeypatch.delenv(key, raising=False)
    installed = tmp_path / 'site-packages' / 'dotmatch'
    monkeypatch.setattr(core, '__file__', str(installed / 'core.py'))
    monkeypatch.setattr(native, '__file__', str(installed / 'native.py'))
    assert all(path.parent == installed for path in core._candidate_paths())
    assert native.native_cli_candidates() == [installed / 'dotmatch-native']


def test_python_module_uses_the_same_cli():
    import os
    import subprocess
    import sys
    from pathlib import Path
    root = Path(__file__).resolve().parents[2]
    env = {**os.environ, 'PYTHONPATH': str(root / 'python')}
    version = subprocess.run([sys.executable, '-m', 'dotmatch', '--version'],
        text=True, capture_output=True, check=True, env=env, timeout=10)
    assert version.stdout.strip() == f'dotmatch {dotmatch.__version__}'
    distance = subprocess.run([sys.executable, '-m', 'dotmatch', 'dist', 'ACGT', 'AGGT'],
        text=True, capture_output=True, check=True, env=env, timeout=10)
    assert distance.stdout.strip() == '1'
