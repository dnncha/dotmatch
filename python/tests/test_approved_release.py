"""Release guards are tested without publishing or contacting external services."""
import hashlib
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('approved_release', ROOT / 'scripts/approved_release.py')
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)


def request(**changes):
    return dict({'version': '0.5.0', 'tag': 'v0.5.0', 'authorized': True}, **changes)


def test_authorized_request_is_limited_to_canonical_stable_release():
    assert release.validate_request(request(), '0.5.0', 'dnncha/dotmatch') == 'v0.5.0'
    for value, version, repository in [
        (request(authorized=False), '0.5.0', 'dnncha/dotmatch'),
        (request(authorized=1), '0.5.0', 'dnncha/dotmatch'),
        (request(), '0.5.1', 'dnncha/dotmatch'),
        (request(tag='v0.4.1'), '0.5.0', 'dnncha/dotmatch'),
        (request(), '0.5.0', 'another/dotmatch'),
        (request(version='0.5.0rc1', tag='v0.5.0rc1'), '0.5.0rc1', 'dnncha/dotmatch'),
    ]:
        with pytest.raises(RuntimeError):
            release.validate_request(value, version, repository)


def completed(sha='abc'):
    return [dict(id=index, path=path, head_sha=sha, status='completed', conclusion='success') for index, path in enumerate(sorted(release.REQUIRED_WORKFLOWS), 1)]


def test_every_required_check_must_pass_for_the_exact_commit():
    assert release.workflow_verdict(completed(), 'abc')
    assert not release.workflow_verdict(completed('other'), 'abc')
    assert not release.workflow_verdict(completed()[:-1], 'abc')
    rows = completed()
    rows[0].update(status='in_progress', conclusion=None)
    assert not release.workflow_verdict(rows, 'abc')


@pytest.mark.parametrize('conclusion', ['failure', 'cancelled', 'timed_out', 'skipped'])
def test_failed_or_skipped_checks_cannot_release(conclusion):
    rows = completed()
    rows[0]['conclusion'] = conclusion
    with pytest.raises(RuntimeError, match='Required workflow failed'):
        release.workflow_verdict(rows, 'abc')


def test_latest_run_wins_and_unrelated_runs_do_not_satisfy_the_gate():
    rows = completed()
    old = dict(rows[0], id=0, conclusion='failure')
    assert release.workflow_verdict(rows + [old], 'abc')
    newer = dict(rows[0], id=100, status='in_progress', conclusion=None)
    assert not release.workflow_verdict(rows + [newer], 'abc')
    assert not release.workflow_verdict([dict(id=1, head_sha='abc', path='unrelated.yml', status='completed', conclusion='success')], 'abc')


def test_release_manifest_checks_every_listed_file(tmp_path):
    data = b'release artifact'
    path = tmp_path / 'dotmatch-0.5.0.tar.gz'
    path.write_bytes(data)
    digest = hashlib.sha256(data).hexdigest()
    (tmp_path / 'SHA256SUMS.txt').write_text(digest + '  ' + path.name + '\n')
    assert release.checksum_manifest(tmp_path) == {path.name: digest}
    path.write_bytes(b'changed bytes')
    with pytest.raises(RuntimeError, match='hash mismatch'):
        release.checksum_manifest(tmp_path)


@pytest.mark.parametrize('entry', ['../escape.whl', '/absolute.whl', 'nested/file.whl'])
def test_release_manifest_rejects_paths_outside_asset_directory(tmp_path, entry):
    (tmp_path / 'SHA256SUMS.txt').write_text('a' * 64 + '  ' + entry + '\n')
    with pytest.raises(RuntimeError, match='Unsafe'):
        release.checksum_manifest(tmp_path)


def test_release_manifest_rejects_duplicate_missing_or_empty_artifacts(tmp_path):
    manifest = tmp_path / 'SHA256SUMS.txt'
    manifest.write_text('')
    with pytest.raises(RuntimeError, match='No release artifacts'):
        release.checksum_manifest(tmp_path)
    manifest.write_text('a' * 64 + '  missing.whl\n')
    with pytest.raises(RuntimeError, match='missing'):
        release.checksum_manifest(tmp_path)
    file = tmp_path / 'artifact.whl'
    file.write_bytes(b'x')
    line = hashlib.sha256(b'x').hexdigest() + '  artifact.whl\n'
    manifest.write_text(line + line)
    with pytest.raises(RuntimeError, match='duplicate'):
        release.checksum_manifest(tmp_path)
