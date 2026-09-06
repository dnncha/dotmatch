#!/usr/bin/env python3
"""Tag an explicitly approved, checked main commit; verify and publish its artifacts."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request

REPOSITORY = 'dnncha/dotmatch'
ROOT = Path(__file__).resolve().parents[1]
REQUIRED_WORKFLOWS = {
    '.github/workflows/ci.yml',
    '.github/workflows/site-validation.yml',
    '.github/workflows/codeql.yml',
    '.github/workflows/workflow-ecosystem.yml',
    '.github/workflows/pages.yml',
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def command(*args, cwd=ROOT):
    return subprocess.check_output(list(args), cwd=cwd, text=True).strip()


def github(path):
    return json.loads(command('gh', 'api', f'repos/{REPOSITORY}/{path}'))


def project_version():
    match = re.search(r'^version\s*=\s*"(\d+\.\d+\.\d+)"', (ROOT / 'pyproject.toml').read_text(), re.M)
    require(match is not None, 'Missing stable project version')
    return match.group(1)


def validate_request(request, version, repository):
    require(repository == REPOSITORY, 'Publication is limited to the canonical repository')
    require(request.get('authorized') is True, 'Release has not been explicitly authorized')
    require(re.fullmatch(r'\d+\.\d+\.\d+', version) is not None, 'Not a stable semantic version')
    require(request.get('version') == version and request.get('tag') == 'v' + version, 'Request, tag and source version differ')
    return request['tag']


def context():
    version = project_version()
    request = json.loads((ROOT / '.github/release-request.json').read_text())
    tag = validate_request(request, version, os.environ.get('GITHUB_REPOSITORY'))
    sha = command('git', 'rev-parse', 'HEAD')
    require(sha == os.environ.get('GITHUB_SHA'), 'Checkout differs from the workflow commit')
    record = json.loads((ROOT / 'docs/distribution-release.json').read_text())
    require(record.get('publication_authorized') is True and record.get('release_version') == version, 'Distribution authorization differs')
    return version, tag, sha, record


def workflow_verdict(runs, sha):
    latest = {}
    for run in sorted(runs, key=lambda item: item['id'], reverse=True):
        if run.get('head_sha') != sha:
            continue
        path = str(run.get('path', '')).split('@')[0]
        if path in REQUIRED_WORKFLOWS:
            latest.setdefault(path, run)
    for path, run in latest.items():
        if run.get('status') == 'completed' and run.get('conclusion') != 'success':
            raise RuntimeError(f'Required workflow failed: {path}: {run.get("conclusion")}')
    return len(latest) == len(REQUIRED_WORKFLOWS) and all(item.get('conclusion') == 'success' for item in latest.values())


def status_verdict(status):
    """Keep release/security checks authoritative, not an unused preview account.

    GitHub Pages is required separately. The exact Vercel account-level block
    is visible in logs and is not misrepresented as a successful deployment.
    A Vercel build error, or any other failed status, still blocks publication.
    """
    ready = True
    for item in status.get('statuses', []):
        if (item.get('context') == 'Vercel' and item.get('state') == 'failure'
                and item.get('description') == 'Account is blocked.'
                and item.get('target_url') == 'https://vercel.com/knowledge/why-is-my-account-deployment-blocked'):
            print('Non-release preview unavailable: Vercel account is blocked. Production GitHub Pages must pass separately.', flush=True)
            continue
        require(item.get('state') not in {'failure', 'error'}, f'Commit status failed: {item.get("context")}')
        ready = ready and item.get('state') == 'success'
    return ready


def tag_release():
    version, tag, sha, _ = context()
    require(os.environ.get('GITHUB_REF') == 'refs/heads/main', 'Tags may only be requested from main')
    for attempt in range(120):
        require(github('git/ref/heads/main')['object']['sha'] == sha, 'Main moved during release validation; refusing to tag a stale candidate')
        runs = github(f'actions/runs?head_sha={sha}&per_page=100')['workflow_runs']
        external_ready = status_verdict(github(f'commits/{sha}/status'))
        if workflow_verdict(runs, sha) and external_ready:
            break
        print('Required checks for the exact main commit are not all complete.', flush=True)
        time.sleep(10)
    else:
        raise RuntimeError('Required checks did not complete within the bounded release gate')
    existing = command('git', 'ls-remote', 'origin', f'refs/tags/{tag}')
    if existing:
        command('git', 'fetch', 'origin', f'refs/tags/{tag}:refs/tags/{tag}')
        require(command('git', 'cat-file', '-t', f'refs/tags/{tag}') == 'tag', 'Existing tag is not annotated')
        require(command('git', 'rev-parse', f'{tag}^{{commit}}') == sha, 'Existing tag points elsewhere; tags are never moved')
    else:
        command('git', 'config', 'user.name', 'github-actions[bot]')
        command('git', 'config', 'user.email', '41898282+github-actions[bot]@users.noreply.github.com')
        command('git', 'tag', '-a', tag, sha, '-m', f'DotMatch {tag}\nValidated main commit: {sha}')
        command('git', 'push', 'origin', f'refs/tags/{tag}')
    runs = github(f'actions/runs?head_sha={sha}&per_page=100')['workflow_runs']
    if any(run.get('path') == '.github/workflows/release.yml' for run in runs):
        print('The immutable tag already has a release run; use a targeted job retry for any failure.')
        return
    # GITHUB_TOKEN tag pushes do not trigger push workflows. Dispatch explicitly.
    command('gh', 'workflow', 'run', 'release.yml', '--ref', tag, '--repo', REPOSITORY)
    print(f'Created/verified immutable annotated tag {tag} and dispatched the existing trusted-publishing workflow.')


def fetch_json(url, headers=None):
    request = urllib.request.Request(url, headers={'User-Agent': 'DotMatch-release-verifier', **(headers or {})})
    with urllib.request.urlopen(request, timeout=30) as response:
        raw = response.read()
        return json.loads(raw), raw, dict(response.headers)


def checksum_manifest(directory):
    expected = {}
    for line in (directory / 'SHA256SUMS.txt').read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        name = name.lstrip('*')
        require(Path(name).name == name and name not in expected, 'Unsafe or duplicate checksum entry')
        require(re.fullmatch(r'[0-9a-f]{64}', digest) is not None, 'Malformed checksum')
        require((directory / name).is_file(), f'Release artifact missing: {name}')
        require(hashlib.sha256((directory / name).read_bytes()).hexdigest() == digest, f'Artifact hash mismatch: {name}')
        expected[name] = digest
    require(bool(expected), 'No release artifacts to verify')
    return expected


def verify_publication():
    version, tag, sha, record = context()
    require(os.environ.get('GITHUB_REF') == 'refs/tags/' + tag, 'Artifact publication must run on the approved tag')
    ref = github('git/ref/tags/' + tag)['object']
    require(ref['type'] == 'tag', 'Release tag is not annotated')
    require(github('git/tags/' + ref['sha'])['object']['sha'] == sha, 'Release tag and built commit differ')
    verification = ROOT / 'release-verification'
    assets = verification / 'github'
    assets.mkdir(parents=True, exist_ok=True)
    command('gh', 'release', 'download', tag, '--repo', REPOSITORY, '--dir', str(assets), '--pattern', 'SHA256SUMS.txt', '--pattern', '*.whl', '--pattern', '*.tar.gz', '--clobber')
    expected = checksum_manifest(assets)
    pypi_url = f'https://pypi.org/pypi/dotmatch/{version}/json'
    for attempt in range(12):
        try:
            pypi, _, _ = fetch_json(pypi_url)
            break
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
            if attempt == 11:
                raise
            time.sleep(10)
    require(pypi['info']['version'] == version, 'Public PyPI version differs')
    published = {item['filename']: item['digests']['sha256'] for item in pypi['urls']}
    portable = {name: digest for name, digest in expected.items() if not re.search(r'-linux_(x86_64|aarch64)\.whl$', name)}
    require(published == portable, 'PyPI file set or hashes differ from the checked release artifacts')
    require(any('macosx_' in name for name in published), 'macOS wheel missing')
    for architecture in ('x86_64', 'aarch64'):
        for family in ('manylinux', 'musllinux'):
            require(any(family in name and name.endswith(architecture + '.whl') for name in published), f'Missing {family}/{architecture} wheel')
    with tempfile.TemporaryDirectory(prefix='dotmatch-public-install-') as temporary:
        directory = Path(temporary)
        environment = directory / 'venv'
        subprocess.run([sys.executable, '-m', 'venv', str(environment)], check=True)
        python = str(environment / 'bin/python')
        for attempt in range(6):
            result = subprocess.run([python, '-m', 'pip', 'install', '--disable-pip-version-check', '--only-binary=:all:', '--index-url', 'https://pypi.org/simple', 'dotmatch==' + version])
            if result.returncode == 0:
                break
            if attempt == 5:
                raise RuntimeError('Clean installation from public PyPI failed')
            time.sleep(10)
        binary = str(environment / 'bin/dotmatch')
        require(command(binary, '--version', cwd=directory) == 'dotmatch ' + version, 'Installed CLI version differs')
        output = directory / 'sensitivity'
        command(binary, 'sensitivity', '--targets', str(ROOT / 'examples/assignment_sensitivity/targets.tsv'), '--reads', str(ROOT / 'examples/assignment_sensitivity/reads.fastq'), '--target-start', '0', '--target-length', '20', '--out-dir', str(output), cwd=directory)
        summary = json.loads((output / 'summary.json').read_text())
        require(summary['read_count'] == 9 and summary['changed_reads'] == 5, 'Published sensitivity workflow differs from the checked fixture')
        require([summary['outcomes'][mode]['unique'] for mode in ('exact', 'radius_k1', 'best_k1')] == [3, 3, 5], 'Published policy counts differ')
    token, _, _ = fetch_json('https://ghcr.io/token?service=ghcr.io&scope=repository:dnncha/dotmatch:pull')
    manifest, raw, headers = fetch_json('https://ghcr.io/v2/dnncha/dotmatch/manifests/' + tag, {'Authorization': 'Bearer ' + token['token'], 'Accept': 'application/vnd.oci.image.index.v1+json, application/vnd.docker.distribution.manifest.list.v2+json'})
    platforms = {item.get('platform', {}).get('os', '') + '/' + item.get('platform', {}).get('architecture', '') for item in manifest.get('manifests', [])}
    require({'linux/amd64', 'linux/arm64'} <= platforms, 'Published GHCR manifest lacks a required platform')
    digest = 'sha256:' + hashlib.sha256(raw).hexdigest()
    header_digest = next((value for key, value in headers.items() if key.lower() == 'docker-content-digest'), None)
    require(header_digest == digest, 'GHCR manifest digest differs from response bytes')
    date = datetime.now(timezone.utc).date().isoformat()
    run_url = f'https://github.com/{REPOSITORY}/actions/runs/{os.environ["GITHUB_RUN_ID"]}'
    record.update(status='partially_verified', release_commit=sha, tag_object=ref['sha'], github_release_url=f'https://github.com/{REPOSITORY}/releases/tag/{tag}', release_workflow_url=run_url)
    for channel in record['channels']:
        if channel['id'] in {'pypi', 'ghcr'}:
            channel.update(status='verified', verified_date=date, public_url=channel['expected_url'], evidence_url=pypi_url if channel['id'] == 'pypi' else run_url)
            channel.pop('blocker', None); channel.pop('next_action', None)
            if channel['id'] == 'pypi':
                channel.update(artifact_sha256=published, linux_wheel_architectures=['x86_64', 'aarch64'])
            else:
                channel.update(manifest_digest=digest, platforms=['linux/amd64', 'linux/arm64'])
        else:
            channel['status'] = 'blocked'
            channel['blocker'] = f'Version {version} publication on this community/archive channel remains separately unverified.'
    record['blockers'] = [item['blocker'] for item in record['channels'] if item['status'] != 'verified']
    record['next_action'] = 'Verify the matching Bioconda/AssayCode build, generated BioContainers images and newly minted version-specific Zenodo archive independently.'
    (verification / 'distribution-release.json').write_text(json.dumps(record, indent=2) + '\n')
    (verification / 'verification.json').write_text(json.dumps({'version': version, 'commit': sha, 'tag': tag, 'pypi_files': len(published), 'artifact_hashes_verified': True, 'clean_public_install': True, 'sensitivity_fixture': {'reads': 9, 'changed_reads': 5, 'unique_counts': [3, 3, 5]}, 'ghcr_manifest_digest': digest, 'ghcr_platforms': ['linux/amd64', 'linux/arm64'], 'workflow': run_url}, indent=2) + '\n')
    command('gh', 'release', 'upload', tag, str(verification / 'verification.json'), '--repo', REPOSITORY, '--clobber')
    command('gh', 'release', 'edit', tag, '--repo', REPOSITORY, '--draft=false', '--latest', '--title', 'DotMatch ' + version, '--notes-file', str(ROOT / 'docs/releases' / (tag + '.md')))
    print(f'Published {tag}: checked GitHub artifacts, identical public PyPI hashes, clean native workflow installation, and both GHCR architectures.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=['tag', 'finalize'])
    args = parser.parse_args()
    if args.operation == 'tag':
        tag_release()
    else:
        verify_publication()


if __name__ == '__main__':
    main()
