#!/usr/bin/env python3
"""Prepare the explicitly authorized 0.5.0 candidate; publication is separate."""
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OLD, NEW = '0.4.1', '0.5.0'
assert subprocess.check_output(['git', 'rev-parse', 'HEAD^'], cwd=ROOT, text=True).strip() == '310c3e0fff133be31b88dd9c307bb36ea2f25fbc', 'Unexpected preparation parent'

def read(name): return (ROOT / name).read_text(encoding='utf-8')
def write(name, text):
    path = ROOT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding='utf-8')
def scalar(name, old, new):
    text = read(name)
    assert old in text, (name, old)
    write(name, text.replace(old, new))
def json_file(name, change):
    data = json.loads(read(name)); change(data)
    write(name, json.dumps(data, indent=2, ensure_ascii=False) + '\n')

for name in ['pyproject.toml', 'package.json', 'package-lock.json', 'codemeta.json', '.zenodo.json', 'CITATION.cff', 'Dockerfile', 'include/qdalign.h', 'packaging/bioconda/meta.yaml', 'packaging/bioconda/assaycode-meta.yaml', 'DESCRIPTION', 'python/dotmatch/__init__.py', 'docs/conf.py', 'docs/registries/biotools.yml', 'docs/methods-and-citation.md']:
    scalar(name, OLD, NEW)
# Do not attribute the prior 0.4.0 archive's DOI to this new release.
write('CITATION.cff', re.sub(r'^doi:.*\n', '', read('CITATION.cff'), flags=re.M))
json_file('.zenodo.json', lambda data: [item.update(identifier='10.5281/zenodo.22214073') for item in data.get('related_identifiers', []) if item.get('relation') == 'isNewVersionOf'])
for name in ['agent-capabilities.json', 'agent-tools.json']:
    json_file(name, lambda data: data.update(generated_for_version=NEW))
json_file('agent-capabilities.json', lambda data: data['install'].update(recommended=f'python3 -m pip install dotmatch=={NEW}'))
json_file('agent-reference-crispr.json', lambda data: data.update(dotmatch_version=NEW))
subprocess.run(['python3', 'scripts/sync_agent_discovery.py'], cwd=ROOT, check=True)

text = read('README.md')
text = text.replace('installed version. Review the [packaging details]', 'installed version. The Bioconda recipe includes `osx-arm64` for Apple Silicon.\nReview the [packaging details]')
assert 'osx-arm64' in text
text = text.replace('Use the [methods and citation guide]', 'Use `dotmatch citation` and [CITATION.cff](https://github.com/dnncha/dotmatch/blob/main/CITATION.cff)\nto record the actual software version. Use the [methods and citation guide]')
write('README.md', text)
text = read('CHANGELOG.md').replace('## Unreleased', '## Unreleased\n\n## 0.5.0 - 2026-09-06', 1)
text = text.replace('### Fixed\n', '''### Fixed

- Defer optional scientific imports so basic matching and CLI use do not load
  dataframe stacks unnecessarily; keep explicit dataframe integrations working.
- Do not expose an ambiguous nearest candidate as an assigned dataframe label.
- Normalize dataframe target case, diagnose missing FASTQ identifiers, and
  retain original decompressed bytes when a caller requests a content digest.
- Correct draft-project review instructions and preserve actionable input errors.
''', 1)
text = text.replace('The new sensitivity command requires a source install until a release is tagged.\n', '')
write('CHANGELOG.md', text)
# Archive previous verified facts instead of re-labelling them as 0.5.0 evidence.
write('docs/releases/distribution-0.4.1.json', read('docs/distribution-release.json'))
channels = []
expected = {
 'pypi': 'https://pypi.org/project/dotmatch/0.5.0/',
 'ghcr': 'https://github.com/dnncha/dotmatch/pkgs/container/dotmatch',
 'bioconda': 'https://anaconda.org/bioconda/dotmatch',
 'bioconda-assaycode': 'https://anaconda.org/bioconda/assaycode',
 'biocontainers': 'https://quay.io/repository/biocontainers/dotmatch',
 'zenodo': 'https://doi.org/10.5281/zenodo.20541628',
}
for name, url in expected.items():
    blocker = 'Version 0.5.0 has not yet been published and independently verified on this channel.'
    action = 'Publish through the approved release workflow, then verify the public version and artifacts.'
    if name == 'pypi': action = 'Publish the source distribution, macOS wheel, and repaired Linux wheels through trusted publishing, then verify their hashes and a clean installation.'
    if name.startswith('bioconda') or name == 'biocontainers': action = 'Submit the immutable release recipe upstream and verify community build propagation; do not infer availability from the source tag.'
    if name == 'zenodo': action = 'Verify the archive generated for the immutable GitHub release; add only the newly minted version DOI.'
    channels.append(dict(id=name, status='prepared', expected_url=url, verification_command='make distribution-channels', blocker=blocker, next_action=action))
manifest = dict(schema_version=1, status='not_released', release_version=NEW, publication_authorized=True, release_tag='v'+NEW, post_release_gate='make distribution-channels', channels=channels, blockers=['Version 0.5.0 is an authorized release candidate; artifact publication and public-channel checks are pending.'], next_action='Pass the release gates, publish the immutable v0.5.0 tag through the existing release workflow, then record each verified channel separately.')
write('docs/distribution-release.json', json.dumps(manifest, indent=2) + '\n')
print('Prepared aligned 0.5.0 metadata; prior public-channel evidence archived; no publication claimed.')
