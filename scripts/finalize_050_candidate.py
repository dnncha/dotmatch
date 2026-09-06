#!/usr/bin/env python3
"""One-shot corrections reviewed during the authorized release preflight."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]

def edit(name, old, new):
    path = ROOT / name
    text = path.read_text(encoding='utf-8')
    assert old in text, (name, old)
    path.write_text(text.replace(old, new), encoding='utf-8')

for name in ['README.md', 'docs/sensitivity.md', 'examples/assignment_sensitivity/README.md']:
    edit(name, '/tmp/dotmatch-sensitivity-example', 'sensitivity-example')
edit('scripts/check_repository_ready.py', 'import argparse\n', 'import argparse\nimport hashlib\n')
edit('scripts/check_repository_ready.py', 'RAW_DATA_SUFFIXES =', '''# The exact reviewed, hand-written nine-read policy fixture is allowed.
# Replacing these bytes requires a new provenance review, not a directory exemption.
APPROVED_SYNTHETIC_RAW_SHA256 = {
    "examples/assignment_sensitivity/reads.fastq":
        "cd775aa791e509ec605a8d543713e43726f6fe0b76255df54f3d46dfd0fa4475",
}

RAW_DATA_SUFFIXES =''')
edit('scripts/check_repository_ready.py', '''        if not relative.startswith(ALLOWED_RAW_DATA_PREFIXES):
            offenders.append(relative)
            continue''', '''        expected = APPROVED_SYNTHETIC_RAW_SHA256.get(relative)
        if expected is not None:
            if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                offenders.append(relative)
                continue
        elif not relative.startswith(ALLOWED_RAW_DATA_PREFIXES):
            offenders.append(relative)
            continue''')
path = ROOT / 'docs/agent-crispr.md'
path.write_text(path.read_text() + '''
## Reference contract and scope

This route provides known-guide counting only: no downstream screen statistics.
Inspect the [checked contract fixture](https://dnncha.github.io/dotmatch/agent-reference-crispr.json)
when integrating the tool envelope. Its intentionally failed verdict exercises
unsafe and low-assignment states; it is not biological validation.
''')
edit('README.md', 'The current public release is 0.4.1 and includes the six', 'Release 0.5.0 includes the six')
edit('README.md', 'dotmatch==0.4.1', 'dotmatch==0.5.0')
edit('README.md', 'pinned published container:', 'pinned release container:')
edit('README.md', 'ghcr.io/dnncha/dotmatch:v0.4.1', 'ghcr.io/dnncha/dotmatch:v0.5.0')
edit('README.md', 'The new source-tree `dotmatch sensitivity` command', 'The `dotmatch sensitivity` command, introduced in 0.5.0,')
edit('README.md', '**This command is not in PyPI 0.4.1.** From a reviewed source checkout:', 'Run the included synthetic example from a checkout of the v0.5.0 release:')
edit('README.md', 'python3 -m pip install .\n', 'python3 -m pip install dotmatch==0.5.0\n')
edit('README.md', 'The six structured agent tools are available in published 0.4.1:', 'The six structured agent tools are included in release 0.5.0:')
edit('docs/index.md', 'The current public PyPI release is 0.4.1 and includes the six', 'Release 0.5.0 includes the six')
edit('docs/index.md', 'reports 0.4.1', 'reports 0.5.0')
edit('docs/agent-guide.md', 'This six-tool interface is included in the current public PyPI release, 0.4.1.', 'This six-tool interface is included in release 0.5.0.')
edit('docs/agent-guide.md', 'until that channel reaches 0.4.1', 'until that channel reaches 0.5.0')
edit('scripts/check_agent_discovery.py', '    readme = _read(root, "README.md")\n', '    release = re.search(r\'^version\\s*=\\s*"([^\"]+)"\', pyproject, re.M).group(1)\n    readme = _read(root, "README.md")\n')
edit('scripts/check_agent_discovery.py', '"The current public release is 0.4.1"', 'f"Release {release}"')
edit('scripts/check_agent_discovery.py', '"The current public PyPI release is 0.4.1"', 'f"Release {release}"')
edit('scripts/check_agent_discovery.py', '"This six-tool interface is included in the current public PyPI release, 0.4.1."', 'f"This six-tool interface is included in release {release}."')
edit('scripts/check_agent_discovery.py', '"until that channel reaches 0.4.1"', 'f"until that channel reaches {release}"')
edit('docs/sensitivity.md', 'The source-tree `dotmatch sensitivity` command', 'The `dotmatch sensitivity` command, introduced in 0.5.0,')
edit('docs/sensitivity.md', '''It is new source functionality for the next release, not part of published 0.4.1.
Install from a reviewed source checkout with `python3 -m pip install .`.''', 'Install the matching release with `python3 -m pip install dotmatch==0.5.0`.')
edit('docs/sensitivity.md', 'python3 -m pip install .\n', 'python3 -m pip install dotmatch==0.5.0\n')
edit('examples/assignment_sensitivity/README.md', 'From a reviewed source checkout (the sensitivity command is not in PyPI 0.4.1):', 'From a checkout of the v0.5.0 release:')
edit('examples/assignment_sensitivity/README.md', 'python3 -m pip install .\n', 'python3 -m pip install dotmatch==0.5.0\n')
print('Portable examples, exact synthetic fixture approval, release documentation and agent boundaries aligned.')
