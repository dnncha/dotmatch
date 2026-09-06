#!/usr/bin/env python3
"""Test an already-built image's public workflows. Does not push an image."""
from __future__ import annotations
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    image = sys.argv[1] if len(sys.argv) == 2 else 'dotmatch:integrity'
    with tempfile.TemporaryDirectory() as temporary:
        work = Path(temporary)
        for name in ('targets.tsv', 'reads.fastq'):
            (work / name).write_bytes((ROOT / 'examples/assignment_sensitivity' / name).read_bytes())
        docker = ['docker', 'run', '--rm', '--network=none', '--read-only',
                  '--user', f'{os.getuid()}:{os.getgid()}',
                  '--tmpfs', '/tmp:rw,nosuid,nodev,size=64m', '-v', f'{work}:/data']
        def run(*args, entrypoint=None):
            command = [*docker]
            if entrypoint: command += ['--entrypoint', entrypoint]
            command += [image, *args]
            result = subprocess.run(command, text=True, capture_output=True, timeout=120)
            if result.returncode:
                raise RuntimeError(f'{args!r} failed:\n{result.stdout}\n{result.stderr}')
            return result.stdout
        assert run('dist', 'ACGT', 'AGGT').strip() == '1'
        assert run('leq', '1', 'ACGT', 'AGGT').strip() == 'true'
        version = run('--version').strip()
        assert version == run('-m', 'dotmatch', '--version', entrypoint='python').strip()
        assert json.loads(run('agent', 'tools', '--json'))
        run('assay', '--help')
        run('sensitivity', '--targets', '/data/targets.tsv', '--reads', '/data/reads.fastq',
            '--target-start', '0', '--target-length', '20', '--out-dir', '/data/sensitivity')
        summary = json.loads((work / 'sensitivity/summary.json').read_text())
        assert summary['read_count'] == 9
        for mode in ('exact', 'radius_k1', 'best_k1'):
            assert sum(summary['outcomes'][mode].values()) == 9
        run('count', '--targets', '/data/targets.tsv', '--reads', '/data/reads.fastq',
            '--sample-label', 'sample', '--target-start', '0', '--target-length', '20',
            '--k', '0', '--metric', 'hamming', '--format', 'mageck', '--out', '/data/counts.tsv')
        assert (work / 'counts.tsv').is_file()
        run('-c', 'import shutil,pathlib; assert not shutil.which("cc"); assert not shutil.which("git"); assert not pathlib.Path("/build").exists(); assert not pathlib.Path("/opt/dotmatch/.git").exists()', entrypoint='python')
        print(json.dumps({'version': version, 'native_commands': 'passed', 'agent_tools': 'passed',
            'sensitivity': 'passed', 'counts': 'passed', 'read_only_nonroot_offline': 'passed',
            'runtime_compiler_and_git': 'absent'}, indent=2))


if __name__ == '__main__':
    main()
