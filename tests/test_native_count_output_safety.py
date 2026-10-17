#!/usr/bin/env python3
"""Count reruns must preserve inputs and previous artifacts on failure."""
import os
import resource
import signal
import shlex
import sys
from pathlib import Path
import subprocess
import tempfile

BIN = Path(os.environ.get('DOTMATCH_BIN', Path(__file__).resolve().parents[1] / 'dotmatch')).resolve()
VALID = b'@ok\nACGTACGT\n+\nIIIIIIII\n'


def main():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        library = root / 'targets.tsv'
        library.write_bytes(b'g\tACGTACGT\n')
        reads = root / 'reads.fastq'
        reads.write_bytes(VALID)
        roles = ['--out', '--assignments', '--ambiguous-out', '--unmatched-out', '--sample-qc', '--target-counts-long', '--summary', '--report']
        paths = [root / f'output-{i}.txt' for i in range(len(roles))]
        base = [str(BIN), 'count', '--targets', str(library), '--reads', str(reads), '--target-start', '0', '--target-length', '8', '--k', '0', '--metric', 'hamming', '--threads', '1']

        def run(outputs):
            return subprocess.run(base + [part for role, path in zip(roles, outputs) for part in (role, str(path))], capture_output=True)

        before = {path: b'previous bytes\n' for path in paths}
        for path, payload in before.items():
            path.write_bytes(payload)
        reads.write_bytes(VALID + b'@broken\nACGTACGT\n+\n')
        result = run(paths)
        assert result.returncode != 0, result.stderr
        for path, payload in before.items():
            assert path.exists() and path.read_bytes() == payload, f'failed rerun destroyed {path.name}'
        fresh = [root / f'fresh-{i}' for i in range(len(paths))]
        assert run(fresh).returncode != 0
        assert not any(path.exists() for path in fresh)
        reads.write_bytes(VALID)
        bad_parent = paths.copy()
        bad_parent[-1] = root / 'missing' / 'report.html'
        assert run(bad_parent).returncode != 0
        for path, payload in before.items():
            assert path.read_bytes() == payload, f'secondary failure changed {path.name}'
        for source in (library, reads):
            for alias_kind in ('same', 'hardlink', 'symlink'):
                alias = root / 'alias'
                alias.unlink(missing_ok=True)
                if alias_kind == 'hardlink':
                    os.link(source, alias)
                elif alias_kind == 'symlink':
                    alias.symlink_to(source)
                else:
                    alias = source
                outputs = paths.copy()
                outputs[0] = alias
                original = source.read_bytes()
                assert run(outputs).returncode != 0, alias_kind
                assert source.read_bytes() == original, f'{alias_kind} destroyed input'
        unrelated = root / 'unrelated.txt'
        unrelated.write_bytes(b'unrelated bytes\n')
        symlink = root / 'symlink-output'
        symlink.symlink_to(unrelated)
        outputs = paths.copy()
        outputs[3] = symlink
        assert run(outputs).returncode != 0
        assert unrelated.read_bytes() == b'unrelated bytes\n'
        duplicate = paths.copy()
        duplicate[1] = duplicate[0]
        assert run(duplicate).returncode != 0
        for index in range(len(paths)):
            outputs = paths.copy()
            outputs[index] = reads
            assert run(outputs).returncode != 0, roles[index]
            assert reads.read_bytes() == VALID
        linked = root / 'linked-output'
        os.link(paths[0], linked)
        outputs = paths.copy()
        outputs[1] = linked
        assert run(outputs).returncode != 0
        linked.unlink()
        outputs[1] = root / '.' / paths[0].name
        assert run(outputs).returncode != 0
        def limit_writes():
            signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
            resource.setrlimit(resource.RLIMIT_FSIZE, (64, 64))
        command = base + [part for role, path in zip(roles, paths) for part in (role, str(path))]
        assert subprocess.run(command, capture_output=True, preexec_fn=limit_writes).returncode != 0
        for path, payload in before.items():
            assert path.read_bytes() == payload, f'write failure changed {path.name}'
        if sys.platform.startswith('linux'):
            shim = root / 'rename_failure.so'
            source = Path(__file__).with_name('count_rename_failure.c')
            subprocess.run(shlex.split(os.environ.get('CC', 'cc')) + ['-shared', '-fPIC', str(source), '-ldl', '-o', str(shim)], check=True)
            for mode in ('edit', 'replace'):
                environment = dict(os.environ, LD_PRELOAD=str(shim), DOTMATCH_MUTATE_TRIGGER='close',
                                   DOTMATCH_MUTATE_PATH=str(paths[3]), DOTMATCH_MUTATE_MODE=mode)
                assert subprocess.run(command, capture_output=True, env=environment).returncode != 0
                for path, payload in before.items():
                    expected = b'external bytes\n' if path == paths[3] else payload
                    assert path.read_bytes() == expected, f'{mode} race overwrote {path.name}'
                paths[3].write_bytes(before[paths[3]])
            fresh_command = base + [part for role, path in zip(roles, fresh) for part in (role, str(path))]
            environment = dict(os.environ, LD_PRELOAD=str(shim), DOTMATCH_MUTATE_TRIGGER='close',
                               DOTMATCH_MUTATE_PATH=str(fresh[3]), DOTMATCH_MUTATE_MODE='create')
            assert subprocess.run(fresh_command, capture_output=True, env=environment).returncode != 0
            assert fresh[3].read_bytes() == b'external bytes\n'
            assert not any(path.exists() for path in fresh if path != fresh[3])
            fresh[3].unlink()
            environment = dict(os.environ, LD_PRELOAD=str(shim), DOTMATCH_FAIL_PUBLICATION_PATH=str(paths[3]))
            assert subprocess.run(command, capture_output=True, env=environment).returncode != 0
            for path, payload in before.items():
                assert path.read_bytes() == payload, f'publication failure changed {path.name}'
            fresh_command = base + [part for role, path in zip(roles, fresh) for part in (role, str(path))]
            environment['DOTMATCH_FAIL_PUBLICATION_PATH'] = str(fresh[3])
            assert subprocess.run(fresh_command, capture_output=True, env=environment).returncode != 0
            assert not any(path.exists() for path in fresh)
            environment = dict(os.environ, LD_PRELOAD=str(shim), DOTMATCH_FAIL_PUBLICATION_PATH=str(paths[3]),
                               DOTMATCH_MUTATE_TRIGGER='failure', DOTMATCH_MUTATE_PATH=str(paths[0]),
                               DOTMATCH_MUTATE_MODE='replace')
            result = subprocess.run(command, capture_output=True, env=environment)
            assert result.returncode != 0
            assert paths[0].read_bytes() == b'external bytes\n'
            for path, payload in before.items():
                if path != paths[0]:
                    assert path.read_bytes() == payload
            backups = list(root.glob('*.dotmatch-backup-*'))
            assert len(backups) == 1, result.stderr
            assert backups[0].read_bytes() == before[paths[0]]
            assert str(backups[0]).encode() in result.stderr
            backups[0].unlink()
            paths[0].write_bytes(before[paths[0]])
        assert run(paths).returncode == 0
        assert paths[0].read_bytes() != before[paths[0]]
        assert run(paths).returncode == 0
        long_name = root / ('x' * min(255, os.pathconf(root, 'PC_NAME_MAX')))
        long_name.write_bytes(b'old long filename bytes\n')
        outputs = paths.copy()
        outputs[0] = long_name
        assert run(outputs).returncode == 0
        assert long_name.read_bytes() != b'old long filename bytes\n'
        samples = root / 'samples.tsv'
        samples.write_text(f'sample\tfastq\nsample\t{reads}\n')
        original = samples.read_bytes()
        command = [str(BIN), 'crispr-count', '--library', str(library), '--samples', str(samples),
                   '--guide-start', '0', '--guide-length', '8', '--k', '0', '--out', str(samples)]
        assert subprocess.run(command, capture_output=True).returncode != 0
        assert samples.read_bytes() == original
        assert not (root / 'sample_qc.tsv').exists()
        guide_library = root / 'guides.tsv'
        guide_library.write_text('guide\tsequence\tgene\ng\tACGTACGT\tGENE\n')
        prefix = root / 'compat'
        auxiliary = [root / f'compat.{suffix}' for suffix in ('counts.txt', 'extended-counts.txt', 'stats.txt')]
        for path in auxiliary:
            path.write_bytes(b'previous compatibility bytes\n')
        reads.write_bytes(VALID + b'@broken\nACGTACGT\n+\n')
        command = [str(BIN), 'guide-counter', 'count', '--library', str(guide_library),
                   '--input', str(reads), '--samples', 'sample', '--output', str(prefix)]
        assert subprocess.run(command, capture_output=True).returncode != 0
        for path in auxiliary:
            assert path.read_bytes() == b'previous compatibility bytes\n'
        reads.write_bytes(VALID)
        assert subprocess.run(command, capture_output=True).returncode == 0
        annotation = auxiliary[2]
        original_annotation = annotation.read_bytes()
        assert subprocess.run(command + ['--essential-genes', str(annotation)], capture_output=True).returncode != 0
        assert annotation.read_bytes() == original_annotation
        assert not list(root.glob('*.dotmatch-*')), 'temporary output leak'
    print('native count output safety passed')


if __name__ == '__main__':
    main()
