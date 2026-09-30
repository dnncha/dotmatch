#!/usr/bin/env python3
"""Verify the distributed AR-001 bundle and reproduce its summary offline.

Usage: python verify_bundle.py --bundle . --out regenerated-summary
Python 3.11+; standard library only. Output must not already exist.
The original raw-read experiments are not rerun by this command.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import stat
import sys
import tempfile
from types import SimpleNamespace
import zipfile


def require(value, message):
    if not value:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def safe_extract(archive, destination):
    destination.mkdir()
    total = 0
    with zipfile.ZipFile(archive) as handle:
        names = set()
        for item in handle.infolist():
            require(item.filename not in names, 'duplicate archive member')
            names.add(item.filename)
            require(not stat.S_ISLNK(item.external_attr >> 16), 'archive symlink rejected')
            path = (destination / item.filename).resolve()
            require(path.is_relative_to(destination.resolve()), 'unsafe archive member path')
            total += item.file_size
            require(total <= 2_000_000_000, 'archive exceeds uncompressed size limit')
        handle.extractall(destination)


def run(bundle, output):
    require(not output.exists(), 'output directory must not exist')
    manifest = json.loads((bundle / 'bundle-manifest.json').read_text())
    require(manifest['files'], 'empty bundle manifest')
    for relative, expected in manifest['files'].items():
        path = (bundle / relative).resolve()
        require(path.is_relative_to(bundle), 'unsafe bundle path')
        require(path.is_file() and path.stat().st_size == expected['bytes'], f'missing or wrong-size file: {relative}')
        require(digest(path) == expected['sha256'], f'file hash mismatch: {relative}')
    ledger = json.loads((bundle / 'archive-ledger.json').read_text())
    roots = {}
    with tempfile.TemporaryDirectory(prefix='dotmatch-ar001-') as tmp:
        work = Path(tmp)
        for item in ledger['archives']:
            archive = bundle / 'evidence' / item['filename']
            require(archive.stat().st_size == item['bytes'] and digest(archive) == item['sha256'], 'archive digest mismatch')
            if item['status'] != 'complete':
                continue
            destination = work / item['lane']
            safe_extract(archive, destination)
            roots[item['lane']] = str(destination / item['output_directory'])
        require(set(roots) == {'pilot', 'yusa_prefix', 'brunello_prefix', 'overlap', 'full_yusa'}, 'successful evidence cohort mismatch')
        sys.path.insert(0, str(bundle / 'source/research/assignment_robustness'))
        import summarize
        summarize.summarize(SimpleNamespace(**roots, out=str(output)))
    expected = json.loads((bundle / 'summary/results.json').read_text())
    actual = json.loads((output / 'results.json').read_text())
    require(actual == expected, 'recomputed descriptive summary differs from distributed results')
    print(f'PASS: {len(manifest["files"])} bundle files verified; all descriptive results reproduced exactly.')
    print('Raw FASTQ assignment is a separate reproduction step. No biological hit calling was performed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', default='.')
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    run(Path(args.bundle).resolve(), Path(args.out).resolve())
