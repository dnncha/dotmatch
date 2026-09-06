#!/usr/bin/env python3
"""Verified public Yusa acquisition and pinned DotMatch runs. No synthetic fallback."""
from __future__ import annotations
import csv, gzip, hashlib, io, json, os, platform, subprocess, sys, time, traceback
import urllib.parse, urllib.request, zipfile
from pathlib import Path

ROOT = Path('research-output')
RUNS = ('ERR376998', 'ERR376999')
LIB_URLS = (
 'https://sourceforge.net/projects/mageck/files/libraries/yusa_library.csv.zip/download',
 'https://downloads.sourceforge.net/project/mageck/libraries/yusa_library.csv.zip',
)

def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()

def save(path, obj):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.pending')
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + '\n')
    tmp.replace(path)

def request(url):
    for attempt in range(3):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'DotMatch-research/0.2'}), timeout=90)
        except Exception:
            if attempt == 2: raise
            time.sleep(2 ** attempt)

def fastq(path):
    with gzip.open(path, 'rb') as f:
        ordinal = 0
        while True:
            h = f.readline()
            if not h: return
            s, p, q = f.readline(), f.readline(), f.readline()
            ordinal += 1
            if not h.startswith(b'@') or not p.startswith(b'+') or not q:
                raise ValueError(f'Malformed FASTQ record {ordinal}')
            seq, qual = s.rstrip(b'\r\n'), q.rstrip(b'\r\n')
            if len(seq) != len(qual) or not seq or any(c < 33 or c > 126 for c in qual):
                raise ValueError(f'Invalid sequence/quality record {ordinal}')
            yield ordinal, (h, s, p, q)

def main():
    ROOT.mkdir(exist_ok=False)
    raw = ROOT / 'raw'; raw.mkdir()
    meta = []
    for accession in RUNS:
        url = 'https://www.ebi.ac.uk/ena/portal/api/filereport?' + urllib.parse.urlencode({
          'accession': accession, 'result': 'read_run', 'format': 'json',
          'fields': 'run_accession,study_accession,sample_accession,sample_title,fastq_ftp,fastq_md5,fastq_bytes,read_count'})
        with request(url) as f: data = f.read(2000001)
        if len(data) > 2000000: raise ValueError('Metadata response too large')
        rows = json.loads(data)
        if len(rows) != 1 or rows[0]['run_accession'] != accession: raise ValueError('Accession mismatch')
        row = rows[0]
        for key in ('fastq_ftp', 'fastq_md5', 'fastq_bytes'):
            if not row.get(key) or ';' in row[key]: raise ValueError('Expected one FASTQ per Yusa run')
        hostpath = row['fastq_ftp'].removeprefix('ftp://')
        if not hostpath.startswith('ftp.sra.ebi.ac.uk/'): raise ValueError('Unexpected archive host')
        row['url'] = 'https://' + hostpath
        row['metadata_url'] = url
        meta.append(row)
    save(ROOT / 'download-plan.json', meta)
    if sum(int(x['fastq_bytes']) for x in meta) > 2000000000: raise ValueError('2 GB download budget exceeded')
    attempts = []
    for url in LIB_URLS:
        try:
            with request(url) as f: data = f.read(16777217)
            if len(data) > 16777216: raise ValueError('Library too large')
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                entry = z.getinfo('yusa_library.csv')
                if entry.file_size > 16777216: raise ValueError('Library member too large')
                source = z.read(entry)
            break
        except Exception as e:
            attempts.append({'url': url, 'error': str(e)})
    else:
        save(ROOT / 'library-failures.json', attempts)
        raise RuntimeError('No real library obtained')
    (raw / 'library.source.csv').write_bytes(source)
    (raw / 'library.source.zip').write_bytes(data)
    rows = list(csv.DictReader(io.StringIO(source.decode('utf-8-sig'))))
    if len(rows) != 87437: raise ValueError('Wrong library size')
    ids = [r['id'] for r in rows]
    if len(set(ids)) != len(ids): raise ValueError('Duplicate target IDs')
    library = raw / 'library.tsv'
    with library.open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(['id', 'sequence', 'gene'])
        for r in rows:
            s = r['gRNA.sequence'].upper()
            if len(s) != 19 or not set(s) <= set('ACGT'): raise ValueError('Unexpected Yusa sequence')
            w.writerow([r['id'], s, r['Gene']])
    save(ROOT / 'library-provenance.json', {'source_url': url, 'failed_attempts': attempts,
      'csv_sha256': digest(raw / 'library.source.csv'), 'derived_sha256': digest(library), 'rows': len(rows)})
    from dotmatch.sensitivity import run_sensitivity
    import dotmatch
    samples = []
    for row in meta:
        accession = row['run_accession']; path = raw / (accession + '.fastq.gz')
        total = 0; md5 = hashlib.md5(usedforsecurity=False)
        with request(row['url']) as response, path.with_suffix('.pending').open('xb') as f:
            for b in iter(lambda: response.read(1048576), b''):
                total += len(b)
                if total > int(row['fastq_bytes']): raise ValueError('Archive size overrun')
                md5.update(b); f.write(b)
        if total != int(row['fastq_bytes']) or md5.hexdigest() != row['fastq_md5'].lower():
            raise ValueError('Archive size or MD5 mismatch')
        path.with_suffix('.pending').replace(path)
        prefix = raw / (accession + '.prefix100k.fastq.gz')
        n = 0
        with prefix.open('xb') as base, gzip.GzipFile(filename='', fileobj=base, mode='wb', mtime=0) as out:
            for n, lines in fastq(path):
                if n <= 100000: out.write(b''.join(lines))
        if n != int(row['read_count']): raise ValueError('Full record count differs from ENA')
        sample = dict(row, observed_records=n, sha256=digest(path), full_archive_md5_verified=True,
                      prefix_sha256=digest(prefix), prefix_records=min(100000,n))
        save(ROOT / (accession + '.provenance.json'), sample)
        samples.append(sample)
        for scope, reads in (('full', path), ('prefix100k', prefix)):
            output = ROOT / accession / scope
            output.parent.mkdir(exist_ok=True)
            result = run_sensitivity(targets=library, reads=reads, target_start=23, target_length=19,
               sample_label=accession, out_dir=output, write_read_changes=scope != 'full')
            if result['read_count'] != (n if scope == 'full' else min(n,100000)):
                raise AssertionError('Native record total differs')
            if result['inputs']['reads']['sha256'] != digest(reads): raise AssertionError('Native hash differs')
        print(accession, n, 'verified and counted', flush=True)
    save(ROOT / 'completion.json', {'completion': 'complete', 'samples': samples,
      'engine_version': dotmatch.__version__, 'python': sys.version, 'platform': platform.platform(),
      'scope': 'full_archive_and_separately_labelled_prefix', 'biological_inference': False,
      'files': {str(p.relative_to(ROOT)): digest(p) for p in sorted(ROOT.rglob('*')) if p.is_file()}})

if __name__ == '__main__':
    try: main()
    except Exception:
        save(ROOT / 'failure.json', {'completion': 'failed', 'traceback': traceback.format_exc(), 'synthetic_fallback': False})
        raise
