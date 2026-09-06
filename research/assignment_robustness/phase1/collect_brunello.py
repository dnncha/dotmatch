#!/usr/bin/env python3
"""Complete, verified Brunello plasmid acquisition; no synthetic fallback."""
from __future__ import annotations
import csv, gzip, hashlib, io, json, platform, sys, traceback, urllib.parse, zipfile
from pathlib import Path
from collect import digest, save, request, fastq
ROOT=Path('brunello-output')
ACCESSION='SRR8297997'
LIB_MEMBER='broadgpp-brunello-library-corrected.txt'
LIB_URL='https://sourceforge.net/projects/mageck/files/libraries/'+LIB_MEMBER+'.zip/download'

def main():
    ROOT.mkdir(exist_ok=False);raw=ROOT/'raw';raw.mkdir()
    url='https://www.ebi.ac.uk/ena/portal/api/filereport?'+urllib.parse.urlencode({
      'accession':ACCESSION,'result':'read_run','format':'json',
      'fields':'run_accession,study_accession,sample_accession,sample_alias,experiment_title,library_layout,fastq_ftp,fastq_md5,fastq_bytes,read_count'})
    with request(url) as f: rows=json.loads(f.read(2000000))
    if len(rows)!=1 or rows[0]['run_accession']!=ACCESSION:raise ValueError('Wrong accession')
    meta=rows[0];meta['metadata_url']=url
    if meta['library_layout']!='SINGLE':raise ValueError('Expected single-end archive')
    for field in ('fastq_ftp','fastq_md5','fastq_bytes'):
        if not meta.get(field) or ';' in meta[field]:raise ValueError('Expected one FASTQ')
    if int(meta['fastq_bytes'])>500000000:raise ValueError('Declared 500 MB acquisition budget exceeded')
    host=meta['fastq_ftp'].removeprefix('ftp://')
    if not host.startswith('ftp.sra.ebi.ac.uk/'):raise ValueError('Unexpected archive host')
    meta['url']='https://'+host;save(ROOT/'download-plan.json',meta)
    with request(LIB_URL) as f: data=f.read(16777217)
    if len(data)>16777216:raise ValueError('Library archive too large')
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        wanted=[n for n in z.namelist() if Path(n).name==LIB_MEMBER]
        if len(wanted)!=1 or z.getinfo(wanted[0]).file_size>16777216:raise ValueError('Invalid library archive')
        source=z.read(wanted[0])
    (raw/'library.source.txt').write_bytes(source);(raw/'library.source.zip').write_bytes(data)
    text=source.decode('utf-8-sig');first=text.splitlines()[0]
    sep='\t' if '\t' in first else ','
    rows=list(csv.reader(io.StringIO(text),delimiter=sep));header=None
    if len(rows[0])<3:raise ValueError('Expected three library columns')
    if len(rows[0][1])!=20 or not set(rows[0][1].upper())<=set('ACGT'):
        header=rows.pop(0)
    if len(rows)!=77441:raise ValueError(f'Wrong reference size: {len(rows)}; header={header}')
    if len({r[0] for r in rows})!=len(rows):raise ValueError('Duplicate target IDs')
    library=raw/'library.tsv'
    with library.open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['id','sequence','gene'])
        for r in rows:
            if len(r)<3 or not r[0] or not r[2] or len(r[1])!=20 or not set(r[1].upper())<=set('ACGT'):
                raise ValueError('Invalid library row')
            w.writerow([r[0],r[1].upper(),r[2]])
    save(ROOT/'library-provenance.json',{'source_url':LIB_URL,'source_header':header,
       'source_sha256':digest(raw/'library.source.txt'),'derived_sha256':digest(library),'rows':len(rows)})
    path=raw/(ACCESSION+'.fastq.gz');h=hashlib.md5(usedforsecurity=False);total=0
    with request(meta['url']) as response,path.with_suffix('.pending').open('xb') as f:
        for b in iter(lambda:response.read(1048576),b''):
            total+=len(b)
            if total>int(meta['fastq_bytes']):raise ValueError('Size overrun')
            h.update(b);f.write(b)
    if total!=int(meta['fastq_bytes']) or h.hexdigest()!=meta['fastq_md5'].lower():raise ValueError('Size/MD5 mismatch')
    path.with_suffix('.pending').replace(path)
    prefix=raw/(ACCESSION+'.prefix100k.fastq.gz');n=0
    with prefix.open('xb') as base,gzip.GzipFile(filename='',fileobj=base,mode='wb',mtime=0) as out:
        for n,lines in fastq(path):
            if n<=100000:out.write(b''.join(lines))
    if n!=int(meta['read_count']):raise ValueError('Record count differs from ENA')
    meta.update(observed_records=n,sha256=digest(path),full_archive_md5_verified=True,
       prefix_sha256=digest(prefix),prefix_records=min(n,100000))
    save(ROOT/(ACCESSION+'.provenance.json'),meta)
    from dotmatch.sensitivity import run_sensitivity
    import dotmatch
    for scope,reads in (('full',path),('prefix100k',prefix)):
        output=ROOT/ACCESSION/scope;output.parent.mkdir(exist_ok=True)
        result=run_sensitivity(targets=library,reads=reads,target_start=21,target_length=20,
          sample_label=ACCESSION,out_dir=output,write_read_changes=scope!='full')
        if result['read_count']!=(n if scope=='full' else min(n,100000)):raise AssertionError('Native read count mismatch')
        if result['inputs']['reads']['sha256']!=digest(reads):raise AssertionError('Native input hash differs')
    save(ROOT/'completion.json',{'completion':'complete','sample':meta,'python':sys.version,
       'platform':platform.platform(),'engine_version':dotmatch.__version__,'biological_inference':False,
       'files':{str(p.relative_to(ROOT)):digest(p) for p in sorted(ROOT.rglob('*')) if p.is_file()}})
    print(ACCESSION,n,'verified and counted',flush=True)

if __name__=='__main__':
    try:main()
    except Exception:
        save(ROOT/'failure.json',{'completion':'failed','traceback':traceback.format_exc(),'synthetic_fallback':False})
        raise
