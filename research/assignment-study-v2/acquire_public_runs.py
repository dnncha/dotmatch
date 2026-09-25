#!/usr/bin/env python3
"""Full public-run acquisition. No biological truth or accuracy claims.
Keeps registry metadata, complete archive digests, all count matrices, and a
hash-selected audit sample. The original raw archive is NOT an output artifact.
"""
from __future__ import annotations
import argparse, collections, csv, datetime, gzip, hashlib, json, os
from pathlib import Path
import shutil, subprocess, sys, time, urllib.parse, urllib.request

LIB_SHA='d41a4122a46fdedb47deca61081dde9037d461a1fce1e6f9744522e21005ceba'
LIB_URL='https://raw.githubusercontent.com/dnncha/dotmatch/0decef38c4380edde5436f47fe84271f5f6d4a46/benchmarks/real/data/yusa_library.csv'
WHEEL='dotmatch-0.6.0-py3-none-manylinux2014_x86_64.manylinux_2_17_x86_64.manylinux_2_28_x86_64.whl'
WHEEL_SHA='2f9fbaffe9965bce646432956920cf9aef4c47ffdc2a6b73d43b2471faf74f53'

def digest(p, kind='sha256'):
    h=hashlib.new(kind)
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()

def write(p,obj):
    Path(p).write_text(json.dumps(obj,indent=2)+'\n')

def download(url,path):
    request=urllib.request.Request(url,headers={'User-Agent':'DotMatch-public-methods-study/2'})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request,timeout=120) as r, Path(path).open('wb') as f:
                shutil.copyfileobj(r,f,1048576)
            return
        except Exception:
            if attempt==2: raise
            time.sleep(3*(attempt+1))

def run(accession,out):
    if accession not in {'ERR376998','ERR376999','ERR377000'}: raise ValueError('Unapproved accession')
    out.mkdir(parents=True,exist_ok=False)
    scratch=out.parent/(out.name+'-scratch'); scratch.mkdir(exist_ok=False)
    lib=out/'yusa_library.csv'; download(LIB_URL,lib)
    if digest(lib)!=LIB_SHA: raise ValueError('Library identity mismatch')
    wheel=scratch/WHEEL
    download('https://github.com/dnncha/dotmatch/releases/download/v0.6.0/'+WHEEL,wheel)
    if digest(wheel)!=WHEEL_SHA: raise ValueError('Wheel identity mismatch')
    subprocess.run([sys.executable,'-m','pip','install','--no-deps',str(wheel)],check=True)
    version=subprocess.check_output(['dotmatch','--version'],text=True).strip()
    if '0.6.0' not in version: raise ValueError('Wrong executable version')
    fields='run_accession,study_accession,sample_accession,experiment_accession,scientific_name,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes'
    url='https://www.ebi.ac.uk/ena/portal/api/filereport?'+urllib.parse.urlencode({'accession':accession,'result':'read_run','fields':fields,'format':'json'})
    download(url,out/'ena_registry.json')
    registry=json.loads((out/'ena_registry.json').read_text())
    if len(registry)!=1 or registry[0]['run_accession']!=accession: raise ValueError('Registry accession mismatch')
    row=registry[0]
    if ';' in row['fastq_ftp']: raise ValueError('Expected one single-end FASTQ archive')
    archive=scratch/(accession+'.fastq.gz')
    fqurl='https://'+row['fastq_ftp'].removeprefix('ftp://').removeprefix('https://')
    download(fqurl,archive)
    md5=digest(archive,'md5'); nbytes=archive.stat().st_size
    if md5!=row['fastq_md5'] or nbytes!=int(row['fastq_bytes']): raise ValueError('Complete archive identity mismatch')
    targets={x['gRNA.sequence'] for x in csv.DictReader(lib.open())}
    n=0; bases=0; sample_n=0; lengths=collections.Counter(); contexts=collections.Counter(); offsets=collections.Counter()
    sample=out/'audit_sample.fastq.gz'
    with gzip.open(archive,'rt') as f, gzip.open(sample,'wt') as audit:
        while True:
            name=f.readline()
            if not name: break
            seq=f.readline().rstrip('\r\n'); plus=f.readline(); quality=f.readline().rstrip('\r\n')
            if not name.startswith('@') or not plus.startswith('+') or len(seq)!=len(quality): raise ValueError('Malformed FASTQ')
            if n<50000:
                for start in range(max(0,len(seq)-18)):
                    if seq[start:start+19] in targets: offsets[start]+=1
            token=hashlib.blake2b(f'20260926:{accession}:{n}'.encode(),digest_size=8).digest()
            if int.from_bytes(token,'big') % 1000 == 0:
                audit.write(name+seq+'\n'+plus+quality+'\n'); sample_n+=1
            lengths[len(seq)]+=1; bases+=len(seq); n+=1
            if seq[23:42] in targets: contexts[(seq[:23],seq[42:])]+=1
    if n!=int(row['read_count']) or bases!=int(row['base_count']): raise ValueError('Registry read/base accounting mismatch')
    write(out/'input_audit.json',{'accession':accession,'records':n,'bases':bases,'archive_bytes':nbytes,'archive_md5':md5,'archive_sha256':digest(archive),'archive_md5_verified':True,'archive_gzip_eof_verified':True,'library_sha256':LIB_SHA,'wheel_sha256':WHEEL_SHA,'version':version,'registry_url':url,'fastq_url':fqurl,'length_counts':dict(lengths),'first_50000_exact_offsets':dict(offsets),'calibrated_offsets_at_1percent':[k for k,v in sorted(offsets.items()) if v/min(n,50000)>=.01],'top_exact23_contexts':[{'prefix':k[0],'suffix':k[1],'records':v} for k,v in contexts.most_common(20)],'audit_sample_records':sample_n,'audit_sample_sha256':digest(sample),'audit_sample_design':'BLAKE2b 64-bit of 20260926:accession:zero-based ordinal, integer modulo 1000 equals zero; no prefix restriction','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
    records=[]
    # Adjacent starts are the predeclared v1 experimental candidate set, not a claim about another tool's defaults.
    for policy,k,ambiguity,starts in [('exact',0,'radius',[22,23,24]),('best1',1,'best',[22,23,24]),('radius1',1,'radius',[23])]:
        for start in starts:
            stem=f'{policy}_offset{start}'; counts=out/(stem+'.counts.tsv'); summary=out/(stem+'.summary.json')
            cmd=['dotmatch','count','--targets',str(lib),'--reads',str(archive),'--target-start',str(start),'--target-length','19','--k',str(k),'--metric','hamming','--ambiguity-policy',ambiguity,'--out',str(counts),'--summary',str(summary),'--sample-label',accession,'--no-progress']
            begin=time.perf_counter(); p=subprocess.run(cmd,text=True,capture_output=True)
            (out/(stem+'.stderr.txt')).write_text(p.stderr)
            if p.returncode: raise RuntimeError('Native command failed: '+p.stderr)
            plain_sha=digest(counts)
            with counts.open('rb') as src,gzip.open(str(counts)+'.gz','wb') as dst: shutil.copyfileobj(src,dst)
            counts.unlink()
            records.append({'policy':policy,'start':start,'command':cmd,'returncode':p.returncode,'wall_seconds':time.perf_counter()-begin,'plain_count_sha256':plain_sha,'summary_sha256':digest(summary)})
            print(accession,stem,'complete',flush=True)
    write(out/'native_runs.json',records)
    write(out/'complete.json',{'complete':True,'accession':accession,'native_runs':len(records),'github_sha':os.getenv('GITHUB_SHA'),'github_run_id':os.getenv('GITHUB_RUN_ID'),'scope':'Actual full-accession count comparison; origin and biological accuracy unknown. No gene significance tests.'})
    shutil.rmtree(scratch)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('accession');p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    run(a.accession,a.out)
