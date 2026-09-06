#!/usr/bin/env python3
"""Execute the new ESC2 archive and actual MAGeCK; no simulated biology."""
from __future__ import annotations
import argparse, csv, hashlib, json, os, platform, subprocess, sys, traceback, random
from collections import Counter
from functools import lru_cache
from pathlib import Path
import numpy as np

POLICIES=('exact','radius_k1','best_k1')

def sha(p,algorithm='sha256'):
    h=hashlib.new(algorithm)
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
def read(p):
    with Path(p).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,header,rows):
    with Path(p).open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(header);w.writerows(rows)
def run(cmd,log,env=None):
    with Path(log).open('w') as f:subprocess.run(cmd,check=True,stdout=f,stderr=subprocess.STDOUT,env=env)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);args=ap.parse_args();root=Path(args.root).resolve();work=root/'work';out=work/'results/replicated';out.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(root/'prior/source/research/assignment_robustness'));import replay
    from dotmatch.sensitivity import run_sensitivity
    resources=root/'resources';lib=root/'raw/library.tsv';reads=resources/'ERR377000.fastq.gz';binary=root/'upstream/guide-counter'
    p=json.loads((resources/'ERR377000.provenance.json').read_text())
    if sha(reads)!=p['sha256'] or sha(reads,'md5')!=p['fastq_md5'] or reads.stat().st_size!=int(p['fastq_bytes']):raise ValueError('Raw archive digest mismatch')
    if sha(lib)!='252e3b81b809c50f5cc347238a52926818027ad78a3ec98686e8012a8a46a896':raise ValueError('Reference digest mismatch')
    if sha(binary)!='96602fd0b9732204b530afb912ff679d48b0ba9e13d32c5eea67c10cbbdbf777':raise ValueError('Upstream binary mismatch')
    rows=replay.library(lib);n=len(rows);print('Verified complete original archive and reference',flush=True)
    summary=run_sensitivity(targets=lib,reads=reads,target_start=23,target_length=19,sample_label='ESC2',out_dir=out/'native')
    if summary['read_count']!=int(p['read_count']):raise AssertionError('Native full-read budget mismatch')
    cmd=[str(binary),'count','--input',str(reads),'--samples','ESC2','--library',str(lib),'--output',str(out/'upstream'),'--offset-sample-size','100000','--offset-min-fraction','0.0025']
    run(cmd,out/'upstream-execution.log');print('Native and actual upstream complete',flush=True)
    index=replay.Index([r['sequence'] for r in rows]);offsets=replay.discover(reads,index,23)['selected']
    calls=lru_cache(maxsize=262144)(index.calls);up=lru_cache(maxsize=262144)(index.upstream)
    counts=np.zeros((4,n),dtype=np.int64);states=[Counter() for _ in POLICIES];eventstats=Counter();classes=Counter()
    checks=set(random.Random(20260906).sample(range(1,int(p['read_count'])+1),200));checked=0
    for ordinal,rid,seq in replay.fastq(reads):
        window=seq[23:42] if len(seq)>=42 else None
        ca=calls(window)
        if ordinal in checks:
            if ca!=replay.oracle(index.sequences,window):raise AssertionError('All-target oracle discrepancy')
            checked+=1
        for j,c in enumerate(ca):
            if c>=0:counts[j,c]+=1;states[j]['unique']+=1
            else:states[j][replay.STATE[c]]+=1
        events=[]
        for off in offsets:
            if len(seq)>=off+19:
                target=up(seq[off:off+19])
                if target is not None:events.append((off,target));counts[3,target]+=1
        replay.add_events(eventstats,events,rows)
        distinct=tuple(sorted({i for _,i in events}))
        if len(distinct)>1:classes[distinct]+=1
        if ordinal%2000000==0:print('Independently reconciled',ordinal,flush=True)
    if ordinal!=int(p['read_count']) or checked!=200:raise AssertionError('Incomplete audit')
    for j,pol in enumerate(POLICIES):
        native=replay.read_counts(out/'native'/f'{pol}.counts.tsv',rows,'ESC2')
        if list(counts[j])!=native:raise AssertionError('Native complete-table discrepancy')
    if list(counts[3])!=replay.read_counts(out/'upstream.counts.txt',rows,'ESC2'):raise AssertionError('Actual upstream complete-table discrepancy')
    replay.check_conservation(eventstats)
    write(out/'all-guide-counts.tsv',['id','gene',*POLICIES,'upstream_count_events'],((r['id'],r['gene'],*(int(counts[j,i]) for j in range(4))) for i,r in enumerate(rows)))
    write(out/'multiple-guide-read-classes.tsv',['target_ids','genes','read_count'],((json.dumps([rows[i]['id'] for i in ids]),json.dumps(sorted({rows[i]['gene'] for i in ids})),num) for ids,num in sorted(classes.items())))
    dump(out/'raw-audit.json',dict(completion='complete',input=p,library_sha256=sha(lib),offsets=offsets,records=ordinal,full_guide_count_cells=4*n,all_target_checked_records=checked,disagreements=0,statistics=eventstats,source_sha256=sha(__file__),upstream_command=cmd))
    if sha(reads)!=p['sha256']:raise AssertionError('Input changed')
    # Assemble same original guide/annotation rows for all three samples.
    old=[]
    for accession in ('ERR376998','ERR376999'):
        folder=root/'prior/evidence/prior'/('assignment-replay-'+accession)/'audit-output'
        m=json.loads((folder/'completion.json').read_text())
        if sha(folder/'all-guide-counts.tsv')!=m['files']['all-guide-counts.tsv']:raise AssertionError('Old counts digest mismatch')
        rr=read(folder/'all-guide-counts.tsv')
        if [(x['id'],x['gene']) for x in rr]!=[(r['id'],r['gene']) for r in rows]:raise AssertionError('Guide identity mismatch')
        old.append(rr)
    mageck=resources/'mageck-source/bin/mageck';env=dict(os.environ);env['PYTHONPATH']=str(resources/'mageck-source');env['PATH']=str(resources/'mageck-source/bin')+os.pathsep+env.get('PATH','')
    commands=[]
    for j,pol in enumerate((*POLICIES,'upstream_count_events')):
        table=out/(pol+'.counts.tsv')
        write(table,['sgRNA','Gene','plasmid','ESC1','ESC2'],((r['id'],r['gene'],int(old[0][i][pol]),int(old[1][i][pol]),int(counts[j,i])) for i,r in enumerate(rows)))
        prefix=out/(pol+'.mageck')
        cmd=[sys.executable,str(mageck),'test','-k',str(table),'-c','plasmid','-t','ESC1,ESC2','--norm-method','median','--adjust-method','fdr','--gene-lfc-method','median','--keep-tmp','-n',str(prefix)]
        run(cmd,out/(pol+'.mageck-execution.log'),env);commands.append(cmd)
        if not Path(str(prefix)+'.gene_summary.txt').exists():raise AssertionError('MAGeCK result missing')
        print('Actual MAGeCK complete:',pol,flush=True)
    dump(out/'completion.json',dict(completion='complete',scope='actual_MAGeCK_0.5.9.5_original_ESC1_and_new_ESC2_vs_plasmid',protocol_commit='c41be48cbea1a330ed93c3dae79922f5824e5f67',python=sys.version,platform=platform.platform(),numpy=np.__version__,commands=commands,mageck_source_archive_sha256=sha(resources/'mageck-v0.5.9.5.tar.gz'),rra_sha256=sha(resources/'mageck-source/bin/RRA'),files={str(f.relative_to(out)):sha(f) for f in sorted(out.rglob('*')) if f.is_file()}))

if __name__=='__main__':main()
