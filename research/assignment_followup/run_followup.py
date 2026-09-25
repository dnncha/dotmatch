#!/usr/bin/env python3
"""Exploratory sequence/quality follow-up with original-count parity."""
import argparse,collections,concurrent.futures as cf,csv,hashlib,json,math,pathlib,random,shutil,subprocess,sys,time
HERE=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'assignment_study'))
import run as core

def table(p):
 with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def main():
 p=argparse.ArgumentParser();p.add_argument('--baseline',type=pathlib.Path,required=True);p.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args();base=a.baseline.resolve();out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
 assert json.loads((base/'COMPLETE.json').read_text())['discrepancies']==0
 for line in (base/'SHA256SUMS').read_text().splitlines():
  h,path=line.split('  ',1)
  if core.digest(base/path)!=h:raise RuntimeError('Baseline artifact checksum mismatch')
 (out/'source').mkdir(exist_ok=True)
 for f in HERE.glob('*'):
  if f.is_file() and f.suffix in ['.py','.cpp','.md']:shutil.copy2(f,out/'source'/f.name)
 manifest=json.loads((base/'input_manifest.json').read_text());lib=base/'metadata'/'library.tsv';guides=table(lib);byid={r['id']:r for r in guides}
 original={m['accession']:{r['guide']:r for r in table(base/'oracle'/(m['accession']+'.counts.tsv'))} for m in manifest}
 plasmid=[m['accession'] for m in manifest if m['sample']=='plasmid'];counts={g:sum(int(original[acc][g]['fixed_exact']) for acc in plasmid) for g in byid}
 gc=lambda seq:sum(b in 'GC' for b in seq)/len(seq)
 target_genes={'RPL4','EMG1'};targets=sorted(g for g,r in byid.items() if r['gene'] in target_genes);assert len(targets)==10
 selected=[dict(id=g,gene=byid[g]['gene'],role='selected_gene',matched_to=g,pooled_plasmid_exact=counts[g],gc=gc(byid[g]['sequence'])) for g in targets]
 rng=random.Random(20260925);used=set(targets);used_genes=set(target_genes)
 for g in targets:
  choices=sorted(h for h,r in byid.items() if h not in used and r['gene'] not in used_genes and 0.5*counts[g]<=counts[h]<=2*counts[g] and abs(gc(r['sequence'])-gc(byid[g]['sequence']))<=2/19+1e-12)
  rng.shuffle(choices);picked=0
  for h in choices:
   if byid[h]['gene'] in used_genes:continue
   used.add(h);used_genes.add(byid[h]['gene']);selected.append(dict(id=h,gene=byid[h]['gene'],role='matched_control',matched_to=g,pooled_plasmid_exact=counts[h],gc=gc(byid[h]['sequence'])));picked+=1
   if picked==3:break
  if picked!=3:raise RuntimeError('Insufficient matched controls')
 assert len(selected)==40
 selection=out/'selection.tsv'
 with selection.open('w') as f:w=csv.DictWriter(f,fieldnames=list(selected[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(selected)
 logs=out/'logs';logs.mkdir(exist_ok=True);records=[]
 def run(label,cmd):
  t=time.monotonic();print('START',label,flush=True)
  with (logs/(label+'.txt')).open('w') as f:
   f.write(json.dumps(list(map(str,cmd)))+'\n');f.flush();r=subprocess.run(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT)
  records.append(dict(label=label,command=list(map(str,cmd)),returncode=r.returncode,wall_seconds=time.monotonic()-t));print('END',label,r.returncode,flush=True)
  if r.returncode:raise RuntimeError('Failed '+label)
 raw=HERE.parents[1]/'research-work'/'raw';raw.mkdir(parents=True,exist_ok=True)
 try:
  run('build',['g++','-O3','-std=c++17',HERE/'variants.cpp','-lz','-o',HERE/'variants'])
  def process(m):
   acc=m['accession'];f=raw/(acc+'.fastq.gz');core.retrieve(m['url'],f,m['expected_md5'],m['expected_bytes']);assert core.digest(f)==m['sha256']
   run(acc,[HERE/'variants',lib,selection,f,out/(acc+'.variants.tsv'),m['expected_reads']])
  with cf.ThreadPoolExecutor(max_workers=3) as ex:
   futures=[ex.submit(process,m) for m in manifest]
   for f in cf.as_completed(futures):f.result()
  checks=[]
  for m in manifest:
   acc=m['accession'];observed=collections.Counter()
   for r in table(out/(acc+'.variants.tsv')):observed[r['guide'],int(r['distance'])]+=int(r['count'])
   for s in selected:
    g=s['id'];old=original[acc][g]
    for distance,expected in [(0,int(old['fixed_exact'])),(1,int(old['fixed_best1'])-int(old['fixed_exact']))]:
     got=observed[g,distance];checks.append(dict(run=acc,guide=g,distance=distance,expected=expected,observed=got,pass_check=got==expected))
  (out/'validation.json').write_text(json.dumps(checks,indent=2)+'\n')
  if not all(x['pass_check'] for x in checks):raise RuntimeError('Original selected-guide count parity failed')
  (out/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
  (out/'COMPLETE.json').write_text(json.dumps(dict(status='complete',checks=len(checks),guides=len(selected),case_guides=10,control_guides=30,total_reads=sum(m['expected_reads'] for m in manifest),scope='post-result exploratory sequence follow-up; same six source files, not new biological replication'),indent=2)+'\n')
  print('COMPLETE',len(checks),'original-count checks,',len(selected),'guides',flush=True)
 finally:
  (out/'commands.json').write_text(json.dumps(records,indent=2)+'\n')
  with (out/'SHA256SUMS').open('w') as f:
   for p in sorted(out.rglob('*')):
    if p.is_file() and p.name!='SHA256SUMS':f.write(core.digest(p)+'  '+str(p.relative_to(out))+'\n')
if __name__=='__main__':main()
