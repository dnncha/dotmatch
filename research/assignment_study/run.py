#!/usr/bin/env python3
"""Complete public-data measurement study with independent accounting and gates."""
from __future__ import annotations
import argparse, concurrent.futures as cf, csv, hashlib, io, json, pathlib
import platform, shutil, subprocess, sys, time, urllib.request, zipfile
BASE='5c3934d648ff55928db5a9219d7703845941d28e'
UPSTREAM='f24de175282064773328e35baf03373e7e3b895d'
LIB_SHA='d41a4122a46fdedb47deca61081dde9037d461a1fce1e6f9744522e21005ceba'
RUNS=[
 ('ERR376998','plasmid','A','SAMEA2188756',10093905,282339185,'7b552f8efd6fe90fc6bc34879c617040'),
 ('ERR376999','ESC1','A','SAMEA2188757',10300758,290645388,'be96811adce2f5d69efd18a02f8e2416'),
 ('ERR377000','ESC2','A','SAMEA2188758',10820594,305692488,'998646750e84033902d2239a329b6745'),
 ('ERR377001','plasmid','B','SAMEA2188756',10246121,291385291,'3abe7d2d1d38cbda76dc7c27f5dc7ea5'),
 ('ERR377002','ESC1','B','SAMEA2188757',10393396,297070739,'456d1f4b02333a990513a1cb2b90f5c0'),
 ('ERR377003','ESC2','B','SAMEA2188758',10976838,314600730,'02b021412e06a57886551fc66b4fd2e6')]
def digest(p,algorithm='sha256'):
 h=hashlib.new(algorithm)
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def retrieve(url,p,expected_md5=None,expected_size=None):
 p=pathlib.Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 if p.exists() and expected_md5 and digest(p,'md5')==expected_md5:return
 for attempt in range(3):
  tmp=p.with_name(p.name+'.partial')
  try:
   req=urllib.request.Request(url,headers={'User-Agent':'DotMatch-research-reproducibility/1.0'})
   with urllib.request.urlopen(req,timeout=120) as r,open(tmp,'wb') as f:shutil.copyfileobj(r,f,1024*1024)
   if expected_size is not None and tmp.stat().st_size!=expected_size:raise ValueError('download length mismatch')
   if expected_md5 and digest(tmp,'md5')!=expected_md5:raise ValueError('download MD5 mismatch')
   tmp.replace(p);return
  except Exception:
   tmp.unlink(missing_ok=True)
   if attempt==2:raise
   time.sleep(3)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=pathlib.Path,required=True);a=ap.parse_args()
 out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
 if (out/'COMPLETE.json').exists():raise RuntimeError('Refusing to overwrite completed study')
 src=pathlib.Path(__file__).resolve().parent;root=src.parents[1]
 work=root/'research-work';raw=work/'raw';raw.mkdir(parents=True,exist_ok=True)
 for sub in ['logs','oracle','native','upstream','metadata','source']:(out/sub).mkdir(exist_ok=True)
 for f in src.glob('*'):
  if f.is_file() and f.suffix in ['.py','.cpp','.md','.json']:shutil.copy2(f,out/'source'/f.name)
 commands=[]
 def run(label,cmd,cwd=None,timeout=1500):
  cmd=[str(x) for x in cmd];t=time.monotonic();print('START',label,flush=True);log=out/'logs'/(label+'.txt')
  with open(log,'w') as f:
   f.write(json.dumps(cmd)+'\n');f.flush();p=subprocess.run(cmd,cwd=cwd or root,stdout=f,stderr=subprocess.STDOUT,timeout=timeout)
  record={'label':label,'argv':cmd,'cwd':str(cwd or root),'returncode':p.returncode,'wall_seconds':time.monotonic()-t}
  commands.append(record);print('END',label,p.returncode,round(record['wall_seconds'],2),flush=True)
  if p.returncode:print(log.read_text()[-12000:],flush=True);raise RuntimeError('Failed: '+label)
 manifest=[]
 for acc,sample,lane,bio,n,size,md5 in RUNS:
  url=f'https://ftp.sra.ebi.ac.uk/vol1/fastq/{acc[:6]}/{acc}/{acc}.fastq.gz'
  manifest.append(dict(accession=acc,sample=sample,technical_set=lane,biosample=bio,expected_reads=n,expected_bytes=size,expected_md5=md5,url=url))
 (out/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 for bio in sorted({x[3] for x in RUNS}):retrieve('https://www.ebi.ac.uk/ena/browser/api/xml/'+bio,out/'metadata'/(bio+'.xml'))
 retrieve('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB4038&result=read_run&fields=run_accession,sample_accession,experiment_accession,read_count,fastq_ftp,fastq_bytes,fastq_md5&format=json',out/'metadata'/'study_runs.json')
 libzip=raw/'yusa_library.csv.zip';retrieve('https://sourceforge.net/projects/mageck/files/libraries/yusa_library.csv.zip/download',libzip)
 with zipfile.ZipFile(libzip) as z:libbytes=z.read('yusa_library.csv')
 if hashlib.sha256(libbytes).hexdigest()!=LIB_SHA:raise RuntimeError('Original library changed')
 (out/'metadata'/'yusa_library.original.csv').write_bytes(libbytes)
 rows=list(csv.reader(io.StringIO(libbytes.decode())))[1:]
 assert len(rows)==87437 and len({r[0] for r in rows})==len(rows) and len({r[1] for r in rows})==len(rows)
 assert all(len(r[1])==19 and set(r[1])<=set('ACGT') and len(r)>=3 for r in rows)
 lib=out/'metadata'/'library.tsv'
 with open(lib,'w') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['id','sequence','gene']);w.writerows([r[:3] for r in rows])
 def download(m):
  p=raw/(m['accession']+'.fastq.gz');retrieve(m['url'],p,m['expected_md5'],m['expected_bytes'])
  m['observed_bytes']=p.stat().st_size;m['sha256']=digest(p);m['md5_verified']=True
  print('VERIFIED INPUT',m['accession'],m['observed_bytes'],flush=True)
 def build():
  run('build-dotmatch',['make','dotmatch'])
  run('build-oracle',['g++','-std=c++17','-O3','-Wall','-Wextra',src/'oracle.cpp','-lz','-o',src/'oracle'])
  run('test-oracle',[sys.executable,src/'test_oracle.py'])
  up=work/'guide-counter'
  if not up.exists():run('clone-guide-counter',['git','clone','https://github.com/fulcrumgenomics/guide-counter.git',up])
  run('pin-guide-counter',['git','checkout','--detach',UPSTREAM],cwd=up)
  run('build-guide-counter',['cargo','build','--release','--locked'],cwd=up)
 try:
  with cf.ThreadPoolExecutor(max_workers=4) as ex:
   tasks=[ex.submit(download,m) for m in manifest]+[ex.submit(build)]
   for f in cf.as_completed(tasks):f.result()
  (out/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
  binaries={'dotmatch':root/'dotmatch','oracle':src/'oracle','guide_counter':work/'guide-counter'/'target'/'release'/'guide-counter'}
  versions={'source_base':BASE,'source_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),'guide_counter_commit':UPSTREAM,'python':sys.version,'platform':platform.platform(),'binary_sha256':{k:digest(p) for k,p in binaries.items()},'library_sha256':LIB_SHA}
  (out/'versions.json').write_text(json.dumps(versions,indent=2)+'\n')
  jobs=[]
  for m in manifest:
   acc=m['accession'];jobs.append(('oracle-'+acc,[binaries['oracle'],lib,raw/(acc+'.fastq.gz'),out/'oracle'/acc,'100000']))
  for label,k,policy in [('exact',0,'radius'),('radius1',1,'radius'),('best1',1,'best')]:
   cmd=[binaries['dotmatch'],'count','--targets',lib]
   for m in manifest:cmd+=['--reads',raw/(m['accession']+'.fastq.gz')]
   cmd+=['--sample-label',','.join(m['accession'] for m in manifest),'--target-start','23','--target-length','19','--metric','hamming','--k',str(k),'--ambiguity-policy',policy,'--ambiguous','discard','--backend','cpu','--format','mageck','--out',out/'native'/(label+'.counts.tsv'),'--summary',out/'native'/(label+'.summary.json')]
   jobs.append(('native-'+label,cmd))
  for label in ['exact','best1']:
   cmd=[binaries['guide_counter'],'count','--input']+[raw/(m['accession']+'.fastq.gz') for m in manifest]+['--samples']+[m['accession'] for m in manifest]+['--library',lib,'--output',out/'upstream'/label,'--offset-sample-size','100000','--offset-min-fraction','0.0025']
   if label=='exact':cmd+=['--exact-match']
   jobs.append(('upstream-'+label,cmd))
  with cf.ThreadPoolExecutor(max_workers=2) as ex:
   fut=[ex.submit(run,label,cmd) for label,cmd in jobs]
   for f in cf.as_completed(fut):f.result()
  def readtable(p):
   with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
  native={p:readtable(out/'native'/(p+'.counts.tsv')) for p in ['exact','radius1','best1']}
  upstream={p:readtable(out/'upstream'/(p+'.counts.txt')) for p in ['exact','best1']}
  validation=[]
  for m in manifest:
   acc=m['accession'];oracle=readtable(out/'oracle'/(acc+'.counts.tsv'));byid={r['guide']:r for r in oracle}
   qc=readtable(out/'oracle'/(acc+'.qc.tsv'));n=next(int(r['value']) for r in qc if r['mode']=='all' and r['metric']=='reads')
   if n!=m['expected_reads']:raise RuntimeError('FASTQ record count mismatch '+acc)
   for kind,tables in [('dotmatch',native),('guide_counter',upstream)]:
    for policy,table in tables.items():
     column=('fixed_' if kind=='dotmatch' else 'multi_')+policy;diffs=[];seen=set()
     for r in table:
      gid=r.get('sgRNA',r.get('guide'));seen.add(gid)
      if gid not in byid or int(r[acc])!=int(byid[gid][column]):diffs.append({'guide':gid,'observed':r[acc],'oracle':byid.get(gid,{}).get(column)})
     if seen!=set(byid):diffs.append({'guide_set_mismatch':True})
     validation.append({'run':acc,'tool':kind,'policy':policy,'guides_checked':len(byid),'differences':len(diffs),'examples':diffs[:20]})
  (out/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
  if any(x['differences'] for x in validation):raise RuntimeError('Independent parity gate failed; biological interpretation blocked')
  total=sum(m['expected_reads'] for m in manifest)
  print('RESULT: COMPLETE RAW DATA ACCOUNTING',total,'reads',len(rows),'guides; all 30 full-matrix parity checks passed',flush=True)
  (out/'COMPLETE.json').write_text(json.dumps({'status':'complete','total_reads':total,'library_guides':len(rows),'comparisons':len(validation),'discrepancies':0,'scope':'six complete sequencing runs; three sample identities; not six biological replicates'},indent=2)+'\n')
 finally:
  (out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
  with open(out/'SHA256SUMS','w') as f:
   for p in sorted(out.rglob('*')):
    if p.is_file() and p.name!='SHA256SUMS':f.write(digest(p)+'  '+str(p.relative_to(out))+'\n')
if __name__=='__main__':main()
