#!/usr/bin/env python3
"""AR004 complete public-run execution. Never substitutes prefixes or synthetic data."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,os,platform,random,shutil,subprocess,sys,time,traceback,urllib.request
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent))
import replay
RUNS={'pDNA':['SRR8297997'],'A':['SRR8297836','SRR8297837'],'B':['SRR8297838','SRR8297839'],'C':['SRR8297840','SRR8297841']}
ALIASES={'pDNA':'Brunello_mod_tracr_pDNA',**{k:'Brunello_mod_tracr_Rep'+k+'_Dropout_A375' for k in 'ABC'}}
POLICIES=('event_exact','event_best','joint_exact','joint_best','joint_radius')
OFFSETS=tuple(range(21,31))
LIB_SHA='0d2906187829ea9f736de94a47369bd94d42cde5f348fea9d12a385625cc2ca1'
def sha(path,algorithm='sha256'):
 h=hashlib.new(algorithm)
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def dump(path,obj):
 Path(path).parent.mkdir(parents=True,exist_ok=True)
 Path(path).write_text(json.dumps(obj,sort_keys=True,indent=2,allow_nan=False)+'\n')
def table(path,key):
 with Path(path).open(newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
 result={r[key]:r for r in rows}
 if not rows or len(result)!=len(rows):raise ValueError('Missing/duplicate table identities')
 return result
def metadata():
 path=ROOT/'curation/PRJNA508200.json';rows=json.loads(path.read_text());lookup={r['run_accession']:r for r in rows};selected=[]
 for sample,accessions in RUNS.items():
  available={r['run_accession'] for r in rows if r['sample_alias']==ALIASES[sample]}
  if available!=set(accessions):raise ValueError('Not all technical runs included: '+sample)
  if len({lookup[a]['sample_accession'] for a in accessions})!=1:raise ValueError('Technical runs have different sample identities')
  for accession in accessions:
   r=dict(lookup[accession],biological_sample=sample)
   if r['library_layout']!='SINGLE' or any(';' in r[k] for k in ('fastq_ftp','fastq_md5','fastq_bytes')):raise ValueError('Unexpected paired or multipart archive')
   selected.append(r)
 if sum(int(r['fastq_bytes']) for r in selected)>8000000000:raise ValueError('8 GB transport budget exceeded')
 return selected

def download(row,path):
 source=row['fastq_ftp'].removeprefix('ftp://')
 if not source.startswith('ftp.sra.ebi.ac.uk/'):raise ValueError('Unexpected archive host')
 url='https://'+source
 if path.exists():
  if path.stat().st_size!=int(row['fastq_bytes']) or sha(path,'md5')!=row['fastq_md5']:raise ValueError('Existing raw input differs')
  return url
 for attempt in range(3):
  pending=path.with_suffix('.pending')
  try:
   total=0;h=hashlib.md5()
   with urllib.request.urlopen(url,timeout=120) as response,pending.open('wb') as out:
    for b in iter(lambda:response.read(1048576),b''):
     total+=len(b)
     if total>int(row['fastq_bytes']):raise ValueError('Archive size overrun')
     out.write(b);h.update(b)
   if total!=int(row['fastq_bytes']) or h.hexdigest()!=row['fastq_md5']:raise ValueError('Archive integrity mismatch')
   pending.rename(path);return url
  except Exception:
   pending.unlink(missing_ok=True)
   if attempt==2:raise
   time.sleep(2**attempt)

def validate_public(binary,lib,reads,N):
 import numpy as np
 import dotmatch
 rows=list(table(lib,'id').values());sequences=[r['sequence'] for r in rows]
 matrix=np.array([list(s.encode('ascii')) for s in sequences],dtype=np.uint8)
 ordinals=set(random.Random(20260906).sample(range(1,N+1),min(200,N)))
 selected=[]
 for ordinal,_,seq in replay.fastq(reads):
  if ordinal in ordinals:selected.append(seq)
  if ordinal==max(ordinals):break
 if len(selected)!=len(ordinals):raise AssertionError('Missing selected validation records')
 result=subprocess.run([str(binary),str(lib),','.join(map(str,OFFSETS)),'--probe'],input='\n'.join(selected)+'\n',text=True,capture_output=True,check=True)
 lines=result.stdout.splitlines()
 if len(lines)!=len(selected):raise AssertionError('Native probe count differs')
 windows=[];expected_calls=[];checks=0
 for seq,line in zip(selected,lines):
  fields=line.split('\t');got=tuple(tuple(map(int,h.split(','))) for h in fields[1].split(';') if h)
  valid=set(seq)<=set('ACGTRYSWKMBDHVN') and len(seq)>=50
  if int(fields[0])!=valid:raise AssertionError('Probe validity mismatch')
  expected=[]
  if valid:
   for offset in OFFSETS:
    w=seq[offset:offset+20];dist=np.count_nonzero(matrix!=np.frombuffer(w.encode('ascii'),dtype=np.uint8),axis=1)
    ids=np.flatnonzero(dist<=1);expected.extend((offset,int(i),int(dist[i])) for i in ids)
    best=int(min(dist)) if len(dist) else 2;nearest=[int(i) for i in ids if dist[i]==best]
    windows.append(w);expected_calls.append((len(ids),best,nearest));checks+=1
  if got!=tuple(expected):raise AssertionError('Exhaustive public-record candidate mismatch')
 with dotmatch.Matcher(sequences) as matcher: native=matcher.assign_hamming(windows,k=1,policy='best')
 if len(native)!=len(windows):raise AssertionError('Pinned DotMatch output length mismatch')
 for result,(n,best,nearest) in zip(native,expected_calls):
  state='none' if not n else 'unique' if len(nearest)==1 else 'ambiguous'
  if result.match_count!=n or dotmatch.status_name(result.status)!=state:raise AssertionError('Pinned DotMatch candidate/status mismatch')
  if n and result.best_distance!=best:raise AssertionError('Pinned distance mismatch')
  if len(nearest)==1 and result.target_index!=nearest[0]:raise AssertionError('Pinned unique identity mismatch')
 return {'selected_records':len(selected),'full_library_position_checks':checks,'pinned_DotMatch_checks':len(windows),'disagreements':0,'seed':20260906,'dotmatch_version':dotmatch.__version__}

def verify_tables(folder,lib):
 reference=table(lib,'id');counts=table(folder/'guide-counts.tsv','id');genes=table(folder/'gene-counts.tsv','gene')
 if set(counts)!=set(reference) or any(counts[i]['gene']!=reference[i]['gene'] for i in reference):raise AssertionError('Count/reference identity mismatch')
 metrics=json.loads((folder/'metrics.json').read_text());N=metrics['records']
 with (folder/'qc.tsv').open() as f:qc=list(csv.DictReader(f,delimiter='\t'))
 if len(qc)!=8:raise AssertionError('Missing read-state categories')
 for row in qc:
  q=[int(row[k]) for k in ('unique','ambiguous','none','invalid')]
  if any(x<0 for x in q) or sum(q)!=N:raise AssertionError('Read-state budget failed')
 for policy in POLICIES:
  values=[int(r[policy]) for r in counts.values()]
  thinned=[int(r['thin_'+policy]) for r in counts.values()]
  if any(a<0 or b<0 or b>a for a,b in zip(values,thinned)):raise AssertionError('Count/thinning budget invalid')
  if policy.startswith('event'):
   e=metrics['policies'][policy]
   if sum(values)!=e['count_events'] or e['count_events']!=e['matched_records']+e['extra_events']:raise AssertionError('Event conservation failed')
  else:
   q=next(r for r in qc if r['policy']==policy and r['resolution']=='guide')
   if sum(values)!=int(q['unique']):raise AssertionError('Joint guide budget failed')
   q=next(r for r in qc if r['policy']==policy and r['resolution']=='gene')
   if sum(int(r[policy+'_lower']) for r in genes.values())!=int(q['unique']):raise AssertionError('Joint gene budget failed')
 return {'guide_rows':len(counts),'gene_annotations':len(genes),'all_count_budgets':True}

def baseline_check(folder,previous):
 old=next(Path(previous).rglob('guide-counts.tsv'));manifest=json.loads((old.parent/'completion.json').read_text())
 for name,h in manifest['files'].items():
  if sha(old.parent/name)!=h:raise AssertionError('AR003 input evidence checksum failed')
 a=table(folder/'guide-counts.tsv','id');b=table(old,'id')
 if set(a)!=set(b):raise AssertionError('Baseline reference differs')
 mapping={'joint_exact':'exact','joint_best':'best_k1','joint_radius':'radius_k1'}
 for p,q in mapping.items():
  if any(a[i][p]!=b[i][q] for i in a):raise AssertionError('Full AR003 guide recount differs: '+p)
 ga=table(folder/'gene-counts.tsv','gene');gb=table(old.parent/'gene-counts-and-bounds.tsv','gene')
 if set(ga)!=set(gb):raise AssertionError('Gene annotation sets differ')
 for p,q in mapping.items():
  for k in ('lower','upper'):
   if any(ga[g][p+'_'+k]!=gb[g][q+'_'+k] for g in ga):raise AssertionError('Full AR003 gene bounds differ')
 return {'full_guide_cells':3*len(a),'full_gene_bound_cells':6*len(ga),'disagreements':0}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--accession',required=True);ap.add_argument('--library',required=True);ap.add_argument('--binary',required=True);ap.add_argument('--work',required=True);ap.add_argument('--out',required=True);ap.add_argument('--previous')
 args=ap.parse_args();out=Path(args.out);work=Path(args.work);work.mkdir(parents=True,exist_ok=True)
 if out.exists():raise FileExistsError(out)
 stage=out.with_name(out.name+'.pending');stage.mkdir(parents=True,exist_ok=False)
 try:
  rows=metadata();row=next(r for r in rows if r['run_accession']==args.accession);lib=Path(args.library).resolve();binary=Path(args.binary).resolve()
  if sha(lib)!=LIB_SHA or len(table(lib,'id'))!=77441:raise AssertionError('Wrong Brunello reference')
  totals={k:sum(int(r['read_count']) for r in rows if r['biological_sample']==k) for k in RUNS}
  rate=1.0 if row['biological_sample']=='pDNA' else min(totals[k] for k in 'ABC')/totals[row['biological_sample']]
  seed=int.from_bytes(hashlib.sha256(('AR004:'+args.accession).encode()).digest()[:8],'big')
  reads=work/(args.accession+'.fastq.gz');url=download(row,reads);raw_sha=sha(reads)
  started=time.monotonic();cmd=[str(binary),str(lib),','.join(map(str,OFFSETS)),str(reads.resolve()),str(stage/'core'),str(rate),str(seed),'131072']
  with (stage/'execution.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
  metrics=json.loads((stage/'core/metrics.json').read_text())
  if metrics['records']!=int(row['read_count']):raise AssertionError('Full FASTQ count differs from ENA')
  valid=verify_tables(stage/'core',lib);valid.update(validate_public(binary,lib,reads,metrics['records']))
  if args.previous:valid['AR003_full_plasmid_replay']=baseline_check(stage/'core',args.previous)
  if sha(reads)!=raw_sha or sha(lib)!=LIB_SHA:raise AssertionError('Input changed during execution')
  # ReCo gets a native sample sheet and unchanged original FASTQ, not extracted windows.
  with (work/'reco-library.tsv').open('w') as f:
   w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerows((r['id'],r['sequence']) for r in table(lib,'id').values())
  header=['Sample name','Sample type','Vector','FastQ read 1','FastQ read 2','Lib 1','Lib 2','Expected reads','Emails','Notes']
  with (work/'reco-samples.tsv').open('w') as f:
   w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(header);w.writerow([args.accession,'single','',str(reads.resolve()),'',str((work/'reco-library.tsv').resolve()),'',metrics['records'],'',''])
  completion={'schema':'dotmatch.ar004.run.v1','completion':'complete','accession':args.accession,'sample':row['biological_sample'],'metadata':row,'raw_url':url,'raw_sha256':raw_sha,'library_sha256':LIB_SHA,'offsets':OFFSETS,'metrics':metrics,'validation':valid,'thinning':{'rate':rate,'seed':seed,'hash':'splitmix64(seed XOR 1-based ordinal); upper 53 bits / 2^53','target_expected_replicate_depth':min(totals[k] for k in 'ABC')},'source':{'runner_sha256':sha(__file__),'counter_sha256':sha(binary),'code_commit':os.environ.get('GITHUB_SHA'),'python':sys.version,'platform':platform.platform()},'command':cmd,'wall_seconds':time.monotonic()-started,'files':{str(p.relative_to(stage)):sha(p) for p in sorted(stage.rglob('*')) if p.is_file()}}
  dump(stage/'completion.json',completion);stage.rename(out);print(args.accession,metrics,valid,flush=True)
 except Exception:
  dump(stage/'FAILED.json',{'completion':'failed','traceback':traceback.format_exc()});raise
if __name__=='__main__':main()
