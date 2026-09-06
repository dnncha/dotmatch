#!/usr/bin/env python3
"""Run actual pinned ReCo on complete raw archives; do not emulate failed runs."""
import argparse,csv,json,os,shutil,subprocess,time,traceback
from pathlib import Path
import run as ar4

def main():
 p=argparse.ArgumentParser();p.add_argument('accession');p.add_argument('--library',required=True);p.add_argument('--source',required=True);a=p.parse_args()
 out=Path('reco-result')/a.accession;out.mkdir(parents=True,exist_ok=False);work=Path('reco-work');work.mkdir()
 status={'accession':a.accession,'completion':'failed','kind':'actual_ReCo_native_workflow','source_commit':'e2daf48b610f8db29ad014bff5be8bb983aaa76f','same_window_comparison':False}
 try:
  lib=Path(a.library).resolve()
  if ar4.sha(lib)!=ar4.LIB_SHA:raise ValueError('Reference digest mismatch')
  source=Path(a.source).resolve();commit=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
  if commit!=status['source_commit']:raise ValueError('Comparator source differs')
  rows=ar4.metadata();row=next(r for r in rows if r['run_accession']==a.accession);raw=work/(a.accession+'.fastq.gz');url=ar4.download(row,raw)
  with (work/'library.tsv').open('w') as f:
   w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerows((r['id'],r['sequence']) for r in ar4.table(lib,'id').values())
  header=['Sample name','Sample type','Vector','FastQ read 1','FastQ read 2','Lib 1','Lib 2','Expected reads','Emails','Notes']
  with (work/'samples.tsv').open('w') as f:
   w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(header);w.writerow([a.accession,'single','',str(raw.resolve()),'',str((work/'library.tsv').resolve()),'',int(row['read_count']),'',''])
  command=['micromamba','run','-n','reco','python','-c','import reco,sys; r=reco.ReCo(sample_sheet_file=sys.argv[1],output_dir=sys.argv[2]); r.run(remove_unused_files=True,cores=2)',str((work/'samples.tsv').resolve()),str((work/'output').resolve())]
  start=time.monotonic()
  with (out/'execution.log').open('w') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=1500)
  final=list((work/'output').rglob('*_final_guidecounts.csv'))
  if len(final)!=1:raise ValueError('Actual ReCo did not produce exactly one count table')
  shutil.copy2(final[0],out/'original-counts.csv')
  for i,log in enumerate(sorted((work/'output').rglob('*.log'))):shutil.copy2(log,out/('upstream-'+str(i)+'.log'))
  for i,report in enumerate(sorted((work/'output').rglob('report.txt'))):shutil.copy2(report,out/('upstream-report-'+str(i)+'.txt'))
  with final[0].open() as f:
   table=list(csv.reader(f));status['output_columns']=table[0];status['output_rows']=len(table)-1
  status.update(completion='complete',metadata=row,raw_url=url,raw_sha256=ar4.sha(raw),library_sha256=ar4.LIB_SHA,command=command,wall_seconds=time.monotonic()-start)
 except Exception:status['traceback']=traceback.format_exc();print(status['traceback'],flush=True)
 status['files']={f.name:ar4.sha(f) for f in out.iterdir() if f.is_file()}
 ar4.dump(out/'status.json',status);print(json.dumps(status,sort_keys=True),flush=True)
if __name__=='__main__':main()
