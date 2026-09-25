#!/usr/bin/env python3
"""Resolve tool contracts and library before outcome inspection."""
import argparse,pathlib,urllib.request,json,subprocess,zipfile,time
p=argparse.ArgumentParser();p.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
urls=['https://sourceforge.net/projects/mageck/files/libraries/yusa_library.csv.zip/download','https://downloads.sourceforge.net/project/mageck/libraries/yusa_library.csv.zip']
for url in urls:
 try:
  with urllib.request.urlopen(url,timeout=90) as r: raw=r.read()
  zpath=a.out/'yusa_library.csv.zip';zpath.write_bytes(raw)
  with zipfile.ZipFile(zpath) as z: data=z.read('yusa_library.csv')
  (a.out/'yusa_library.csv').write_bytes(data); print('LIBRARY',data[:1000],len(data),flush=True);break
 except Exception as e: print(repr(e),flush=True)
else: raise RuntimeError('No original library')
for i,cmd in enumerate([['make','dotmatch'],['./dotmatch','count','--help'],['python','-m','pip','install','.']]):
 r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True); print('COMMAND',cmd,'RC',r.returncode,r.stdout[-20000:],flush=True)
 (a.out/('tool-probe-'+str(i)+'.txt')).write_text(r.stdout)
