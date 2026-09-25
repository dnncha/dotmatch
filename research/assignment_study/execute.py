#!/usr/bin/env python3
"""Network adapter for the frozen study. No assignment/statistical rules change.
SourceForge sometimes denies a regional redirect; public mirrors are accepted
only when the extracted original library has the frozen SHA256.
"""
import hashlib, pathlib, urllib.request, zipfile
import run
original_retrieve=run.retrieve
LIBRARY_URL='https://sourceforge.net/projects/mageck/files/libraries/yusa_library.csv.zip/download'
MIRRORS=[
 'https://downloads.sourceforge.net/project/mageck/libraries/yusa_library.csv.zip',
 'https://master.dl.sourceforge.net/project/mageck/libraries/yusa_library.csv.zip',
 LIBRARY_URL,
 'https://netix.dl.sourceforge.net/project/mageck/libraries/yusa_library.csv.zip',
 'https://phoenixnap.dl.sourceforge.net/project/mageck/libraries/yusa_library.csv.zip']
def retrieve(url,p,expected_md5=None,expected_size=None):
 if url!=LIBRARY_URL:return original_retrieve(url,p,expected_md5,expected_size)
 p=pathlib.Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 failures=[]
 for mirror in MIRRORS:
  try:
   with urllib.request.urlopen(mirror,timeout=90) as r:data=r.read()
   tmp=p.with_name(p.name+'.partial');tmp.write_bytes(data)
   with zipfile.ZipFile(tmp) as z:payload=z.read('yusa_library.csv')
   if hashlib.sha256(payload).hexdigest()!=run.LIB_SHA:raise ValueError('library checksum changed')
   tmp.replace(p);print('VERIFIED LIBRARY MIRROR',mirror,flush=True);return
  except Exception as e:failures.append({'url':mirror,'error':repr(e)});print('LIBRARY MIRROR FAILURE',mirror,repr(e),flush=True)
 raise RuntimeError('All public library mirrors failed: '+repr(failures))
run.retrieve=retrieve
if __name__=='__main__':run.main()
