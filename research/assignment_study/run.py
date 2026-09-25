#!/usr/bin/env python3
"""Metadata-only design resolution; not an outcome analysis."""
import argparse, pathlib, urllib.request, json, zipfile, io
p=argparse.ArgumentParser();p.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
url='https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB4038&result=read_run&fields=run_accession,sample_accession,experiment_accession,sample_title,library_name,read_count,fastq_ftp,fastq_bytes,fastq_md5&format=json'
with urllib.request.urlopen(url,timeout=90) as r: raw=r.read()
(a.out/'study_metadata.json').write_bytes(raw)
rows=json.loads(raw)
for row in rows: print(json.dumps(row),flush=True)
for row in rows:
 if row['run_accession'] in ['ERR376998','ERR376999','ERR377000','ERR377001']:
  for key in ['sample_accession','experiment_accession']:
   acc=row[key]
   with urllib.request.urlopen('https://www.ebi.ac.uk/ena/browser/api/xml/'+acc,timeout=90) as r: raw=r.read()
   (a.out/(acc+'.xml')).write_bytes(raw)
   print(raw.decode(),flush=True)
