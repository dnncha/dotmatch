#!/usr/bin/env python3
"""Separate native-workflow comparison. Never treats incomplete ReCo as counts."""
import argparse,csv,importlib.util,json,os,traceback
from collections import Counter
from pathlib import Path
spec=importlib.util.spec_from_file_location('ar004_statistics',Path(__file__).with_name('statistics.py'))
st=importlib.util.module_from_spec(spec);spec.loader.exec_module(st)

def main():
 p=argparse.ArgumentParser();p.add_argument('--statistics',required=True);p.add_argument('--reco',required=True);p.add_argument('--out',required=True);args=p.parse_args()
 out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
 base=Path(args.statistics);completion=json.loads((base/'completion.json').read_text())
 if completion['completion']!='complete':raise ValueError('Incomplete primary statistics')
 for name,h in completion['files'].items():
  item=base/name
  if not item.resolve().is_relative_to(base.resolve()) or st.sha(item)!=h:raise AssertionError('Primary result checksum mismatch')
 info=st.read(base/'eligibility.tsv');eligible=sorted(i for i,r in info.items() if r['decision']=='retained')
 expected={a for group in st.RUNS.values() for a in group};found={};failures=[]
 for file in Path(args.reco).rglob('status.json'):
  m=json.loads(file.read_text());a=m.get('accession')
  if a not in expected:continue
  if a in found:raise ValueError('Duplicate comparator archive')
  found[a]=(file.parent,m)
  if m.get('completion')!='complete':failures.append({'accession':a,'failure':m.get('traceback','No completion')})
 missing=sorted(expected-set(found))
 if failures or missing:
  st.dump(out/'status.json',{'completion':'comparator_incomplete','failed_runs':failures,'missing_runs':missing,'substituted_counts':False})
  (out/'REPORT.md').write_text('# ReCo comparison status\n\nThe actual comparator did not complete every required original archive. No partial cohort, fabricated zero counts or emulated results were substituted. Full failures are retained in status.json. The separate matched-position primary analysis is unaffected.\n')
  return
 primary_summary=json.loads((base/'summary.json').read_text());matrices={};run_rows=[]
 for sample,accessions in st.RUNS.items():
  values=Counter()
  for a in accessions:
   folder,m=found[a]
   if m['library_sha256']!=st.LIB_SHA or m['source_commit']!='e2daf48b610f8db29ad014bff5be8bb983aaa76f':raise AssertionError('Wrong actual comparator source/reference')
   for name,h in m['files'].items():
    if st.sha(folder/name)!=h:raise AssertionError('Comparator artifact hash differs')
   with (folder/'original-counts.csv').open(newline='') as f:
    reader=csv.reader(f);head=next(reader)
    if len(head)!=2 or head[0]!='Guide':raise ValueError('Unexpected actual ReCo count schema')
    counts={}
    for row in reader:
     if len(row)!=2 or row[0] in counts:raise ValueError('Duplicate or malformed comparator row')
     counts[row[0]]=st.natural(row[1])
   if set(counts)!=set(info):raise AssertionError('ReCo omitted or added reference guides; no implicit zero-filling')
   n=int(m['metadata']['read_count'])
   if n!=primary_summary['run_manifests'][a]['metrics']['records']:raise AssertionError('Original archive record identities differ')
   raw_table=st.read(base/'run-accounting.tsv','run')
   if raw_table[a]['raw_sha256']!=m['raw_sha256']:raise AssertionError('Actual raw archive bytes differ')
   if sum(counts.values())>n:raise AssertionError('Native comparator count total exceeds raw records; requires investigation')
   values.update(counts);run_rows.append((a,sample,n,sum(counts.values()),m['raw_sha256']))
  matrices[sample]=values
 st.tsv(out/'native-reco-counts.tsv',['sgRNA','Gene',*st.ORDER],((i,info[i]['gene'],*(matrices[s][i] for s in st.ORDER)) for i in eligible))
 st.tsv(out/'run-accounting.tsv',['run','biological_sample','raw_records','ReCo_assigned_count','raw_sha256'],run_rows)
 command=st.execute_mageck(out/'native-reco-counts.tsv',out/'mageck'/'native_reco','')
 repeat=st.execute_mageck(out/'native-reco-counts.tsv',out/'mageck'/'native_reco_repeat','')
 reco=st.gene_summary(out/'mageck/native_reco.gene_summary.txt')
 comparisons=[];detail=[]
 for policy in ('joint_best','event_best','joint_exact'):
  ref=st.gene_summary(base/'mageck'/('full_'+policy+'.gene_summary.txt'))
  s,d=st.compare(ref,reco,policy,'actual_ReCo','native_workflow');comparisons.append(s);detail+=d
 s,d=st.compare(reco,st.gene_summary(out/'mageck/native_reco_repeat.gene_summary.txt'),'actual_ReCo','actual_ReCo','identical_input_repeat');comparisons.append(s);detail+=d
 header=['scope','left','right','gene','left_guides','right_guides','left_lfc','right_lfc','delta_lfc','left_fdr','right_fdr','left_rank','right_rank','left_call_0_05','right_call_0_05','left_distance_from_0_05','right_distance_from_0_05']
 st.tsv(out/'all-gene-comparisons.tsv',header,detail);st.tsv(out/'all-discordant-calls.tsv',header,(r for r in detail if r[13]!=r[14]))
 text='# Actual ReCo native-workflow comparison\n\nThe unchanged ReCo source at `e2daf48b610f8db29ad014bff5be8bb983aaa76f`, Cutadapt 2.8 and Bowtie 2.3.0 processed all seven complete original archives. This comparison changes extraction/context rules and assignment; it does not isolate one algorithmic choice or establish true biological accuracy.\n\nThe same baseline-only guide population and MAGeCK 0.5.9.5 settings used in the matched-position study were applied to its counts. Technical files were combined within the three biological replicates.\n\n'
 text+='| Comparison | Common genes | Effect changes >=0.5 | Only left at FDR0.05 | Only ReCo at FDR0.05 | Rank correlation |\n|---|---:|---:|---:|---:|---:|\n'
 for s in comparisons:
  h=next(r for r in s['thresholds'] if r['fdr']==.05)
  text+=f"| {s['scope']}: {s['left']} vs {s['right']} | {s['common_genes']} | {s['genes_abs_lfc_difference_ge_0_5']} | {h['only_left']} | {h['only_right']} | {s['spearman_negative_rank']:.6f} |\n"
 text+='\nAll gene rows, FDR0.01/0.05/0.10 comparisons, individual archive count budgets, original commands and repeated-input checks are retained. Discordant calls are not labelled true/false discoveries. No current-version replacement or handwritten ReCo imitation was used. Original release acquisition/dependency failures are retained in their earlier workflow logs and do not enter these results.\n'
 (out/'REPORT.md').write_text(text)
 st.dump(out/'summary.json',{'completion':'complete','kind':'actual_native_workflow_comparator','raw_records':sum(r[2] for r in run_rows),'eligible_guides':len(eligible),'comparisons':comparisons,'command':command,'repeat_command':repeat,'source_sha256':st.sha(__file__),'source_commit':os.environ.get('GITHUB_SHA'),'known_origin_accuracy_established':False})
 st.dump(out/'completion.json',{'completion':'complete','files':{str(f.relative_to(out)):st.sha(f) for f in sorted(out.rglob('*')) if f.is_file()}})
 print(json.dumps(comparisons,sort_keys=True),flush=True)
if __name__=='__main__':main()
