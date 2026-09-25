#!/usr/bin/env python3
"""Recompute the principal published descriptive results with Python stdlib.
Input: the unmodified output directory from primary run 36191069868, or a
fresh successful execution of the frozen study. Does not download raw reads.
No causal, molecular-origin or statistical-significance interpretation.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,pathlib,statistics
from collections import defaultdict

POLICIES=['fixed_exact','fixed_radius1','fixed_best1','multi_exact','distinct_exact','multi_best1','distinct_best1','joint_best1']
COMPARISONS={'deduplicate_best':('distinct_best1','multi_best1'),'deduplicate_exact':('distinct_exact','multi_exact'),'fixed_mismatch':('fixed_best1','fixed_exact'),'radius_vs_best':('fixed_radius1','fixed_best1'),'multi_vs_fixed':('multi_best1','fixed_best1'),'joint_vs_multi':('joint_best1','multi_best1')}
EXPECTED={'deduplicate_best':[0,0],'deduplicate_exact':[0,0],'fixed_mismatch':[178,199],'radius_vs_best':[43,41],'multi_vs_fixed':[324,379],'joint_vs_multi':[489,547]}
def rows(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--primary',type=pathlib.Path,required=True);p.add_argument('--out',type=pathlib.Path);a=p.parse_args();root=a.primary
 for line in (root/'SHA256SUMS').read_text().splitlines():
  h,path=line.split('  ',1)
  if sha(root/path)!=h:raise RuntimeError('Checksum failure: '+path)
 complete=json.loads((root/'COMPLETE.json').read_text());assert complete['discrepancies']==0 and complete['total_reads']==62831612
 manifest=json.loads((root/'input_manifest.json').read_text());first=rows(root/'oracle'/(manifest[0]['accession']+'.counts.tsv'))
 ids=[r['guide'] for r in first];genes=[r['gene'] for r in first];assert len(ids)==87437
 ci={g:i for i,g in enumerate(ids)};samples=['plasmid','ESC1','ESC2'];counts={p:[[0]*len(ids) for _ in samples] for p in POLICIES}
 qc=defaultdict(int);same_guide_differences=0
 for m in manifest:
  acc=m['accession'];sample=samples.index(m['sample']);rs=rows(root/'oracle'/(acc+'.counts.tsv'));assert {r['guide'] for r in rs}==set(ids)
  for r in rs:
   i=ci[r['guide']]
   for policy in POLICIES:counts[policy][sample][i]+=int(r[policy])
   same_guide_differences+=int(r['distinct_exact'])!=int(r['multi_exact'])
   same_guide_differences+=int(r['distinct_best1'])!=int(r['multi_best1'])
  for r in rows(root/'oracle'/(acc+'.qc.tsv')):qc[r['mode']+'_'+r['metric']]+=int(r['value'])
 assert same_guide_differences==0
 x=counts['fixed_exact'];positive=[i for i in range(len(ids)) if all(x[s][i]>0 for s in range(3))]
 gm={i:math.exp(sum(math.log(x[s][i]) for s in range(3))/3) for i in positive}
 sf=[statistics.median(x[s][i]/gm[i] for i in positive) for s in range(3)];scale=math.exp(sum(map(math.log,sf))/3);sf=[v/scale for v in sf]
 eligible=[i for i in range(len(ids)) if x[0][i]>=60];bygene=defaultdict(list)
 for i in eligible:bygene[genes[i]].append(i)
 bygene={g:is_ for g,is_ in bygene.items() if len(is_)>=3};assert len(eligible)==73567 and len(bygene)==16654
 effects={}
 for policy in POLICIES:
  z=counts[policy];effects[policy]={}
  for g,indices in bygene.items():effects[policy][g]=[statistics.median(math.log2((z[s][i]/sf[s]+1)/(z[0][i]/sf[0]+1)) for i in indices) for s in [1,2]]
 observed={name:[sum(abs(effects[left][g][s]-effects[right][g][s])>=.5 for g in bygene) for s in [0,1]] for name,(left,right) in COMPARISONS.items()}
 assert observed==EXPECTED,(observed,EXPECTED)
 assert qc['multi_best1_multiple_hit_reads']==4371978
 assert qc['multi_best1_same_gene_cross_guide_reads']==4324884
 assert qc['multi_best1_cross_gene_reads']==47094
 assert qc['multi_best1_window_contributions']==61871230 and qc['multi_best1_matched_reads']==57278925
 assert qc['fixed_best1_unique']-qc['fixed_exact_unique']==2642177
 report=dict(status='pass',total_reads=62831612,library_targets=len(ids),eligible_guides=len(eligible),eligible_genes=len(bygene),normalization='common fixed-exact median-ratio, geometric recentering; pseudocount1; gene median of eligible guide log2 ratios',same_guide_count_differences=same_guide_differences,multi_hit_reads=4371978,same_gene_multi_hit_reads=4324884,cross_gene_multi_hit_reads=47094,flag_counts_ESC1_ESC2=observed,cases={g:{p:effects[p][g] for p in ['fixed_exact','fixed_best1']} for g in ['RPL4','EMG1']},scope='Descriptive measurement sensitivity in three biosamples resequenced twice; no FDR or novel dependency claim')
 text=json.dumps(report,indent=2)+'\n'
 if a.out:a.out.write_text(text)
 print(text)
if __name__=='__main__':main()
