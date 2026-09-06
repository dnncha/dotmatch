#!/usr/bin/env python3
"""AR004 supplementary baseline ascertainment, not a changed primary cohort."""
from __future__ import annotations
import argparse,csv,importlib.util,json,math,os
from collections import defaultdict
from pathlib import Path
spec=importlib.util.spec_from_file_location('ar004_stats',Path(__file__).with_name('statistics.py'))
st=importlib.util.module_from_spec(spec);spec.loader.exec_module(st)
def median(values):
 s=sorted(values);n=len(s)
 if not n:raise ValueError('Empty median')
 return s[n//2] if n%2 else (s[n//2-1]+s[n//2])/2

def verify_gene_lfc(folder):
 """Independent reconstruction from rounded original per-guide outputs.
 MAGeCK prints five significant digits; tolerance covers both output roundings.
 This verifies gene aggregation, not its statistical calibration.
 """
 checks=[]
 for gene_file in sorted(Path(folder).glob('*.gene_summary.txt')):
  prefix=gene_file.name.removesuffix('.gene_summary.txt')
  sg=st.read(gene_file.with_name(prefix+'.sgrna_summary.txt'),'sgrna')
  grouped=defaultdict(list)
  for row in sg.values():grouped[row['Gene']].append(st.floatval(row,'LFC'))
  genes=st.gene_summary(gene_file);maximum=0.0
  for gene,row in genes.items():
   if gene not in grouped:raise AssertionError('Missing guide support for gene '+gene)
   observed=st.floatval(row,'neg|lfc');expected=median(grouped[gene]);error=abs(observed-expected)
   tolerance=1.1e-4*max(1.0,abs(expected),max(abs(x) for x in grouped[gene]))
   if error>tolerance:raise AssertionError(f'Independent gene median differs: {prefix} {gene} {observed} {expected}')
   if len(grouped[gene])!=int(row['num']):raise AssertionError('Gene guide number differs from original guide output')
   maximum=max(maximum,error)
  checks.append({'analysis':prefix,'gene_medians_checked':len(genes),'guide_rows':len(sg),'max_absolute_rounding_difference':maximum,'tolerance':'1.1e-4 * max(1, |median|, maximum |guide LFC|); original outputs use five significant digits','disagreements':0})
 return checks

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--primary',required=True);p.add_argument('--out',required=True);a=p.parse_args();base=Path(a.primary);out=Path(a.out)
 if out.exists():raise FileExistsError(out)
 out.mkdir(parents=True)
 completion=json.loads((base/'completion.json').read_text());assert completion['completion']=='complete'
 for name,h in completion['files'].items():
  f=base/name
  if not f.resolve().is_relative_to(base.resolve()) or st.sha(f)!=h:raise AssertionError('Primary result integrity failed')
 data=st.read(base/'all-aggregated-guide-counts.tsv');elig=st.read(base/'eligibility.tsv')
 if set(data)!=set(elig):raise AssertionError('Guide identities differ')
 selected=sorted(i for i,r in data.items() if r['gene'] not in st.EXCLUDED and st.natural(r['pDNA_event_exact'])>=30)
 original={i for i,r in elig.items() if r['decision']=='retained'};extra=sorted(set(selected)-original)
 if (len(original),len(selected),len(extra))!=(74293,74359,66) or original-set(selected):raise AssertionError('Previously declared baseline populations differ')
 st.tsv(out/'additional-baseline-guides.tsv',['id','gene','event_exact_plasmid','joint_exact_plasmid'],((i,data[i]['gene'],data[i]['pDNA_event_exact'],data[i]['pDNA_joint_exact']) for i in extra))
 commands={}
 for policy in st.PRIMARY:
  table=out/(policy+'.tsv')
  st.tsv(table,['sgRNA','Gene',*st.ORDER],((i,data[i]['gene'],*(st.natural(data[i][s+'_'+policy]) for s in st.ORDER)) for i in selected))
  commands[policy]=st.execute_mageck(table,out/'mageck'/policy,'')
 left=st.gene_summary(out/'mageck/event_best.gene_summary.txt');right=st.gene_summary(out/'mageck/joint_best.gene_summary.txt')
 summary,details=st.compare(left,right,'event_best','joint_best','supplementary_event_exact_baseline')
 header=['scope','left','right','gene','left_guides','right_guides','left_lfc','right_lfc','delta_lfc','left_fdr','right_fdr','left_rank','right_rank','left_call_0_05','right_call_0_05','left_distance_from_0_05','right_distance_from_0_05']
 st.tsv(out/'all-gene-comparisons.tsv',header,details);st.tsv(out/'all-discordant-calls.tsv',header,(r for r in details if r[13]!=r[14]))
 validation={'primary':verify_gene_lfc(base/'mageck'),'supplementary':verify_gene_lfc(out/'mageck')}
 h=next(x for x in summary['thresholds'] if x['fdr']==.05)
 text='# Baseline-identifiability sensitivity\n\nThis is the declared supplementary analysis, not a replacement of the primary population. The original 74,293-guide population conditions on >=30 uniquely assigned joint-exact plasmid counts. Using >=30 event-exact plasmid counts admits 66 additional guide IDs from 61 gene annotation labels, for 74,359 guides. Both arms use exactly that same broader baseline-only population. No cellular outcome enters selection.\n\n'
 text+=f"With unchanged MAGeCK settings, event counting called **{h['left_calls']}** negative-selection gene annotations at FDR<=0.05 and joint best-distance counting called **{h['right_calls']}**: **{h['both']}** common, **{h['only_left']}** event-only, **{h['only_right']}** joint-only. Among **{summary['common_genes']}** shared reported genes, **{summary['genes_abs_lfc_difference_ge_0_5']}** had absolute effect differences >=0.5 log2. Negative-selection rank correlation was **{summary['spearman_negative_rank']:.6f}**.\n\n"
 text+='This quantifies conditioning on baseline identifiability. A different gene-call set is not a label of correct or incorrect biology. Both the primary and supplementary universes remain available, with every newly admitted guide identified. Missing statistical gene outputs are reported explicitly rather than zero-filled.\n\n'
 text+='The original per-guide LFC outputs were independently regrouped and their medians compared with all original gene-level LFCs, allowing only for the five-significant-digit formatting precision. This validates the stated median aggregation, not the FDR model or molecular origins. Detailed checks and tolerances are in summary.json.\n'
 (out/'REPORT.md').write_text(text)
 st.dump(out/'summary.json',{'completion':'complete','comparison':summary,'original_guides':len(original),'supplementary_guides':len(selected),'additional_guides':len(extra),'additional_original_gene_labels':len({data[i]['gene'] for i in extra}),'commands':commands,'independent_gene_LFC_checks':validation,'declaration_commit':'a60ad8f06dde2235d555a26d23d45c5dc4d25e69','source_sha256':st.sha(__file__),'source_commit':os.environ.get('GITHUB_SHA')})
 st.dump(out/'completion.json',{'completion':'complete','files':{str(f.relative_to(out)):st.sha(f) for f in sorted(out.rglob('*')) if f.is_file()}})
 print(json.dumps(summary,sort_keys=True),flush=True)
if __name__=='__main__':main()
