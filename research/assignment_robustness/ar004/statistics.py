#!/usr/bin/env python3
"""AR004: verified complete-run aggregation and pinned MAGeCK RRA sensitivity.
One shared plasmid, three biological replicates; not biological ground truth.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,os,platform,shutil,sqlite3,subprocess,sys,time,unittest
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
POLICIES=('event_exact','event_best','joint_exact','joint_best','joint_radius')
RUNS={'pDNA':['SRR8297997'],'A':['SRR8297836','SRR8297837'],'B':['SRR8297838','SRR8297839'],'C':['SRR8297840','SRR8297841']}
ORDER=('pDNA','A','B','C')
EXCLUDED=frozenset(('Non-Targeting-Control','NA','na',''))
LIB_SHA='0d2906187829ea9f736de94a47369bd94d42cde5f348fea9d12a385625cc2ca1'
COUNT_COMMIT='862bdc1e0c026b2ce7475870af1a4a2171c98cfb'
PRIMARY=('event_best','joint_best')
COMPARISONS=(PRIMARY,('event_exact','joint_exact'),('joint_exact','joint_best'),('joint_best','joint_radius'))
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def dump(path,value):
 Path(path).parent.mkdir(parents=True,exist_ok=True)
 Path(path).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')
def tsv(path,header,rows):
 with Path(path).open('w',newline='',encoding='utf-8') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(header);w.writerows(rows)
def read(path,key='id'):
 with Path(path).open(newline='',encoding='utf-8-sig') as f:
  reader=csv.DictReader(f,delimiter='\t')
  if not reader.fieldnames or len(set(reader.fieldnames))!=len(reader.fieldnames):raise ValueError('Invalid header')
  rows=list(reader)
 result={}
 for row in rows:
  if None in row or any(v is None for v in row.values()) or row[key] in result:raise ValueError('Malformed or duplicate row')
  result[row[key]]=row
 if not result:raise ValueError('Empty table')
 return result
def natural(value):
 if not isinstance(value,str) or not value.isascii() or not value.isdecimal():raise ValueError('Invalid integer: '+repr(value))
 return int(value)
def load_evidence(root):
 found={};expected={a for group in RUNS.values() for a in group}
 for path in sorted(Path(root).rglob('completion.json')):
  m=json.loads(path.read_text());a=m.get('accession')
  if a not in expected:continue
  if a in found or m.get('completion')!='complete':raise ValueError('Repeated or incomplete archive '+a)
  if m['library_sha256']!=LIB_SHA or tuple(m['offsets'])!=tuple(range(21,31)):raise ValueError('Reference or search-domain mismatch')
  if m['source']['code_commit']!=COUNT_COMMIT:raise ValueError('Unexpected raw counting source')
  if m['metrics']['records']!=int(m['metadata']['read_count']) or m['validation']['disagreements']!=0:raise AssertionError('Input/validation failed')
  for name,checksum in m['files'].items():
   item=path.parent/name
   if not item.resolve().is_relative_to(path.parent.resolve()) or sha(item)!=checksum:raise AssertionError('Artifact checksum mismatch '+str(item))
  group=next(k for k,v in RUNS.items() if a in v)
  if m['sample']!=group:raise ValueError('Run/replicate assignment differs')
  found[a]=(path.parent,m,sha(path))
 if set(found)!=expected:raise ValueError('All seven complete original runs are required; missing '+str(expected-set(found)))
 for group,accessions in RUNS.items():
  if len({found[a][1]['metadata']['sample_accession'] for a in accessions})!=1:raise ValueError('Technical run sample IDs differ')
 if len({found[a][1]['raw_sha256'] for a in found})!=7:raise ValueError('Duplicate raw files cannot create replication')
 return found

def aggregate(evidence,out):
 matrices={};identity=None;exposure={};run_rows=[];guide_cells=0
 db=sqlite3.connect(':memory:')
 db.execute('CREATE TABLE raw(sample TEXT,guide TEXT,gene TEXT,policy TEXT,n INTEGER NOT NULL,PRIMARY KEY(sample,guide,policy))')
 for sample,accessions in RUNS.items():
  group={};exposure[sample]=dict(full=0,thin=0)
  for accession in accessions:
   folder,m,_=evidence[accession];rows=read(folder/'core/guide-counts.tsv');ident={i:r['gene'] for i,r in rows.items()}
   if identity is None:identity=ident
   if ident!=identity:raise AssertionError('All archive guide/gene identities must match exactly')
   for policy in (*POLICIES,*("thin_"+p for p in POLICIES)):
    values={i:natural(r[policy]) for i,r in rows.items()}
    if policy not in group:group[policy]=Counter()
    group[policy].update(values)
    db.executemany('INSERT INTO raw VALUES(?,?,?,?,?) ON CONFLICT(sample,guide,policy) DO UPDATE SET n=n+excluded.n',((sample,i,identity[i],policy,n) for i,n in values.items()))
    guide_cells+=len(values)
   exposure[sample]['full']+=m['metrics']['records'];exposure[sample]['thin']+=m['metrics']['thinned_records']
   run_rows.append((accession,sample,m['metadata']['sample_accession'],m['metrics']['records'],m['metrics']['thinned_records'],m['metrics']['policies']['event_best']['matched_records'],m['metrics']['policies']['event_best']['count_events'],m['metrics']['policies']['event_best']['extra_events'],m['raw_sha256']))
  matrices[sample]=group
 for sample,i,p,n in db.execute('SELECT sample,guide,policy,n FROM raw'):
  if matrices[sample][p][i]!=n:raise AssertionError('Independent SQL technical-run aggregation differs')
 if len(identity)!=77441:raise AssertionError('Reference population differs')
 sql_genes={(s,g,p):n for s,g,p,n in db.execute('SELECT sample,gene,policy,SUM(n) FROM raw GROUP BY sample,gene,policy')};db.close()
 if sum(evidence[a][1]['metrics']['records'] for a in evidence)!=246950411:raise AssertionError('Locked cohort size differs')
 eligible=[];exclusions=[]
 for i,gene in sorted(identity.items()):
  base=matrices['pDNA']['joint_exact'][i]
  reason='annotation_control_or_missing' if gene in EXCLUDED else 'baseline_below_30' if base<30 else 'retained'
  if reason=='retained':eligible.append(i)
  exclusions.append((i,gene,base,reason))
 if not eligible:raise AssertionError('No eligible guides')
 tsv(out/'eligibility.tsv',['id','gene','joint_exact_plasmid','decision'],exclusions)
 tsv(out/'run-accounting.tsv',['run','sample','sample_accession','records','thinned_records','event_best_matched','event_best_events','event_best_extra','raw_sha256'],run_rows)
 countdir=out/'count-matrices';countdir.mkdir()
 for scope in ('full','thin'):
  for p in POLICIES:
   column=p if scope=='full' else 'thin_'+p
   tsv(countdir/(scope+'_'+p+'.tsv'),['sgRNA','Gene',*ORDER],((i,identity[i],*(matrices[s][column][i] for s in ORDER)) for i in eligible))
 tsv(out/'all-aggregated-guide-counts.tsv',['id','gene',*(s+'_'+p for s in ORDER for p in (*POLICIES,*("thin_"+p for p in POLICIES)))],((i,identity[i],*(matrices[s][p][i] for s in ORDER for p in (*POLICIES,*("thin_"+p for p in POLICIES)))) for i in sorted(identity)))
 genecounts={}
 for s in ORDER:
  for p in (*POLICIES,*("thin_"+p for p in POLICIES)):
   values=Counter()
   for i,n in matrices[s][p].items():values[identity[i]]+=n
   if any(sql_genes[(s,g,p)]!=n for g,n in values.items()):raise AssertionError('Independent SQL/Python gene sums differ')
   genecounts[(s,p)]=values
 return matrices,identity,eligible,exposure,genecounts,dict(raw_guide_count_cells_checked=guide_cells,technical_sample_groups=4,biological_cellular_replicates=3,shared_plasmid_samples=1,SQL_guide_and_gene_aggregation=True)

def floatval(row,key):
 x=float(row[key])
 if not math.isfinite(x):raise ValueError('Nonfinite '+key)
 return x
def gene_summary(path):
 rows=read(path)
 for g,r in rows.items():
  for k in ('neg|fdr','neg|p-value','neg|lfc','neg|score','neg|rank'):floatval(r,k)
  if not 0<=float(r['neg|fdr'])<=1 or not 0<=float(r['neg|p-value'])<=1:raise ValueError('Invalid probability')
 return rows
def execute_mageck(matrix,prefix,extra):
 prefix.parent.mkdir(parents=True,exist_ok=True)
 cmd=['mageck','test','-k',str(matrix),'-t','A,B,C','-c','pDNA','--norm-method','median','--remove-zero','none','--gene-lfc-method','median','--keep-tmp','--normcounts-to-file','-n',str(prefix)]
 if extra:cmd+=['--additional-rra-parameters',extra]
 start=time.monotonic()
 with Path(str(prefix)+'.command.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=900)
 path=Path(str(prefix)+'.gene_summary.txt');rows=gene_summary(path)
 return {'command':cmd,'matrix_sha256':sha(matrix),'gene_summary_sha256':sha(path),'genes':len(rows),'wall_seconds':time.monotonic()-start}
def compare(left,right,leftname,rightname,scope):
 from scipy.stats import spearmanr
 universe=sorted(set(left)&set(right));missing_left=sorted(set(right)-set(left));missing_right=sorted(set(left)-set(right));diffs=[];detail=[]
 for g in universe:
  a,b=left[g],right[g];lfca=floatval(a,'neg|lfc');lfcb=floatval(b,'neg|lfc');fa=floatval(a,'neg|fdr');fb=floatval(b,'neg|fdr');delta=lfcb-lfca;diffs.append(abs(delta))
  detail.append((scope,leftname,rightname,g,natural(a['num']),natural(b['num']),lfca,lfcb,delta,fa,fb,floatval(a,'neg|rank'),floatval(b,'neg|rank'),int(fa<=.05),int(fb<=.05),abs(fa-.05),abs(fb-.05)))
 rho=float(spearmanr([float(left[g]['neg|rank']) for g in universe],[float(right[g]['neg|rank']) for g in universe]).statistic) if len(universe)>1 else None
 if rho is not None and not math.isfinite(rho):rho=None
 ordered=sorted(diffs)
 summary={'scope':scope,'left':leftname,'right':rightname,'common_genes':len(universe),'missing_left':missing_left,'missing_right':missing_right,'median_abs_lfc_difference':ordered[len(ordered)//2] if ordered else None,'max_abs_lfc_difference':max(ordered,default=None),'genes_abs_lfc_difference_ge_0_5':sum(d>=.5 for d in diffs),'spearman_negative_rank':rho,'thresholds':[]}
 for threshold in (.01,.05,.10):
  a={g for g in universe if float(left[g]['neg|fdr'])<=threshold};b={g for g in universe if float(right[g]['neg|fdr'])<=threshold}
  summary['thresholds'].append(dict(fdr=threshold,left_calls=len(a),right_calls=len(b),both=len(a&b),only_left=len(a-b),only_right=len(b-a),neither=len(universe)-len(a|b),jaccard=len(a&b)/len(a|b) if a|b else None))
 return summary,detail

def descriptive(matrices,identity,eligible,exposure):
 allowed=defaultdict(list)
 for i in eligible:allowed[identity[i]].append(i)
 rows=[];summaries=[]
 for left,right in COMPARISONS:
  shifts_all=0;sign_consistent=0
  for gene,ids in sorted(allowed.items()):
   ds=[]
   for s in 'ABC':
    c0=sum(matrices['pDNA'][left][i] for i in ids);t0=sum(matrices[s][left][i] for i in ids);c1=sum(matrices['pDNA'][right][i] for i in ids);t1=sum(matrices[s][right][i] for i in ids)
    a=math.log2((t0+.5)/(c0+.5)*exposure['pDNA']['full']/exposure[s]['full']);b=math.log2((t1+.5)/(c1+.5)*exposure['pDNA']['full']/exposure[s]['full'])
    rows.append((left,right,gene,s,c0,t0,c1,t1,a,b,b-a));ds.append(b-a)
   shifts_all+=all(abs(d)>=.5 for d in ds)
   sign_consistent+=all(d>0 for d in ds) or all(d<0 for d in ds)
  summaries.append(dict(left=left,right=right,genes=len(allowed),all_three_replicates_abs_difference_ge_0_5=shifts_all,all_three_policy_delta_same_sign=sign_consistent))
 return rows,summaries

def make_report(out,summary,comparisons):
 primary=next(r for r in comparisons if r['scope']=='full' and (r['left'],r['right'])==PRIMARY);hit=next(r for r in primary['thresholds'] if r['fdr']==.05)
 text='# Replicated-screen assignment sensitivity: AR004\n\n**Executed technical study, 6 September 2026. Not peer reviewed.**\n\n'
 text+=f"Seven complete public archives supplied **{summary['reads']:,} sequencing records**: one matched modified-tracr plasmid baseline and three A375 dropout biological replicates, each sequenced in two runs. Technical runs were summed within biological samples before statistics.\n\n"
 text+='## Primary matched-position result\n\n'
 text+=f"At MAGeCK negative-selection FDR <=0.05, per-position best-distance event counting called **{hit['left_calls']:,}** gene annotations and joint best-distance counting called **{hit['right_calls']:,}**. **{hit['both']:,}** were common, **{hit['only_left']:,}** were event-only and **{hit['only_right']:,}** were joint-only. These are method-dependent calls, not verified true or false discoveries.\n\n"
 text+=f"Among {primary['common_genes']:,} common reported genes, **{primary['genes_abs_lfc_difference_ge_0_5']:,}** had absolute MAGeCK log2-effect differences >=0.5. Negative-selection rank correlation was **{primary['spearman_negative_rank']:.6f}**. Both counting arms used the same ten allowed positions, one-mismatch model, guide reference, common baseline-filtered guide subset and downstream settings.\n\n"
 text+='## All prespecified contrasts\n\n| Scope | Left | Right | Common genes | Effect changes >=0.5 | Only left at 0.05 | Only right at 0.05 |\n|---|---|---|---:|---:|---:|---:|\n'
 for row in comparisons:
  h=next(t for t in row['thresholds'] if t['fdr']==.05)
  text+=f"| {row['scope']} | {row['left']} | {row['right']} | {row['common_genes']:,} | {row['genes_abs_lfc_difference_ge_0_5']:,} | {h['only_left']:,} | {h['only_right']:,} |\n"
 text+='\nThe thinned analysis uses deterministic outcome-independent Bernoulli selection targeting the minimum expected cellular-replicate exposure. Realized counts are recorded and need not be exactly equal. Identical selected records are used across policies. Repeated-input comparisons measure variation from running the same statistical program again; they are not independent experiments.\n\n'
 text+='## Reproducibility and validation\n\n'
 text+=f"The baseline-only common filter retained **{summary['eligible_guides']:,} guides** spanning **{summary['eligible_genes']:,} original gene labels**. Eligibility requires >=30 joint-exact plasmid counts; non-targeting and absent-annotation groups are explicitly excluded from gene inference. All reference rows remain in original and aggregated count tables, with every exclusion in `eligibility.tsv`.\n\n"
 text+='Input archive byte lengths, MD5, SHA-256 and complete record totals were verified. The optimized C++ counter first reproduced every AR003 Brunello plasmid joint guide count and gene bound. Every new archive was checked with 200 seeded records against all reference targets at every allowed position and against pinned DotMatch 0.5.0. Complete count-state budgets and technical-run/annotation aggregation were verified, including independent SQL sums. Sampled all-target validation is not represented as exhaustive validation of all 247 million records.\n\n'
 text+='MAGeCK 0.5.9.5 was built from the original source. Median normalization, no zero-count removal, median gene log-fold changes and the same A/B/C versus shared plasmid contrast were used for all arms. Full commands, executable/source hashes, environment, original gene outputs, normalized counts and intermediate files are retained. The primary FDR threshold is 0.05; 0.01/0.10 are sensitivity thresholds, not separate experiments.\n\n'
 text+='## Scientific interpretation\n\nThis stage answers a narrower and stronger question than a mapping-rate comparison: do read-accounting choices change statistical conclusions in a replicated screen when the search domain is held fixed? The answers above quantify that sensitivity. They do **not** identify which changed calls are biologically correct, establish false-discovery calibration, or demonstrate general superiority.\n\n'
 text+='The cellular outcomes were unseen before the AR004 protocol, but this remains the same Sanson study and Brunello library as the prior plasmid audit, not independent-study replication. The three biological replicates share one plasmid baseline. Gene annotations are retained as supplied. No copy-number correction, independent gene-function validation or experimentally known molecular origin is added. Event counting is the explicitly implemented per-position rule, not an undeclared emulation presented as actual guide-counter execution. ReCo is a separate native-workflow comparison with different extraction and must be reported separately.\n\n'
 text+='All significant/non-significant and discrepant rows are retained. Genes near an FDR boundary should not be called discoveries merely because one configuration moves them across it. Stable effects alongside changed p-values require inspection of guide evidence, normalization and statistical variability rather than an automatic biological narrative.\n\n'
 text+='## Sources and status\n\nPrimary study: Sanson et al., *Nature Communications* (2018), doi:10.1038/s41467-018-07901-8. Original run metadata: ENA PRJNA508200, archived before outcomes. MAGeCK source/documentation: https://sourceforge.net/projects/mageck/ and https://sourceforge.net/p/mageck/wiki/usage/. Frozen AR004 protocol commit: `dc1971a9a9b3938dea17e6961ff16a515a52cd96`. Complete counting commit: `862bdc1e0c026b2ce7475870af1a4a2171c98cfb`, workflow `34034542345`.\n\nNo production defaults, manuscript submission, author outreach or new biological mechanism is claimed.\n'
 (out/'REPORT.md').write_text(text)

def build(args):
 out=Path(args.out)
 if out.exists():raise FileExistsError(out)
 out.mkdir(parents=True);evidence=load_evidence(args.evidence)
 matrices,identity,eligible,exposure,genes,validation=aggregate(evidence,out)
 version=subprocess.check_output(['mageck','--version'],text=True,stderr=subprocess.STDOUT).strip()
 if '0.5.9.5' not in version:raise ValueError('Wrong downstream version '+version)
 testhelp=subprocess.check_output(['mageck','test','--help'],text=True);(out/'mageck-help.txt').write_text(testhelp)
 rra=subprocess.run(['RRA'],text=True,capture_output=True);rrahelp=rra.stdout+rra.stderr;(out/'rra-help.txt').write_text(rrahelp)
 extra='--seed 20260906' if '--seed' in rrahelp else ''
 tasks=[]
 for scope in ('full','thin'):
  for policy in POLICIES:tasks.append((scope,policy,out/'count-matrices'/(scope+'_'+policy+'.tsv'),out/'mageck'/(scope+'_'+policy)))
 for policy in PRIMARY:tasks.append(('repeat',policy,out/'count-matrices'/('full_'+policy+'.tsv'),out/'mageck'/('repeat_'+policy)))
 commands={}
 def work(task):
  scope,p,matrix,prefix=task
  return scope+'_'+p,execute_mageck(matrix,prefix,extra)
 with ThreadPoolExecutor(max_workers=2) as pool:
  for key,result in pool.map(work,tasks):commands[key]=result;print(key,'completed',result['genes'],'gene rows',flush=True)
 outputs={key:gene_summary(out/'mageck'/(key+'.gene_summary.txt')) for key in commands};summaries=[];details=[]
 for scope in ('full','thin'):
  for left,right in COMPARISONS:
   s,d=compare(outputs[scope+'_'+left],outputs[scope+'_'+right],left,right,scope);summaries.append(s);details+=d
 for policy in PRIMARY:
  s,d=compare(outputs['full_'+policy],outputs['repeat_'+policy],policy,policy,'identical_input_repeat');summaries.append(s);details+=d
 header=['scope','left','right','gene','left_guides','right_guides','left_lfc','right_lfc','delta_lfc','left_fdr','right_fdr','left_rank','right_rank','left_call_0_05','right_call_0_05','left_distance_from_0_05','right_distance_from_0_05']
 tsv(out/'all-gene-comparisons.tsv',header,details)
 tsv(out/'all-discordant-calls.tsv',header,(r for r in details if r[13]!=r[14]))
 tsv(out/'all-large-effect-differences.tsv',header,(r for r in details if abs(r[8])>=.5))
 desc,descsum=descriptive(matrices,identity,eligible,exposure)
 tsv(out/'replicate-effect-sensitivity.tsv',['left','right','gene','biological_replicate','left_pDNA','left_cellular','right_pDNA','right_cellular','left_log2','right_log2','delta_log2'],desc)
 sums={'completion':'complete','schema':'dotmatch.ar004.statistics.v1','reads':sum(v[1]['metrics']['records'] for v in evidence.values()),'biological_cellular_replicates':3,'technical_run_archives':7,'shared_plasmid_baselines':1,'eligible_guides':len(eligible),'eligible_genes':len({identity[i] for i in eligible}),'excluded_annotation_groups':sorted(EXCLUDED),'exposures':exposure,'run_manifests':{a:{'completion_sha256':v[2],'sample':v[1]['sample'],'validation':v[1]['validation'],'metrics':v[1]['metrics']} for a,v in evidence.items()},'validation':validation,'mageck_version':version,'mageck_executable_sha256':sha(shutil.which('mageck')),'RRA_executable_sha256':sha(shutil.which('RRA')),'seed_option':extra or 'No seed option advertised; same-input repeats retained','commands':commands,'comparisons':summaries,'replicate_descriptives':descsum,'source_sha256':sha(__file__),'source_commit':os.environ.get('GITHUB_SHA'),'python':sys.version,'platform':platform.platform()}
 dump(out/'summary.json',sums);make_report(out,sums,summaries)
 dump(out/'completion.json',{'completion':'complete','files':{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}})
 print(json.dumps({'primary':summaries[0],'eligible_guides':len(eligible),'genes':sums['eligible_genes'],'reads':sums['reads'],'repeat_checks':summaries[-2:]},sort_keys=True),flush=True)

class Tests(unittest.TestCase):
 def test_strict_counts(self):
  for x in ('-1','2.3',' 3','nan','\u0661'):
   with self.assertRaises(ValueError):natural(x)
  self.assertEqual(natural('0'),0)
 def test_comparisons_and_boundaries(self):
  def row(f,lfc,r):return {'num':'4','neg|lfc':str(lfc),'neg|fdr':str(f),'neg|rank':str(r)}
  a={'a':row(.049,-1,1),'b':row(.06,-.2,2)};b={'a':row(.051,-1,2),'b':row(.02,-.9,1)}
  s,d=compare(a,b,'event','joint','synthetic_unit_test');h=s['thresholds'][1]
  self.assertEqual((h['both'],h['only_left'],h['only_right'],h['neither']),(0,1,1,0));self.assertEqual(s['genes_abs_lfc_difference_ge_0_5'],1);self.assertAlmostEqual(s['spearman_negative_rank'],-1)
 def test_absent_genes_explicit(self):
  r={'num':'1','neg|lfc':'0','neg|fdr':'1','neg|rank':'1'};s,_=compare({'a':r},{'a':r,'b':r},'a','b','test')
  self.assertEqual(s['missing_left'],['b']);self.assertEqual(s['common_genes'],1)
def main():
 if '--test' in sys.argv:unittest.main(argv=[sys.argv[0]],verbosity=2);return
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',required=True);p.add_argument('--out',required=True);build(p.parse_args())
if __name__=='__main__':main()
