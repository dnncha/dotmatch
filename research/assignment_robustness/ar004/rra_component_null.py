#!/usr/bin/env python3
"""Actual RRA component with exact marginal-null p values and empirical sharing.
Not the full MAGeCK pipeline and not biological false-discovery estimation.
"""
from __future__ import annotations
import argparse,csv,gzip,json,os,subprocess,time
from pathlib import Path
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
from scipy import sparse,stats
import read_dependence as core


def guide_classes(folder,data):
    rows=data['guide_rows'];idx={r['id']:i for i,r in enumerate(rows)}
    sets=[];numbers=[]
    for i,n in enumerate(data['singletons']):
        if n:sets.append((i,));numbers.append(int(n))
    for r in core.read(folder/'multiple-guide-read-classes.tsv'):
        sets.append(tuple(idx[i] for i in json.loads(r['target_ids'])));numbers.append(int(r['read_count']))
    rr=[];cc=[]
    for j,xs in enumerate(sets):
        rr.extend([j]*len(xs));cc.extend(xs)
    A=sparse.csr_matrix((np.ones(len(rr),dtype=np.int64),(rr,cc)),shape=(len(sets),len(rows)))
    ns=np.array(numbers,dtype=np.int64);totals=np.asarray(ns@A).reshape(-1)
    expected=np.array([int(r['upstream_count_events']) for r in rows])
    if not np.array_equal(totals,expected):raise AssertionError('All-guide incidence reconstruction failed')
    return A,ns,totals


def execute(args):
    out=Path(args.out);out.mkdir(parents=True,exist_ok=False);binary=Path(args.binary).resolve();summaries=[]
    for a,accession in enumerate(('ERR376998','SRR8297997')):
        folder=Path(args.prior)/('assignment-replay-'+accession)/'audit-output'
        data=core.load_classes(folder);A,ns,n=guide_classes(folder,data)
        rows=data['guide_rows'];active=np.flatnonzero(n>0);genes={rows[i]['gene'] for i in active};target=out/accession;target.mkdir()
        seed=20260910+a
        def task(rep,arm):
            rng=np.random.default_rng(np.random.SeedSequence([seed,rep,0 if arm=='shared_records' else 1]))
            k=np.asarray(rng.binomial(ns,.5)@A).reshape(-1) if arm=='shared_records' else rng.binomial(n,.5)
            nn=n[active];kk=k[active]
            p=stats.binom.cdf(kk-1,nn,.5)+rng.random(len(active))*stats.binom.pmf(kk,nn,.5)
            if np.any(p<0) or np.any(p>1) or not np.isfinite(p).all():raise AssertionError('Invalid exact randomized p')
            prefix=target/f'{arm}-{rep:03d}';inputfile=Path(str(prefix)+'.input.tsv');outputfile=Path(str(prefix)+'.output.tsv');log=Path(str(prefix)+'.log')
            with inputfile.open('w') as f:
                f.write('sgrna\tsymbol\tpool\tp.low\tprob\tchosen\n')
                for j in np.argsort(p,kind='stable'):
                    i=active[j];f.write(f"{rows[i]['id']}\t{rows[i]['gene']}\tlist\t{p[j]:.17g}\t1\t1\n")
            command=[str(binary),'-i',str(inputfile),'-o',str(outputfile),'-p','0.1','--permutation','100']
            with log.open('w') as f:subprocess.run(command,check=True,stdout=f,stderr=subprocess.STDOUT,timeout=300)
            result=core.read(outputfile)
            if {r['group_id'] for r in result}!=genes or len(result)!=len(genes):raise AssertionError('RRA annotation population changed')
            found=sorted(r['group_id'] for r in result if float(r['FDR'])<=.05)
            summary=dict(accession=accession,replicate=rep,arm=arm,seed=seed,active_guides=len(active),gene_groups=len(genes),discoveries=len(found),raw_gene_p_le05=sum(float(r['p'])<=.05 for r in result),guide_p_le01=float((p<=.01).mean()),guide_p_le05=float((p<=.05).mean()),guide_p_le10=float((p<=.1).mean()),selected_genes=json.dumps(found),command=json.dumps(command))
            for path in (inputfile,outputfile,log):
                b=path.read_bytes()
                with path.with_suffix(path.suffix+'.gz').open('wb') as raw,gzip.GzipFile(fileobj=raw,mode='wb',filename='',mtime=0) as f:f.write(b)
                path.unlink()
            return summary
        results=[]
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            jobs=[ex.submit(task,rep,arm) for rep in range(1,101) for arm in ('shared_records','independent_events')]
            for f in as_completed(jobs):
                results.append(f.result())
                if len(results)%20==0:print(accession,len(results),'of 200 actual RRA executions',flush=True)
        results.sort(key=lambda r:(r['replicate'],r['arm']));core.tsv(target/'all-experiments.tsv',results)
        armresults={};gene_rows=[]
        for arm in ('shared_records','independent_events'):
            r=[v for v in results if v['arm']==arm];v=np.array([x['discoveries'] for x in r]);success=int((v>0).sum());ci=stats.binomtest(success,100).proportion_ci(.95,method='exact')
            freq=Counter(g for rr in r for g in json.loads(rr['selected_genes']))
            for g in sorted(genes):gene_rows.append(dict(arm=arm,gene=g,null_replicates_selected=freq[g]))
            armresults[arm]=dict(experiments=100,fraction_with_any_false_discovery=success/100.,binomial_95ci=[ci.low,ci.high],mean_discoveries=float(v.mean()),max_discoveries=int(v.max()),mean_guide_p_le05=float(np.mean([rr['guide_p_le05'] for rr in r])),mean_raw_gene_p_le05=float(np.mean([rr['raw_gene_p_le05']/len(genes) for rr in r])))
        core.tsv(target/'all-gene-selection-frequencies.tsv',gene_rows)
        summary=dict(accession=accession,scope='actual_RRA_component_only_with_exact_randomized_marginal_null_pvalues',source_record_count=data['manifest']['statistics']['reads'],arms=armresults,active_guides=len(active),gene_groups=len(genes),seed=seed,alpha_percentile=.1,permutations=100,binary_sha256=core.sha(binary),amendment_commit='9b090cf74123fa881440157b27a5f3e5248e5edd')
        core.dump(target/'summary.json',summary);summaries.append(summary);print(json.dumps(summary),flush=True)
    core.dump(out/'completion.json',dict(completion='complete',source_sha256=core.sha(__file__),summaries=summaries,files={str(p.relative_to(out)):core.sha(p) for p in sorted(out.rglob('*')) if p.is_file()}))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--prior',required=True);ap.add_argument('--binary',required=True);ap.add_argument('--out',required=True);ap.add_argument('--workers',type=int,default=4);execute(ap.parse_args())
