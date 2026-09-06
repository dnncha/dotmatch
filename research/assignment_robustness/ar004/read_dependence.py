#!/usr/bin/env python3
"""AR004 conditional technical-null experiment on complete read-event classes.

Not MAGeCK, not biological replication and not a test of true molecular origin.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, platform, sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
import numpy as np
from scipy import sparse, stats, signal

ACCESSIONS = ('ERR376998', 'ERR376999', 'SRR8297997')
PROTOCOL = 'c41be48cbea1a330ed93c3dae79922f5824e5f67'


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def dump(p, x):
    Path(p).write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')


def tsv(p, rows, fields=None):
    rows=list(rows)
    fields=fields or list(rows[0])
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)


def read(p):
    with Path(p).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def natural(s):
    if not isinstance(s,str) or not s.isascii() or not s.isdecimal():raise ValueError('Invalid count')
    return int(s)


def bh_count(p):
    a=np.sort(p,axis=1)
    allowed=np.arange(1,a.shape[1]+1)[None,:]*(0.05/a.shape[1])
    return np.max(np.where(a<=allowed,np.arange(1,a.shape[1]+1)[None,:],0),axis=1)


def load_classes(folder):
    folder=Path(folder)
    m=json.loads((folder/'completion.json').read_text())
    if m['completion']!='complete':raise ValueError('Incomplete source')
    for name,value in m['files'].items():
        p=folder/name
        if not p.resolve().is_relative_to(folder.resolve()) or sha(p)!=value:raise ValueError('Source hash mismatch '+name)
    if m['statistics'].get('repeated_same_target_events',-1)!=0:
        raise ValueError('Need full per-target multiplicities; distinct-ID classes are insufficient')
    rows=read(folder/'all-guide-counts.tsv')
    if len({r['id'] for r in rows})!=len(rows):raise ValueError('Duplicate IDs')
    names=sorted({r['gene'] for r in rows});glookup={g:i for i,g in enumerate(names)}
    idmap={r['id']:i for i,r in enumerate(rows)}
    gidx=np.array([glookup[r['gene']] for r in rows],dtype=int)
    counts=np.array([natural(r['upstream_count_events']) for r in rows],dtype=np.int64)
    left=counts.copy();raw_multi=[];pairs=Counter();seen=set()
    for r in read(folder/'multiple-guide-read-classes.tsv'):
        ids=tuple(idmap[i] for i in json.loads(r['target_ids']))
        if len(ids)<2 or len(set(ids))!=len(ids) or tuple(sorted(ids)) in seen:raise ValueError('Invalid/repeated class')
        seen.add(tuple(sorted(ids)))
        n=natural(r['read_count'])
        if n<1:raise ValueError('Empty class')
        if set(json.loads(r['genes']))!={rows[i]['gene'] for i in ids}:raise ValueError('Gene labels mismatch')
        left[list(ids)]-=n
        raw_multi.append((ids,n))
        for a,b in combinations(sorted(ids),2):pairs[(a,b)]+=n
    if np.any(left<0):raise ValueError('Negative singleton residue')
    matched=int(left.sum())+sum(n for _,n in raw_multi)
    events=int(counts.sum());extra=sum(n*(len(ids)-1) for ids,n in raw_multi)
    if (matched,events,extra)!=(m['statistics']['matched_reads'],m['statistics']['count_events'],m['statistics']['extra_events']):raise ValueError('Read budget mismatch')
    if sum(n for _,n in raw_multi)!=m['statistics']['multiple_distinct_guides']:raise ValueError('Multiple class budget')
    # Coalesce identical gene-weight patterns, preserving every induced dependence.
    classes=Counter()
    for i,n in enumerate(left):
        if n:classes[((int(gidx[i]),1),)]+=int(n)
    for ids,n in raw_multi:
        key=tuple(sorted(Counter(int(gidx[i]) for i in ids).items()))
        classes[key]+=n
    keys=sorted(classes); ns=np.array([classes[k] for k in keys],dtype=np.int64)
    rr=[];cc=[];vv=[];hist=[Counter() for _ in names]
    for j,k in enumerate(keys):
        for g,weight in k:
            rr.append(j);cc.append(g);vv.append(weight);hist[g][weight]+=int(ns[j])
    B=sparse.csr_matrix((np.array(vv,dtype=np.int64),(rr,cc)),shape=(len(keys),len(names)))
    C=np.asarray(ns@B).reshape(-1); Q=np.asarray(ns@B.multiply(B)).reshape(-1)
    U=np.asarray(ns@(B!=0)).reshape(-1)
    C_check=np.bincount(gidx,weights=counts,minlength=len(names)).astype(np.int64)
    Q_check=C_check.copy();active=np.bincount(gidx,weights=counts>0,minlength=len(names)).astype(int)
    r2=np.zeros(len(names))
    pairrows=[]
    for (i,j),n in sorted(pairs.items()):
        if not 0<n<=min(counts[i],counts[j]):raise ValueError('Impossible overlap')
        corr=n/math.sqrt(int(counts[i])*int(counts[j]))
        same=gidx[i]==gidx[j]
        if same:
            Q_check[gidx[i]]+=2*n;r2[gidx[i]]+=2*corr*corr
        pairrows.append(dict(guide_a=rows[i]['id'],guide_b=rows[j]['id'],gene_a=rows[i]['gene'],gene_b=rows[j]['gene'],same_gene=int(same),shared_records=n,count_a=int(counts[i]),count_b=int(counts[j]),null_correlation=corr,jaccard=n/(int(counts[i])+int(counts[j])-n)))
    if not np.array_equal(C,C_check) or not np.array_equal(Q,Q_check):raise AssertionError('Independent pair/class covariance check failed')
    if int(ns.sum())!=matched:raise AssertionError('Coalescing lost records')
    effective=np.divide(active*active,active+r2,out=np.zeros_like(r2),where=active>0)
    return dict(manifest=m,names=names,C=C,Q=Q,U=U,active=active,keff=effective,B=B,ns=ns,hist=hist,pairs=pairrows,guide_rows=rows,singletons=left)


def exact_law(hist):
    """Convolve the full weighted-binomial distribution; C - 2T is symmetric."""
    pmf=np.array([1.])
    for w,n in sorted(hist.items()):
        padded=np.zeros(w*n+1)
        padded[::w]=stats.binom.pmf(np.arange(n+1),n,0.5)
        if len(pmf)*len(padded)<500000:pmf=np.convolve(pmf,padded)
        else:pmf=signal.fftconvolve(pmf,padded)
    if np.min(pmf)<-1e-12:raise AssertionError('Negative convolution probability')
    pmf=np.maximum(pmf,0);pmf/=pmf.sum()
    C=sum(w*n for w,n in hist.items());Q=sum(w*w*n for w,n in hist.items())
    d=2*np.arange(len(pmf))-C
    if abs(pmf@d)>1e-7 or not np.isclose(pmf@(d*d),Q,rtol=1e-9):raise AssertionError('Exact distribution moments failed')
    absd=np.abs(d);maxd=int(absd.max())
    mass=np.bincount(absd,weights=pmf,minlength=maxd+1)
    tail=np.cumsum(mass[::-1])[::-1]
    return dict(C=C,Q=Q,pmf=pmf,d=d,p_exact=tail[absd])


def experiment(data, out, replicates=2000, seed=20260906):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    names=data['names'];C=data['C'];Q=data['Q'];B=data['B'];ns=data['ns']
    elig=C>=100;infl=np.divide(Q,C,out=np.ones(len(C)),where=C>0)
    affected=elig&(infl>=1.25); severe=elig&(infl>=1.5)
    groups={'all_eligible':elig,'variance_ratio_ge1.25':affected,'variance_ratio_ge1.5':severe,'no_within_gene_reuse':elig&(Q==C)}
    count_names=('naive','record_covariance','independent_events')
    hits={k:np.zeros(len(C),dtype=np.int64) for k in count_names}
    rng=np.random.default_rng(seed);ind_rng=np.random.default_rng(seed+100000)
    all_splits=[]
    # Structural selection fixed before any random outcomes are generated.
    maxratio=max(infl[elig]);witness=int(next(i for i,g in enumerate(names) if elig[i] and infl[i]==maxratio))
    control=int(min(np.flatnonzero(elig&(Q==C)),key=lambda i:(abs(int(C[i])-int(C[witness])),names[i])))
    witness_D={witness:[],control:[]}
    for start in range(0,replicates,25):
        b=min(25,replicates-start)
        coin=2*rng.binomial(ns,0.5,size=(b,len(ns)))-ns
        D=np.asarray(coin@B)
        independent=2*ind_rng.binomial(C,0.5,size=(b,len(C)))-C
        p={
          'naive':2*stats.norm.sf(np.abs(D)/np.sqrt(np.maximum(C,1))),
          'record_covariance':2*stats.norm.sf(np.abs(D)/np.sqrt(np.maximum(Q,1))),
          'independent_events':2*stats.norm.sf(np.abs(independent)/np.sqrt(np.maximum(C,1)))
        }
        for k in count_names:hits[k]+=(p[k]<0.05).sum(axis=0)
        for i in witness_D:witness_D[i].extend(int(x) for x in D[:,i])
        bh={k:bh_count(p[k][:,elig]) for k in count_names}
        for r in range(b):
            for k in count_names:
                v=dict(split=start+r+1,method=k,bh_discoveries=int(bh[k][r]))
                for g,mask in groups.items():v[g+'_fraction']=float((p[k][r,mask]<.05).mean()) if mask.any() else None
                all_splits.append(v)
    per_gene=[]
    for i,g in enumerate(names):
        per_gene.append(dict(gene=g,eligible=int(elig[i]),count_events=int(C[i]),distinct_supporting_records=int(data['U'][i]),null_variance=int(Q[i]),variance_ratio=float(infl[i]),effective_record_weight=float(C[i]*C[i]/Q[i]) if Q[i] else 0.,active_guide_rows=int(data['active'][i]),effective_technical_guide_dimensions=float(data['keff'][i]),**{k+'_false_positive_rate':float(hits[k][i]/replicates) for k in count_names},record_multiplicity_histogram=json.dumps(dict(sorted(data['hist'][i].items())),sort_keys=True)))
    tsv(out/'all-gene-dependence.tsv',per_gene)
    tsv(out/'all-guide-pair-overlaps.tsv',data['pairs'])
    tsv(out/'all-null-splits.tsv',all_splits)
    selected=[]
    for i in (witness,control):
        law=exact_law(data['hist'][i]);pmf=law['pmf'];d=law['d'];p_exact=law['p_exact']
        Cg=law['C'];Qg=law['Q']
        pn=2*stats.norm.sf(abs(d)/math.sqrt(Cg));pc=2*stats.norm.sf(abs(d)/math.sqrt(Qg))
        result=dict(gene=names[i],role='largest_structural_variance_ratio' if i==witness else 'unaffected_coverage_matched_control',count_events=Cg,null_variance=Qg,variance_ratio=Qg/Cg,exact_law_nominal_naive_false_positive_rate=float(pmf[pn<.05].sum()),exact_law_nominal_covariance_false_positive_rate=float(pmf[pc<.05].sum()),exact_law_exact_test_false_positive_rate=float(pmf[p_exact<=.05].sum()),observed_naive_rate=float(hits['naive'][i]/replicates),observed_covariance_rate=float(hits['record_covariance'][i]/replicates),histogram=dict(data['hist'][i]))
        selected.append(result)
        tsv(out/(result['role']+'-distribution.tsv'),(dict(difference=int(v),probability=float(pr),naive_p=float(p1),covariance_p=float(p2),exact_p=float(p3)) for v,pr,p1,p2,p3 in zip(d,pmf,pn,pc,p_exact)))
    metrics={}
    for k in count_names:
        rr=[r for r in all_splits if r['method']==k];metrics[k]={}
        for g,mask in groups.items():
            vals=[r[g+'_fraction'] for r in rr if r[g+'_fraction'] is not None]
            metrics[k][g]=dict(annotations=int(mask.sum()),mean_false_positive_fraction=float(np.mean(vals)) if vals else None,split_std=float(np.std(vals,ddof=1)) if vals else None,monte_carlo_se=float(np.std(vals,ddof=1)/math.sqrt(replicates)) if vals else None)
        b=np.array([r['bh_discoveries'] for r in rr]);metrics[k]['bh']=dict(mean_discoveries=float(b.mean()),fraction_splits_any_discovery=float((b>0).mean()),maximum=int(b.max()))
    summary=dict(accession=data['manifest']['accession'],protocol_commit=PROTOCOL,scope='complete_archived_upstream_read_event_classes_and_new_random_label_experiment',source_reads=data['manifest']['statistics']['reads'],matched_records=int(ns.sum()),count_events=int(C.sum()),reference_guides=len(data['guide_rows']),original_gene_annotations=len(names),eligible_gene_annotations=int(elig.sum()),replicates=replicates,seed=seed,nominal_alpha=.05,variance_ratio_thresholds={str(t):int((elig&(infl>=t)).sum()) for t in (1.25,1.5,2)},guide_pairs_with_shared_records=len(data['pairs']),maximum_variance_ratio=float(maxratio),metrics=metrics,selected_examples=selected)
    dump(out/'summary.json',summary)
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--prior',required=True);p.add_argument('--out',required=True);p.add_argument('--replicates',type=int,default=2000);a=p.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=False);summaries=[]
    for idx,accession in enumerate(ACCESSIONS):
        candidates=[x.parent for x in Path(a.prior).rglob('completion.json') if x.parent.name=='audit-output' and accession in str(x)]
        if len(candidates)!=1:raise ValueError('Exactly one archive evidence required')
        data=load_classes(candidates[0]);print(accession,'reconciled',len(data['names']),'gene annotations',flush=True)
        result=experiment(data,out/accession,a.replicates,20260906+idx)
        summaries.append(result);print(json.dumps(result,sort_keys=True),flush=True)
    dump(out/'completion.json',dict(completion='complete',protocol_commit=PROTOCOL,summaries=summaries,environment=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__),source_sha256=sha(__file__),files={str(x.relative_to(out)):sha(x) for x in sorted(out.rglob('*')) if x.is_file()}))

if __name__=='__main__':main()
