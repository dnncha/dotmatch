#!/usr/bin/env python3
"""Descriptive paired-policy analysis of actual MAGeCK outputs, not truth calls."""
import argparse,itertools,json,math
from pathlib import Path
import numpy as np
from scipy import stats
import read_dependence as c

ARMS=('exact','radius_k1','best_k1','upstream_count_events')
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
    inp=Path(a.input);out=Path(a.out);out.mkdir(parents=True,exist_ok=False)
    m=json.loads((inp/'completion.json').read_text());assert m['completion']=='complete'
    for name,h in m['files'].items():assert c.sha(inp/name)==h,name
    tables={p:c.read(inp/(p+'.mageck.gene_summary.txt')) for p in ARMS};maps={}
    for p,t in tables.items():
        maps[p]={r['id']:r for r in t}
        assert len(maps[p])==len(t)
        for r in t:
            for k,v in r.items():
                if k!='id':assert math.isfinite(float(v)),(p,k,v)
    names=sorted(maps['exact']);assert all(set(v)==set(names) for v in maps.values())
    sizes=[];comparisons=[];allrows=[]
    for p in ARMS:
        sizes.append(dict(policy=p,gene_groups=len(names),negative_fdr_le05=sum(float(r['neg|fdr'])<=.05 for r in tables[p]),positive_fdr_le05=sum(float(r['pos|fdr'])<=.05 for r in tables[p])))
    for a,b in itertools.combinations(ARMS,2):
        for d in ('neg','pos'):
            ha={g for g in names if float(maps[a][g][d+'|fdr'])<=.05};hb={g for g in names if float(maps[b][g][d+'|fdr'])<=.05}
            ranksa=[float(maps[a][g][d+'|rank']) for g in names];ranksb=[float(maps[b][g][d+'|rank']) for g in names]
            lfc_a=np.array([float(maps[a][g][d+'|lfc']) for g in names]);lfc_b=np.array([float(maps[b][g][d+'|lfc']) for g in names])
            delta=lfc_b-lfc_a
            rho=float(stats.spearmanr(ranksa,ranksb).statistic)
            comparisons.append(dict(policy_a=a,policy_b=b,direction=d,gene_groups=len(names),a_hits=len(ha),b_hits=len(hb),shared_hits=len(ha&hb),a_only=len(ha-hb),b_only=len(hb-ha),symmetric_difference=len(ha^hb),jaccard=len(ha&hb)/len(ha|hb) if ha|hb else 1.,rank_spearman=rho,median_absolute_lfc_delta=float(np.median(abs(delta))),lfc_delta_ge_half=int((abs(delta)>=.5).sum()),strong_direction_reversal=int(((lfc_a*lfc_b<0)&(abs(lfc_a)>=.5)&(abs(lfc_b)>=.5)).sum())))
            for i,g in enumerate(names):
                x=maps[a][g];y=maps[b][g]
                allrows.append(dict(gene=g,policy_a=a,policy_b=b,direction=d,guide_rows_a=int(x['num']),guide_rows_b=int(y['num']),rank_a=int(x[d+'|rank']),rank_b=int(y[d+'|rank']),fdr_a=float(x[d+'|fdr']),fdr_b=float(y[d+'|fdr']),p_a=float(x[d+'|p-value']),p_b=float(y[d+'|p-value']),lfc_a=float(x[d+'|lfc']),lfc_b=float(y[d+'|lfc']),lfc_delta=float(delta[i]),good_guides_a=int(x[d+'|goodsgrna']),good_guides_b=int(y[d+'|goodsgrna']),hit_a=int(g in ha),hit_b=int(g in hb),changed_hit=int((g in ha)!=(g in hb))))
    c.tsv(out/'hit-totals.tsv',sizes);c.tsv(out/'all-policy-comparisons.tsv',comparisons);c.tsv(out/'all-gene-policy-comparisons.tsv',allrows)
    changed=[r for r in allrows if r['changed_hit']];c.tsv(out/'all-changed-hit-classifications.tsv',changed)
    summary=dict(completion='complete',scope='actual_full_MAGeCK_0.5.9.5_on_original_plasmid_ESC1_ESC2',independent_external_studies=0,biological_truth_estimated=False,p_value_threshold=.05,threshold_type='MAGeCK_reported_FDR',contrasts=dict(control=['ERR376998'],cellular_replicates=['ERR376999','ERR377000']),totals=sizes,comparisons=comparisons,source_sha256=c.sha(__file__),raw_extension=json.loads((inp/'raw-audit.json').read_text()))
    c.dump(out/'summary.json',summary)
    c.dump(out/'completion.json',dict(completion='complete',files={p.name:c.sha(p) for p in sorted(out.iterdir()) if p.is_file()}))
    print(json.dumps(dict(totals=sizes,comparisons=comparisons),indent=2))
if __name__=='__main__':main()
