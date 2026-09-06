"""Reproducible, exploratory EditWitness assay and captured-read audit.

No function treats amplified reads as independent original molecules. No function
asserts that an exact primer match predicts laboratory amplification efficiency.
The read classifier identifies a restricted alignment pattern, not a validated
biological editing outcome. All fractions name their observable denominator.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, math, re, statistics
from collections import Counter, defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET

HUMAN = 'human_ATXN2_GRCh38_local'
CUTS = (4916, 5383)
ANCHORS = ((4700,4800),(5500,5600))

def reverse_complement(sequence: str) -> str:
    if not sequence or set(sequence.upper()) - set('ACGTN'):
        raise ValueError('nonempty ACGTN sequence required')
    return sequence.upper().translate(str.maketrans('ACGTN','TGCAN'))[::-1]

def sites(sequence: str, oligo: str) -> list[int]:
    if not oligo or set(oligo)-set('ACGT'):
        raise ValueError('oligo must contain only ACGT')
    return [i for i in range(len(sequence)-len(oligo)+1) if sequence.startswith(oligo,i)]

def products(sequence: str, forward: str, reverse: str) -> list[dict]:
    """Enumerate every exact inward-facing product in BOTH orientations.

    The two named oligos are not assumed to bind only their original positions.
    No length/yield/specificity claim is made. Coordinates are zero-based half-open.
    """
    answer=[]
    for strand,left,right in [('+',forward,reverse_complement(reverse)),
                              ('-',reverse,reverse_complement(forward))]:
        for start in sites(sequence,left):
            for r_start in sites(sequence,right):
                if start+len(left) <= r_start:
                    answer.append({'start0':start,'end0':r_start+len(right),
                                   'length':r_start+len(right)-start,'orientation':strand})
    return answer

def cigar_blocks(hit: dict) -> tuple[list[tuple[int,int]], list[tuple[int,int]]]:
    """Return reference-aligned (M/= /X) blocks and true D gaps (not N skips)."""
    cigar=hit['cigar']; tokens=re.findall(r'(\d+)([MIDNSHP=X])',cigar)
    if not tokens or ''.join(n+op for n,op in tokens)!=cigar:
        raise ValueError('malformed CIGAR')
    pos=hit['r_start']; blocks=[]; deletions=[]
    for n,op in tokens:
        size=int(n)
        if size<=0: raise ValueError('zero-length CIGAR operation')
        if op in 'M=X': blocks.append((pos,pos+size)); pos+=size
        elif op=='D': deletions.append((pos,pos+size)); pos+=size
        elif op=='N': pos+=size
    if pos!=hit['r_end']:
        raise ValueError('CIGAR/reference endpoint inconsistency')
    return blocks,deletions

def anchor_coverage(blocks: list[tuple[int,int]], anchor: tuple[int,int]) -> int:
    # CIGAR blocks do not overlap; M is alignment coverage, NOT guaranteed identity.
    return sum(max(0,min(end,anchor[1])-max(start,anchor[0])) for start,end in blocks)

def prepare_read(record: dict) -> dict:
    hits=[h for h in record['alignments'] if h['contig']==HUMAN and h['is_primary']]
    if not hits:
        return {'status':'no_primary_human_hit','deletion_compatible':False}
    hit=max(hits,key=lambda h:(h['match_bases'],h['mapq'],h['block_length']))
    blocks,gaps=cigar_blocks(hit)
    return {'status':'prepared','mapq':hit['mapq'],'strand':hit['strand'],
            'anchor_coverage':[anchor_coverage(blocks,a) for a in ANCHORS], 'gaps':gaps}

def classify_prepared(prepared: dict, mapq: int=20, endpoint_tolerance: int=30) -> dict:
    if mapq < 0 or endpoint_tolerance < 0: raise ValueError('thresholds must be nonnegative')
    if prepared['status']=='no_primary_human_hit':return prepared
    if prepared['mapq']<mapq:
        return {'status':'below_mapq','deletion_compatible':False}
    cov=prepared['anchor_coverage']
    if min(cov)<80:
        return {'status':'insufficient_flanking_alignment','deletion_compatible':False,
                'anchor_coverage':cov}
    compatible=[(start,end) for start,end in prepared['gaps']
                if abs(start-CUTS[0])<=endpoint_tolerance
                and abs(end-CUTS[1])<=endpoint_tolerance]
    return {'status':'eligible_two_flank_read','deletion_compatible':bool(compatible),
            'anchor_coverage':cov,'compatible_gaps':compatible,
            'mapq':prepared['mapq'],'strand':prepared['strand']}

def classify_read(record: dict, mapq: int=20, endpoint_tolerance: int=30) -> dict:
    return classify_prepared(prepare_read(record), mapq, endpoint_tolerance)

def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows: raise ValueError('refusing ambiguous empty CSV')
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader();writer.writerows(rows)

def save(path: Path, value: object) -> None:
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def verify_receipts(root: Path) -> list[dict]:
    checks=[]
    for folder in [root/'inputs',root/'inputs/assay_sources',root/'inputs/calibration']:
        for rec in json.loads((folder/'receipts.json').read_text()):
            if not rec.get('ok'):
                checks.append({'path':str(folder.relative_to(root)/rec['name']),
                               'download_ok':False,'hash_ok':None,'content_valid':False,
                               'reason':rec.get('error','download failed')});continue
            data=(folder/rec['name']).read_bytes()
            valid_hash=hashlib.sha256(data).hexdigest()==rec['sha256'] and len(data)==rec['bytes']
            if not valid_hash: raise ValueError('input checksum mismatch: '+rec['name'])
            content_valid=not rec['name'].endswith('.pdf') or data.startswith(b'%PDF-')
            checks.append({'path':str(folder.relative_to(root)/rec['name']),
                           'download_ok':True,'hash_ok':True,'content_valid':content_valid,
                           'reason':'validated source bytes' if content_valid else 'HTTP 200 HTML masquerading as PDF; rejected'})
    save(root/'results/input_validation.json',checks)
    return checks

def audit_metadata(root: Path) -> dict:
    table=list(csv.DictReader((root/'inputs/PRJNA916868.tsv').open(),delimiter='\t'))
    if len({r['run_accession'] for r in table})!=len(table): raise ValueError('duplicate run id')
    result=[]
    for row in table:
        exp=ET.parse(root/'inputs'/(row['experiment_accession']+'.xml')).find('.//EXPERIMENT')
        sample=ET.parse(root/'inputs/assay_sources'/(row['sample_accession']+'.xml')).find('.//SAMPLE')
        attrs={e.findtext('TAG'):e.findtext('VALUE') for e in sample.findall('.//SAMPLE_ATTRIBUTE')}
        alias=exp.attrib['alias']
        mode='PCR' if '_PCR_' in alias else 'nCATS'
        processed='PCR-basecalling-unspecified' if mode=='PCR' else alias.split('_')[1]
        result.append({'run':row['run_accession'],'biosample':row['sample_accession'],
                       'specimen_alias':sample.attrib['alias'],'run_alias':alias,'assay':mode,
                       'processing':processed,'treatment':attrs['treatment'],'strain':attrs['strain'],
                       'read_count':int(row['read_count']),'base_count':int(row['base_count']),
                       'fastq_bytes':int(row['fastq_bytes']),
                       'source_experiment':row['experiment_accession']})
    grouped=defaultdict(list)
    for r in result: grouped[r['biosample']].append(r)
    pairs=[]
    for biosample,group in sorted(grouped.items()):
        pcr=[r for r in group if r['assay']=='PCR']; sup=[r for r in group if r['processing']=='SUP']
        if pcr and sup:
            if len(pcr)!=1 or len(sup)!=1: raise ValueError('ambiguous paired run choice')
            pairs.append({'biosample':biosample,'specimen':pcr[0]['specimen_alias'],
                          'treatment':pcr[0]['treatment'],'PCR_run':pcr[0]['run'],
                          'nCATS_SUP_run':sup[0]['run'],
                          'pairing_evidence':'same deposited BioSample; exact aliquot chain unverified'})
    profile={'runs':len(result),'biosamples':len(grouped),'same_biosample_pairs':len(pairs),
             'processing_counts':dict(Counter(r['processing'] for r in result)),
             'compressed_bytes_all_runs':sum(r['fastq_bytes'] for r in result),
             'independent_units_warning':'HAC/SUP/5mC are processing outputs, not biological replicates',
             'original_molecule_inference_eligible':False,
             'failed_gates':['per-state capture calibration absent','specimen haplotypes unverified',
                             'same-extraction/aliquot chain not independently established',
                             'local ddPCR is not a full-deletion denominator',
                             'read-to-target-copy mapping is unresolved for a multicopy transgene']}
    write_csv(root/'results/run_manifest.csv',result)
    write_csv(root/'results/specimen_pairs.csv',pairs)
    save(root/'results/metadata_profile.json',profile)
    return profile

def audit_geometry(root: Path) -> dict:
    manifest=json.loads((root/'protocol/assay_geometry.json').read_text())
    native=json.loads((root/'inputs/assay_sources/hg38-ATXN2-region.json').read_text())['dna'].upper()
    if hashlib.sha256(native.encode()).hexdigest()!=manifest['native_reference']['dna_sha256']:
        raise ValueError('reference DNA digest differs')
    local=reverse_complement(native[146000:155000])
    oligos=manifest['oligos']
    placements={k:{'+':sites(local,v),'-':sites(local,reverse_complement(v))} for k,v in oligos.items()}
    if any(sum(len(a) for a in pair.values())!=1 for pair in placements.values()):
        raise ValueError('a source oligo lacks a unique exact local-reference match')
    # Explicitly constructed 72-codon repeat, NOT a recovered specimen haplotype.
    nominal=local[:5194]+'CAG'*72+local[5263:]
    deletion=nominal[:4916]+nominal[5530:]
    states={'native_23_codon_reference':local,'nominal_72_codon_model':nominal,
            'nominal_precise_dual_cut_deletion':deletion}
    product_table=[]
    for state,sequence in states.items():
        for assay, fwd,rev in [('outer_PCR','pcr_forward','pcr_reverse'),
                               ('g4_local_ddPCR','dd4_forward','dd4_reverse'),
                               ('g5_local_ddPCR','dd5_forward','dd5_reverse')]:
            found=products(sequence,oligos[fwd],oligos[rev])
            product_table.append({'state':state,'assay':assay,'exact_products':len(found),
                                  'product_lengths_bp':';'.join(str(p['length']) for p in found)})
    result={'placements_native0':placements,'native_intercut_bp':467,
            'nominal72_intercut_bp':614,'products':product_table,
            'nominal72_deletion_removes':['dd4_reverse','dd4_probe_ref','dd5_forward','dd5_probe_ref'],
            'conclusion':'The precise dual-cut deletion is outside both short local drop-off amplicons under exact primer matching; local indel fractions are not full-deletion ground truth.',
            'scope':'sequence geometry on native flanks plus a nominal repeat model, not experimental efficiency, genotype truth, or a new discovery of PCR bias'}
    for missing in result['nominal72_deletion_removes']:
        if sites(deletion,oligos[missing]) or sites(deletion,reverse_complement(oligos[missing])):
            raise ValueError('claimed missing site remains in deletion model')
    save(root/'results/assay_geometry_result.json',result)
    write_csv(root/'results/assay_products.csv',product_table)
    (root/'results/native_local_reference.fa').write_text('>native_GRCh38_ATXN2_local_negative_orientation\n'+local+'\n')
    return result

def audit_published_values(root: Path) -> list[dict]:
    # Values exported from the publisher XLSX through artifact_tool, with cell coordinates.
    values=json.loads((root/'results/figure3_cells.json').read_text())
    if isinstance(values,dict): values=values['values']
    rows=[]; guide=treatment=None
    for excelrow in range(47,59):
        r=values[excelrow-1]
        if r[0]: guide=r[0]
        if r[1]: treatment=r[1]
        exclusive=[float(x or 0) for x in r[3:10]]
        rows.append({'excel_row':excelrow,'input_guides':str(guide),'treatment':str(treatment),
                     'replicate':int(r[2]),'exclusive_mutation_sum_percent':sum(exclusive),
                     'large_deletion_percent':float(r[6]),'integration_annotation_percent':float(r[10]),
                     'provenance':'CRISPRLungo source data Fig. 3, Main_Fig3!D'+str(excelrow)+':K'+str(excelrow),
                     'warning':'integration annotation overlaps mutation classes; do not add K to D:J'})
    write_csv(root/'results/published_mutation_values.csv',rows)
    return rows

def audit_calibration(root: Path) -> dict:
    path=root/'inputs/calibration/file-preprocessed-data__19_spikes_mean-count.csv'
    values=list(csv.DictReader(path.open()))
    groups=defaultdict(list)
    for r in values: groups[(r['cycle'],r['mol'])].append(r)
    stats=[]
    for (cycle,concentration),rs in sorted(groups.items()):
        if len({r['spike'] for r in rs})!=len(rs):raise ValueError('duplicate spike in condition')
        shares=[float(r['mean']) for r in rs]
        # Descriptive dispersion only. A condition mean is NOT an independent raw replicate.
        stats.append({'cycles':cycle,'concentration':concentration,'classes':len(rs),
                      'sum_mean_percent':sum(shares),'minimum_mean_percent':min(shares),
                      'maximum_mean_percent':max(shares),'max_min_share_ratio':max(shares)/min(shares),
                      'source_unit':'published class means; raw replicate counts not included in this table'})
    write_csv(root/'results/calibration_condition_profile.csv',stats)
    result={'rows':len(values),'conditions':len(groups),'classes':sorted({r['spike'] for r in values}),
            'full_original_molecule_calibration':False,
            'reasons':['source table reports means and standard deviations, not all replicate counts',
                       'nine constructs in mean table; eight rows per state in separate correction example',
                       'plasmid response is not automatically transferable to mouse-brain nCATS or long amplicons',
                       'no independent held-out truth cohort analysed'],
            'scope':'primary-source calibration audit, not evaluation of CRISPR-A accuracy'}
    save(root/'results/calibration_profile.json',result)
    return result

def audit_reads(root: Path) -> list[dict]:
    source=root/'inputs/reads'
    if not source.exists(): return []
    summaries=[]
    for runpath in sorted(source.iterdir()):
        if not runpath.is_dir():continue
        rec=json.loads((runpath/'receipt.json').read_text())
        if rec['status']!='complete':raise ValueError('read input is not complete: '+runpath.name)
        compact_path=runpath/'compact_manifest.json'
        omissions={}
        if compact_path.exists():
            compact=json.loads(compact_path.read_text())
            if not compact.get('all_source_file_hashes_verified'):raise ValueError('unverified compact source')
            omissions={x['name']:x for x in compact['omitted_files']}
        for line in (runpath/'SHA256SUMS').read_text().splitlines():
            digest,name=line.split('  ',1)
            if Path(name).name!=name:raise ValueError('unsafe artifact path')
            if not (runpath/name).exists():
                omit=omissions.get(name,{})
                if name!='target_candidates.fastq.gz' or omit.get('sha256')!=digest or not omit.get('hash_verified_on_runner'):
                    raise ValueError('unexplained missing source artifact: '+name)
                continue
            if hashlib.sha256((runpath/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('read artifact checksum mismatch')
        counters={(q,t):Counter() for q in [0,20,40] for t in [10,30,100]}
        n=0;seen=set()
        with gzip.open(runpath/'target_candidates.jsonl.gz','rt') as f, gzip.open(root/'results'/(runpath.name+'-audit.jsonl.gz'),'wt') as out:
            for line in f:
                obj=json.loads(line);n+=1
                if obj['read_id'] in seen:raise ValueError('duplicate read id within run')
                seen.add(obj['read_id'])
                prepared=prepare_read(obj)
                for (q,t),counts in counters.items():
                    classified=classify_prepared(prepared,q,t);counts[classified['status']]+=1
                    if classified['deletion_compatible']:counts['deletion_compatible']+=1
                out.write(json.dumps({'read_id':obj['read_id'],**classify_prepared(prepared)},separators=(',',':'))+'\n')
        if n!=rec['target_candidates']:raise ValueError('candidate receipt count mismatch')
        for (q,t),c in counters.items():
            denominator=c['eligible_two_flank_read'];numerator=c['deletion_compatible']
            summaries.append({'run':runpath.name,'mapq_threshold':q,'endpoint_tolerance_bp':t,
                               'total_reads_processed':rec['reads_processed'],'target_candidates':n,
                               'eligible_two_flank_reads':denominator,'deletion_compatible_reads':numerator,
                               'conditional_read_fraction':numerator/denominator if denominator else None,
                               'excluded_no_primary_human':c['no_primary_human_hit'],
                               'excluded_below_mapq':c['below_mapq'],
                               'excluded_flanking_alignment':c['insufficient_flanking_alignment'],
                               'interpretation':'captured-read alignment pattern; not original-molecule editing fraction'})
            if sum(c[s] for s in ['eligible_two_flank_read','no_primary_human_hit','below_mapq','insufficient_flanking_alignment'])!=n:
                raise ValueError('read disposition does not reconcile')
    if summaries: write_csv(root/'results/read_sensitivity.csv',summaries)
    return summaries

def main(root: Path) -> None:
    (root/'results').mkdir(exist_ok=True)
    checks=verify_receipts(root);metadata=audit_metadata(root);geometry=audit_geometry(root)
    published=audit_published_values(root);calibration=audit_calibration(root);reads=audit_reads(root)
    save(root/'results/analysis_receipt.json',{'schema':'editwitness.empirical-audit.v1',
         'metadata':metadata,'source_checks':len(checks),'invalid_content_count':sum(not c['content_valid'] for c in checks),
         'geometry_products_checked':len(geometry['products']),'published_table_rows':len(published),
         'calibration_conditions':calibration['conditions'],'read_sensitivity_rows':len(reads),
         'original_molecule_validation_complete':False})

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    main(p.parse_args().root)
