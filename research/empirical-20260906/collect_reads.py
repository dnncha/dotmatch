"""Exploratory captured-read audit; not an editing caller or biological truth model.

Fetch one prespecified ENA run, verify all compressed bytes, then map every read
to a restricted human ATXN2 reference plus a mouse Atxn2 decoy. Preserve target-
compatible reads and all returned alignments for independent denominator audits.
Mapping qualities are relative to this restricted reference, NOT genome-wide.
No read-level confidence interval or original-molecule fraction is produced.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, io, json, pathlib, platform, resource, time
import urllib.parse, urllib.request
from collections import Counter
from datetime import datetime, timezone
import mappy as mp

RUNS = {'SRR22937744','SRR22937743','SRR22937742','SRR22937757','SRR22937754','SRR22937753'}
NATIVE_SHA = 'fb7de2d6c66dd142917ad6b1c5f9f7e5dd8a886f769a454de7f21e822b65ed74'
MAX_COMPRESSED_BYTES = 1800000000

def get(url: str, limit: int = 5000000) -> bytes:
    request = urllib.request.Request(url, headers={'User-Agent':'EditWitness-research-audit/0.0.2'})
    with urllib.request.urlopen(request, timeout=90) as response:
        data = response.read(limit + 1)
    if len(data) > limit:
        raise ValueError('retrieval size budget exceeded')
    return data

def collect(run: str, output: pathlib.Path) -> None:
    if run not in RUNS:
        raise ValueError('run is not in the frozen exploratory cohort')
    output.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    receipt = {'run':run, 'started_utc':datetime.now(timezone.utc).isoformat(),
               'status':'started','mappy_version':mp.__version__,'python':platform.python_version(),
               'scope':'restricted-reference captured-read audit; no original-molecule or clinical inference'}
    try:
        fields = 'run_accession,sample_accession,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes'
        url = 'https://www.ebi.ac.uk/ena/portal/api/filereport?' + urllib.parse.urlencode(
            {'accession':run,'result':'read_run','format':'tsv','fields':fields})
        metadata = get(url)
        (output/'run_metadata.tsv').write_bytes(metadata)
        table = list(csv.DictReader(io.StringIO(metadata.decode()),delimiter='\t'))
        if len(table)!=1 or table[0]['run_accession']!=run:
            raise ValueError('unexpected ENA metadata grain')
        row=table[0]
        if ';' in row['fastq_ftp']:
            raise ValueError('single-file ONT run expected')
        expected_bytes=int(row['fastq_bytes'])
        if expected_bytes > MAX_COMPRESSED_BYTES:
            raise ValueError('run exceeds compressed-data budget')
        receipt['ena']=row
        ref_url='https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr12;start=111449213;end=111618315'
        ref_data=get(ref_url)
        (output/'human_reference_source.json').write_bytes(ref_data)
        dna=json.loads(ref_data)['dna'].upper()
        if hashlib.sha256(dna.encode()).hexdigest()!=NATIVE_SHA:
            raise ValueError('human reference DNA checksum mismatch')
        local=mp.revcomp(dna[146000:155000])
        fasta='>human_ATXN2_GRCh38_local\n'+local+'\n'
        lookup_url='https://rest.ensembl.org/lookup/symbol/mus_musculus/Atxn2?content-type=application/json'
        lookup_data=get(lookup_url)
        (output/'mouse_reference_lookup.json').write_bytes(lookup_data)
        lookup=json.loads(lookup_data)
        mouse_url='https://rest.ensembl.org/sequence/id/'+urllib.parse.quote(lookup['id'])+'?type=genomic;content-type=text/x-fasta;expand_5prime=5000;expand_3prime=5000'
        mouse=get(mouse_url)
        (output/'mouse_reference_original.fa').write_bytes(mouse)
        if not mouse.startswith(b'>'):
            raise ValueError('mouse decoy is not FASTA')
        mouse_seq=''.join(mouse.decode().splitlines()[1:]).upper()
        fasta+='>mouse_Atxn2_decoy\n'+mouse_seq+'\n'
        refpath=output/'restricted_reference.fa'
        refpath.write_text(fasta)
        receipt['references']={'human_url':ref_url,'human_dna_sha256':NATIVE_SHA,
            'mouse_lookup_url':lookup_url,'mouse_sequence_url':mouse_url,
            'mouse_assembly':lookup.get('assembly_name'),'mouse_dna_sha256':hashlib.sha256(mouse_seq.encode()).hexdigest(),
            'restricted_reference_sha256':hashlib.sha256(fasta.encode()).hexdigest()}
        receipt['mapping']={'preset':'map-ont','best_n':5,'n_threads':1,'human_local_cut_boundaries0':[4916,5383],
            'candidate_rule':'at least one returned human alignment with >=80 matching bases and reference overlap with [3800,6800)',
            'limitations':['not whole-genome alignment','native human reference has 23 repeat codons, not a verified 72Q haplotype','no AAV/vector reference','not validated structural-variant classification']}
        rawpath=output.parent/(run+'.fastq.gz.partial')
        raw_url='https://'+row['fastq_ftp']
        md5=hashlib.md5();sha=hashlib.sha256();nbytes=0
        with urllib.request.urlopen(urllib.request.Request(raw_url,headers={'User-Agent':'EditWitness-research-audit/0.0.2'}),timeout=180) as response,rawpath.open('wb') as target:
            while True:
                block=response.read(8*1024*1024)
                if not block:break
                nbytes+=len(block)
                if nbytes>MAX_COMPRESSED_BYTES:raise ValueError('download exceeds hard byte cap')
                md5.update(block);sha.update(block);target.write(block)
        if nbytes!=expected_bytes or md5.hexdigest()!=row['fastq_md5']:
            raise ValueError('ENA byte-count or MD5 mismatch; analysis refused')
        receipt['download']={'url':raw_url,'bytes':nbytes,'md5':md5.hexdigest(),'sha256':sha.hexdigest(),'validated':True}
        print('VERIFIED_DOWNLOAD',run,nbytes,flush=True)
        aligner=mp.Aligner(str(refpath),preset='map-ont',best_n=5,n_threads=1)
        if not aligner:raise RuntimeError('could not build restricted reference index')
        hist=Counter();nreads=nbases=selected=0
        pafpath=output/'target_candidates.jsonl.gz'
        fastqpath=output/'target_candidates.fastq.gz'
        with gzip.open(pafpath,'wt',compresslevel=6) as records,gzip.open(fastqpath,'wt',compresslevel=6) as fastq:
            for name,sequence,quality,comment in mp.fastx_read(str(rawpath),read_comment=True):
                if quality is None or len(sequence)!=len(quality):raise ValueError('invalid FASTQ record')
                nreads+=1;nbases+=len(sequence);hist[len(sequence)]+=1
                hits=list(aligner.map(sequence))
                if not any(h.ctg=='human_ATXN2_GRCh38_local' and h.mlen>=80 and h.r_st<6800 and h.r_en>3800 for h in hits):continue
                selected+=1
                alignments=[{'contig':h.ctg,'r_start':h.r_st,'r_end':h.r_en,'q_start':h.q_st,'q_end':h.q_en,
                    'strand':h.strand,'mapq':h.mapq,'match_bases':h.mlen,'block_length':h.blen,
                    'NM':h.NM,'is_primary':h.is_primary,'cigar':h.cigar_str} for h in hits]
                records.write(json.dumps({'read_id':name,'comment':comment,'length':len(sequence),'alignments':alignments},separators=(',',':'))+'\n')
                fastq.write('@'+name+(' '+comment if comment else '')+'\n'+sequence+'\n+\n'+quality+'\n')
        if nreads!=int(row['read_count']) or nbases!=int(row['base_count']):
            raise ValueError('parsed FASTQ counts do not match ENA; do not treat output as complete')
        with (output/'read_length_histogram.tsv').open('w') as f:
            f.write('length\tcount\n')
            for length,count in sorted(hist.items()):f.write(f'{length}\t{count}\n')
        receipt.update(status='complete',reads_processed=nreads,bases_processed=nbases,target_candidates=selected,
            wall_seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        rawpath.unlink()
    except Exception as error:
        receipt.update(status='failed',error=repr(error),wall_seconds=time.monotonic()-start)
        raise
    finally:
        receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
        (output/'receipt.json').write_text(json.dumps(receipt,indent=2))
        (output/'SHA256SUMS').write_text('\n'.join(hashlib.sha256(f.read_bytes()).hexdigest()+'  '+f.name for f in sorted(output.iterdir()) if f.is_file() and f.name!='SHA256SUMS')+'\n')
        print(json.dumps(receipt),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run',choices=sorted(RUNS));parser.add_argument('--output',type=pathlib.Path,required=True)
    args=parser.parse_args();collect(args.run,args.output)
