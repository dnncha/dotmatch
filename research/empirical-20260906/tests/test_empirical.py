import importlib.util, json, os
from pathlib import Path
import pytest

CODE_ROOT=Path(__file__).resolve().parents[1]
ROOT=Path(os.environ.get('EW_AUDIT_DATA_ROOT',str(CODE_ROOT)))
spec=importlib.util.spec_from_file_location('empirical',CODE_ROOT/'scripts/empirical.py')
e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)

@pytest.mark.parametrize('s',['ACGT','AAAACCCG','CGTNN','T','ACGT'*30])
def test_reverse_complement_involution(s): assert e.reverse_complement(e.reverse_complement(s))==s
@pytest.mark.parametrize('s',['','ACGX','augt','123'])
def test_bad_reference_rejected(s):
    with pytest.raises(ValueError):e.reverse_complement(s)

def test_overlapping_oligo_matches():assert e.sites('AAAAA','AAA')==[0,1,2]
@pytest.mark.parametrize('s',['','NN','AGX'])
def test_bad_oligo_rejected(s):
    with pytest.raises(ValueError):e.sites('ACGT',s)

def test_enumerates_multiple_products():
    found=e.products('AAAGGGAAAGGGTTTCCC','AAA','GGG')
    assert len(found)==4
    assert sorted(p['orientation'] for p in found)==['+','+','-','-']

def test_orientation_equivariance():
    seq='AAACCTGACGTGTCG';f='AAAC';r='CGAC'
    a=e.products(seq,f,r);b=e.products(e.reverse_complement(seq),f,r)
    assert len(a)==len(b)==1
    assert a[0]['length']==b[0]['length']
    assert a[0]['orientation']!=b[0]['orientation']

@pytest.mark.parametrize('reference_start,cigar,reference_end,blocks,gaps',[
 (100,'20M5I30M',150,[(100,120),(120,150)],[]),
 (100,'20M5D30M',155,[(100,120),(125,155)],[(120,125)]),
 (100,'10=2X5N3M',120,[(100,110),(110,112),(117,120)],[]),
 (100,'5S20M3H',120,[(100,120)],[]),
])
def test_cigar_reference_accounting(reference_start,cigar,reference_end,blocks,gaps):
    assert e.cigar_blocks({'r_start':reference_start,'r_end':reference_end,'cigar':cigar})==(blocks,gaps)
@pytest.mark.parametrize('cigar',['','10Z','0M','M10','10Mbad'])
def test_bad_cigar_refused(cigar):
    with pytest.raises(ValueError):e.cigar_blocks({'r_start':0,'r_end':10,'cigar':cigar})
def test_inconsistent_end_refused():
    with pytest.raises(ValueError):e.cigar_blocks({'r_start':0,'r_end':11,'cigar':'10M'})

def record(cigar='416M467D217M',start=4500,end=5600,mapq=60,primary=True,contig=e.HUMAN):
    return {'read_id':'test','alignments':[{'contig':contig,'is_primary':primary,'mapq':mapq,
      'match_bases':633,'block_length':1100,'r_start':start,'r_end':end,'strand':1,'cigar':cigar}]}

def test_nominal_deletion_pattern():
    c=e.classify_read(record());assert c['status']=='eligible_two_flank_read';assert c['deletion_compatible']
def test_non_deletion_is_eligible():
    c=e.classify_read(record('1100M'));assert c['status']=='eligible_two_flank_read';assert not c['deletion_compatible']
def test_reference_skip_is_not_deletion():
    assert not e.classify_read(record('416M467N217M'))['deletion_compatible']
def test_single_anchor_excluded():
    assert e.classify_read(record('633M',end=5133))['status']=='insufficient_flanking_alignment'
def test_low_mapq_excluded():assert e.classify_read(record(mapq=19))['status']=='below_mapq'
def test_secondary_excluded():assert e.classify_read(record(primary=False))['status']=='no_primary_human_hit'
def test_mouse_decoy_excluded():assert e.classify_read(record(contig='mouse_Atxn2_decoy'))['status']=='no_primary_human_hit'
def test_split_alignments_not_double_counted():
    r=record();r['alignments'].append(dict(r['alignments'][0]));assert e.classify_read(r)['deletion_compatible']
def test_endpoint_tolerance_not_length_only():
    assert not e.classify_read(record('516M467D117M'))['deletion_compatible']
def test_m_coverage_is_not_identity():
    assert e.anchor_coverage([(1,5),(7,10)],(3,9))==4

def test_real_geometry_lengths():
    result=e.audit_geometry(ROOT)
    lookup={(r['state'],r['assay']):r for r in result['products']}
    assert lookup[('nominal_72_codon_model','outer_PCR')]['product_lengths_bp']=='2317'
    assert lookup[('nominal_precise_dual_cut_deletion','outer_PCR')]['product_lengths_bp']=='1703'
    assert lookup[('nominal_precise_dual_cut_deletion','g4_local_ddPCR')]['exact_products']==0
    assert lookup[('nominal_precise_dual_cut_deletion','g5_local_ddPCR')]['exact_products']==0
    assert lookup[('native_23_codon_reference','g4_local_ddPCR')]['product_lengths_bp']=='142'
    assert lookup[('native_23_codon_reference','g5_local_ddPCR')]['product_lengths_bp']=='148'

def test_real_metadata_grain():
    result=e.audit_metadata(ROOT)
    assert result['runs']==42;assert result['biosamples']==15;assert result['same_biosample_pairs']==3
    assert not result['original_molecule_inference_eligible']
def test_published_categories_not_double_counted():
    r=e.audit_published_values(ROOT)
    row=next(x for x in r if x['excel_row']==51)
    assert row['exclusive_mutation_sum_percent']==pytest.approx(14.77)
    assert row['integration_annotation_percent']==pytest.approx(2.02)
def test_false_pdf_is_rejected():
    checks=e.verify_receipts(ROOT)
    wrong=next(c for c in checks if c['path'].endswith('assay_sources/Simpson-supplement.pdf'))
    assert wrong['hash_ok'];assert not wrong['content_valid']
def test_calibration_does_not_claim_validation():
    result=e.audit_calibration(ROOT)
    assert not result['full_original_molecule_calibration'];assert len(result['classes'])==9
