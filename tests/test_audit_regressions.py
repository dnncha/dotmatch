import copy

import pytest
from pydantic import ValidationError

from editwitness import analyze
from editwitness.compare import compare_models
from editwitness.generate import expand_deletions
from editwitness.io import InputError
from editwitness.models import Manifest
from editwitness.scan import scan_deletions


def model(manifest, **changes):
    data = manifest.model_dump(mode='json')
    data.update(changes)
    return Manifest.model_validate(data)


def test_alias_hypotheses_are_not_counterexamples(demo):
    data = demo.model_dump(mode='json')
    alias = copy.deepcopy(data['alleles'][1])
    alias['id'] = 'identical_intended'
    data['alleles'].append(alias)
    data['hypotheses'] = data['hypotheses'][:1]+[{'id':'alias', 'alleles':['identical_intended', 'intended']}]
    result = analyze(Manifest.model_validate(data))
    assert result.witnesses == ()
    assert result.distinct_alternatives == 0
    assert result.conclusion == 'no_distinct_alternatives'
    assert result.hypotheses[1].same_local_genotype_as_expected


def test_baseline_without_positive_expected_signal_is_not_reassuring(demo):
    data = demo.model_dump(mode='json')
    data['hypotheses'][0]['alleles'] = ['window_deleted', 'window_deleted']
    assert analyze(Manifest.model_validate(data)).conclusion == 'baseline_uninformative'


@pytest.mark.parametrize('function', [analyze, scan_deletions, expand_deletions, compare_models])
def test_unchecked_api_model_updates_are_revalidated(function, demo):
    malformed = demo.model_copy(update={'expected_hypothesis':'does_not_exist'})
    with pytest.raises(ValidationError):
        function(malformed)


def test_no_implicit_migration_of_observation_model(demo):
    raw = demo.model_dump(mode='json')
    raw.pop('observation_model')
    assert Manifest.model_validate(raw).observation_model == 'original-sites-presence-v1'


def test_exact_model_and_geometry_scan_report_different_model_identifiers(demo):
    exact = model(demo, observation_model='exact-local-sites-presence-v2')
    assert analyze(exact).model_version == 'exact-local-sites-presence-v2'
    assert scan_deletions(exact).model_version == 'original-sites-presence-v1'


def test_broad_replacement_changes_legacy_result_but_not_exact(demo):
    data = demo.model_dump(mode='json')
    intended = data['alleles'][1]['edits'][0]
    start, end = intended['start'], intended['end']
    ref = data['reference']['sequence']
    data['alleles'].append({'id':'broad_intended', 'edits':[{'start':100, 'end':end,
                           'sequence':ref[100:start]+intended['sequence']}]})
    data['hypotheses'].append({'id':'broad_plus_reference', 'alleles':['broad_intended','reference']})
    result = compare_models(Manifest.model_validate(data))
    assert result['models'][0]['model_version'] != result['models'][1]['model_version']
    assert result['input_manifest_sha256']
    # The comparison is evidence of response-model sensitivity, not empirical validation.
    assert 'not validation' in result['caveat']


def test_generation_complete_reproducible_and_explicit(demo):
    grid = dict(start_min=180, start_max=220, end_min=480, end_max=500, step=20)
    source = model(demo, observation_model='exact-local-sites-presence-v2', deletion_scan=grid)
    before = source.model_dump(mode='json')
    generated = expand_deletions(source)
    assert generated.generation.enumerated_deletions == 6
    assert generated.generation.added_alleles == 5
    assert generated.generation.duplicate_local_sequences == 1  # an existing declared alternative
    assert generated.generation.paired_with_expected_allele == 'intended'
    assert len(generated.hypotheses) == len(source.hypotheses)+5
    assert generated == expand_deletions(source)
    assert source.model_dump(mode='json') == before
    assert analyze(generated).generation == generated.generation
    with pytest.raises(InputError, match='already'):
        expand_deletions(generated)


def test_generation_refuses_large_grid_no_truncation(demo):
    with pytest.raises(InputError, match='128 alleles'):
        expand_deletions(demo)
    huge = model(demo, deletion_scan=dict(start_min=0,start_max=200,end_min=400,end_max=600))
    with pytest.raises(InputError, match='5,000'):
        expand_deletions(huge)
    with pytest.raises(InputError, match='deletion_scan'):
        expand_deletions(model(demo, deletion_scan=None))


def test_generation_refuses_silent_phase_choice(demo):
    data = demo.model_dump(mode='json')
    data['hypotheses'][0]['alleles'] = ['intended', 'reference']
    with pytest.raises(InputError, match='homozygous'):
        expand_deletions(Manifest.model_validate(data))


def test_generation_no_valid_events_is_explicit(demo):
    grid = dict(start_min=100,start_max=100,end_min=50,end_max=50)
    with pytest.raises(InputError, match='no valid'):
        expand_deletions(model(demo, deletion_scan=grid))


def test_multiple_products_are_used_by_hypothesis_engine(demo):
    data = demo.model_dump(mode='json')
    data['observation_model'] = 'exact-local-sites-presence-v2'
    ref = data['reference']['sequence']
    data['alleles'].append({'id':'duplicate_site','edits':[{'start':300,'end':300,'sequence':ref[200:220]}]})
    data['hypotheses'].append({'id':'multi','alleles':['duplicate_site','duplicate_site']})
    result = analyze(Manifest.model_validate(data))
    obs = next(o for o in result.hypothesis_observations if o.hypothesis_id=='multi' and o.assay_id=='inner')
    assert len(obs.signal_ids) == 2
    assert 'MULTIPLE_LOCAL_PRODUCTS' in {n.code for n in result.notices}


def test_duplicate_alternative_genotypes_do_not_inflate_evidence(demo):
    data = demo.model_dump(mode='json')
    data['hypotheses'].append({'id':'zzz_alias', 'alleles':list(reversed(data['hypotheses'][2]['alleles']))})
    before, after = analyze(demo), analyze(Manifest.model_validate(data))
    assert before.distinct_alternatives == after.distinct_alternatives
    assert before.witnesses == after.witnesses
    assert 'DUPLICATE_LOCAL_GENOTYPES' in {n.code for n in after.notices}


def test_comparison_exposes_representation_dependent_dropout(demo):
    data = demo.model_dump(mode='json')
    ref = data['reference']['sequence']
    alt = next(base for base in 'ACGT' if base != ref[750])
    data['alleles'].append({'id':'broad_distal_change','edits':[
        {'start':100,'end':751,'sequence':ref[100:750]+alt}]})
    data['hypotheses'].append({'id':'distal','alleles':['intended','broad_distal_change']})
    compared = compare_models(Manifest.model_validate(data))
    assert 'distal' in compared['original_only']
    assert compared['witnesses_changed'] is True


def test_report_does_not_mislabel_multiple_signals_as_absence(demo):
    from editwitness.report import render_report
    data = demo.model_dump(mode='json')
    data['observation_model'] = 'exact-local-sites-presence-v2'
    ref = data['reference']['sequence']
    data['alleles'][1]['edits'].insert(0, {'start':300,'end':300,'sequence':ref[200:220]})
    result = analyze(Manifest.model_validate(data))
    report = render_report(result)
    assert 'product 2, read 1:' in report
    assert 'Inspect the genomic alternatives' in report
