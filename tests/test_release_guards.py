"""Regression guards for publication-facing errors and bounded evidence views."""
import copy

import pytest

from editwitness import analyze
from editwitness.generate import expand_deletions
from editwitness.io import InputError
from editwitness.models import Manifest


def generated_input(demo):
    data = demo.model_dump(mode='json')
    data['deletion_scan'] = dict(start_min=181, start_max=181, end_min=480, end_max=480)
    return data


def test_generation_refuses_allele_id_collision(demo):
    data = generated_input(demo)
    data['alleles'][0]['id'] = 'del_181_480'
    for h in data['hypotheses']:
        h['alleles'] = ['del_181_480' if a == 'reference' else a for a in h['alleles']]
    with pytest.raises(InputError, match='allele id collides'):
        expand_deletions(Manifest.model_validate(data))


def test_generation_refuses_hypothesis_id_collision(demo):
    data = generated_input(demo)
    data['hypotheses'][1]['id'] = 'expected_plus_del_181_480'
    with pytest.raises(InputError, match='hypothesis id collides'):
        expand_deletions(Manifest.model_validate(data))


def test_generation_reuses_existing_sequence_without_inventing_duplicate_allele(demo):
    data = generated_input(demo)
    data['alleles'].append({'id':'existing_deletion','edits':[{'start':181,'end':480,'sequence':''}]})
    expanded = expand_deletions(Manifest.model_validate(data))
    assert expanded.generation.added_alleles == 0
    assert expanded.hypotheses[-1].alleles == ('intended','existing_deletion')


def test_generation_rejects_hypothesis_overflow(demo):
    data = generated_input(demo)
    data['hypotheses'] = [{'id':f'h{i}', 'alleles':['intended','intended']} for i in range(1000)]
    data['expected_hypothesis'] = 'h0'
    with pytest.raises(InputError, match='1,000 hypotheses'):
        expand_deletions(Manifest.model_validate(data))


def test_public_api_does_not_treat_dictionary_as_validated_model(demo):
    with pytest.raises(TypeError, match='Manifest'):
        analyze(demo.model_dump(mode='json'))


def test_witness_cli_resolves_alternative_alias(demo, tmp_path, capsys):
    import json
    from editwitness.cli import main
    data = demo.model_dump(mode='json')
    alias = copy.deepcopy(next(h for h in data['hypotheses'] if h['id']=='hidden_primer_deletion'))
    alias['id'] = 'zzz_alias'
    data['hypotheses'].append(alias)
    source = tmp_path/'input.json'
    source.write_text(json.dumps(data),encoding='utf-8')
    assert main(['witness',str(source),'--hypothesis','zzz_alias']) == 0
    witness = json.loads(capsys.readouterr().out)
    assert witness['requested_hypothesis'] == 'zzz_alias'
    assert witness['witness']['hypothesis_id'] == 'hidden_primer_deletion'
