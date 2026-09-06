import json
from pathlib import Path

from editwitness.cli import main
from editwitness.models import Manifest


def test_new_cli_workflow(tmp_path, capsys, demo):
    source = tmp_path/'source.json'
    data = demo.model_dump(mode='json')
    data['observation_model'] = 'exact-local-sites-presence-v2'
    data['deletion_scan'] = dict(start_min=180,start_max=220,end_min=480,end_max=500,step=20)
    source.write_text(json.dumps(data))
    expanded = tmp_path/'expanded.json'
    assert main(['expand-deletions',str(source),'-o',str(expanded)]) == 0
    assert Manifest.model_validate_json(expanded.read_text()).generation.enumerated_deletions == 6
    assert main(['compare-models',str(expanded)]) == 0
    assert json.loads(capsys.readouterr().out)['kind'] == 'editwitness.model_comparison'
    assert main(['analyze',str(expanded),'--compact']) == 0
    assert json.loads(capsys.readouterr().out)['distinct_alternatives'] > 0


def test_demo_model_selection_is_explicit(capsys):
    assert main(['demo']) == 0
    assert json.loads(capsys.readouterr().out)['observation_model'] == 'exact-local-sites-presence-v2'
    assert main(['demo','--legacy-model']) == 0
    assert json.loads(capsys.readouterr().out)['observation_model'] == 'original-sites-presence-v1'


def test_init_adds_a_bounded_grid(demo, tmp_path, capsys):
    from editwitness.sequence import reverse_complement
    fasta = tmp_path/'ref.fasta'
    fasta.write_text('>test\n'+demo.reference.sequence+'\n')
    alt = next(b for b in 'ACGT' if b != demo.reference.sequence[450])
    args = ['init','--fasta',str(fasta),'--left-primer',demo.reference.sequence[200:220],
            '--right-primer',reverse_complement(demo.reference.sequence[680:700]),
            '--edit-position','450','--alternate',alt,'--deletion-radius','300','--deletion-step','50']
    assert main(args) == 0
    generated = json.loads(capsys.readouterr().out)
    assert generated['deletion_scan']['start_min'] == 150
    assert generated['observation_model'] == 'exact-local-sites-presence-v2'
    assert main(args[:-1]+['0']) == 2
    assert json.loads(capsys.readouterr().err)['kind'] == 'editwitness.error'


def test_old_evidence_requires_producing_version(tmp_path, capsys):
    path = tmp_path/'old.json'
    path.write_text('{"kind":"editwitness.analysis","schema_version":"1.0"}')
    assert main(['verify',str(path)]) == 2
    assert 'producing package version' in json.loads(capsys.readouterr().err)['message']
