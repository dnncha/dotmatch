import importlib.util
from pathlib import Path
import pytest

spec = importlib.util.spec_from_file_location('approved_release_status', Path(__file__).resolve().parents[2] / 'scripts/approved_release.py')
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)

BLOCKED_PREVIEW = {'context': 'Vercel', 'state': 'failure', 'description': 'Account is blocked.', 'target_url': 'https://vercel.com/knowledge/why-is-my-account-deployment-blocked'}


def test_exact_nonproduction_account_block_is_visible_and_pages_is_required(capsys):
    assert release.status_verdict({'statuses': [BLOCKED_PREVIEW]})
    assert 'Non-release preview unavailable' in capsys.readouterr().out
    assert '.github/workflows/pages.yml' in release.REQUIRED_WORKFLOWS


@pytest.mark.parametrize('changes', [
    {'description': 'Build failed.'},
    {'context': 'security'},
    {'state': 'error'},
    {'target_url': 'https://example.com/unrelated'},
])
def test_no_build_security_or_unknown_failure_is_exempt(changes):
    with pytest.raises(RuntimeError, match='Commit status failed'):
        release.status_verdict({'statuses': [dict(BLOCKED_PREVIEW, **changes)]})


def test_pending_checks_wait_and_successful_checks_pass():
    assert not release.status_verdict({'statuses': [{'context': 'security', 'state': 'pending'}]})
    assert release.status_verdict({'statuses': [{'context': 'security', 'state': 'success'}]})
    with pytest.raises(RuntimeError):
        release.status_verdict({'statuses': [BLOCKED_PREVIEW, {'context': 'security', 'state': 'failure'}]})
