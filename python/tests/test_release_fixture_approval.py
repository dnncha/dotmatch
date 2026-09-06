import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def checker():
    spec = importlib.util.spec_from_file_location('repository_ready_fixture_check', ROOT / 'scripts/check_repository_ready.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_approved_synthetic_fixture_requires_exact_reviewed_bytes(tmp_path):
    module = checker()
    relative = 'examples/assignment_sensitivity/reads.fastq'
    target = tmp_path / relative
    target.parent.mkdir(parents=True)
    target.write_bytes((ROOT / relative).read_bytes())
    result = module.AuditResult()
    module.check_no_private_raw_data(tmp_path, result)
    assert result.ok, result.failures
    target.write_bytes(target.read_bytes() + b'\n')
    result = module.AuditResult()
    module.check_no_private_raw_data(tmp_path, result)
    assert not result.ok and any(relative in error for error in result.failures)


def test_synthetic_approval_does_not_allow_other_files_in_directory(tmp_path):
    module = checker()
    target = tmp_path / 'examples/assignment_sensitivity/unreviewed.fastq'
    target.parent.mkdir(parents=True)
    target.write_text('@r\nACGT\n+\nIIII\n')
    result = module.AuditResult()
    module.check_no_private_raw_data(tmp_path, result)
    assert not result.ok and any('unreviewed.fastq' in error for error in result.failures)
