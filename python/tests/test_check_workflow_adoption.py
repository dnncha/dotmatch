import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_workflow_adoption.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_workflow_adoption", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _manifest(status: str = "not_ready", integrations=None) -> dict:
    manifest = {
        "schema_version": 1,
        "status": status,
        "integrations": [] if integrations is None else integrations,
    }
    if status == "not_ready":
        manifest["blockers"] = ["No external workflow adoption has been recorded yet."]
        manifest["next_action"] = "Land an external workflow integration and record stable public links."
    return manifest


def _integration(**overrides) -> dict:
    item = {
        "id": "nfcore_dotmatch_crispr_count",
        "type": "nf_core_module",
        "name": "nf-core dotmatch/crispr_count module",
        "status": "accepted",
        "adoption_url": "https://github.com/nf-core/modules/tree/master/modules/nf-core/dotmatch/crispr_count",
        "evidence_url": "https://github.com/nf-core/modules/pull/12345",
        "validation_notes": "Accepted upstream module with tests and version reporting.",
        "recorded_date": "2026-05-07",
    }
    item.update(overrides)
    return item


def _write_manifest(root: Path, manifest: dict) -> None:
    path = root / "docs" / "workflow-adoption.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def test_workflow_adoption_accepts_ready_manifest_with_reachable_external_links(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_manifest(tmp_path, _manifest(status="ready", integrations=[_integration()],))
    monkeypatch.setattr(checker, "url_ok", lambda url: url.startswith("https://github.com/nf-core/modules"))

    result = checker.audit(tmp_path)

    assert result.failures == []
    assert any("external workflow adoption recorded" in item for item in result.passed)


def test_workflow_adoption_rejects_not_ready_without_external_records(tmp_path):
    checker = _load_checker()
    _write_manifest(tmp_path, _manifest())

    result = checker.audit(tmp_path)

    assert any("external workflow adoption is not recorded" in failure for failure in result.failures)


def test_workflow_adoption_rejects_ready_status_without_integrations(tmp_path):
    checker = _load_checker()
    _write_manifest(tmp_path, _manifest(status="ready", integrations=[]))

    result = checker.audit(tmp_path)

    assert any("ready status requires at least one external integration" in failure for failure in result.failures)


def test_workflow_adoption_rejects_ready_status_with_stale_blockers(tmp_path, monkeypatch):
    checker = _load_checker()
    manifest = _manifest(status="ready", integrations=[_integration()])
    manifest["blockers"] = ["No external adoption yet."]
    manifest["next_action"] = "Keep working."
    _write_manifest(tmp_path, manifest)
    monkeypatch.setattr(checker, "url_ok", lambda url: True)

    result = checker.audit(tmp_path)

    assert any("ready workflow adoption manifest must not declare blockers" in failure for failure in result.failures)
    assert any("ready workflow adoption manifest must not declare next_action" in failure for failure in result.failures)


def test_workflow_adoption_rejects_unreachable_external_links(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_manifest(tmp_path, _manifest(status="ready", integrations=[_integration()]))
    monkeypatch.setattr(checker, "url_ok", lambda url: False)

    result = checker.audit(tmp_path)

    assert any("adoption_url is not reachable" in failure for failure in result.failures)
    assert any("evidence_url is not reachable" in failure for failure in result.failures)


def test_workflow_adoption_rejects_unsupported_integration_type(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_manifest(tmp_path, _manifest(status="ready", integrations=[_integration(type="blog_post")]))
    monkeypatch.setattr(checker, "url_ok", lambda url: True)

    result = checker.audit(tmp_path)

    assert any("unsupported integration type" in failure for failure in result.failures)


def test_workflow_adoption_rejects_duplicate_integration_ids(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_manifest(tmp_path, _manifest(status="ready", integrations=[_integration(), _integration()]))
    monkeypatch.setattr(checker, "url_ok", lambda url: True)

    result = checker.audit(tmp_path)

    assert any("duplicate workflow adoption integration id: nfcore_dotmatch_crispr_count" in failure for failure in result.failures)


def test_workflow_adoption_rejects_placeholder_urls(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_manifest(
        tmp_path,
        _manifest(
            status="ready",
            integrations=[
                _integration(
                    adoption_url="https://example.org/stable-adoption-page",
                    evidence_url="https://example.org/review-or-release-record",
                )
            ],
        ),
    )
    monkeypatch.setattr(checker, "url_ok", lambda url: True)

    result = checker.audit(tmp_path)

    assert any("adoption_url must not use placeholder domains" in failure for failure in result.failures)
    assert any("evidence_url must not use placeholder domains" in failure for failure in result.failures)


def test_workflow_adoption_rejects_invalid_status_and_recorded_date(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_manifest(
        tmp_path,
        _manifest(
            status="ready",
            integrations=[_integration(status="draft", recorded_date="May 7 2026")],
        ),
    )
    monkeypatch.setattr(checker, "url_ok", lambda url: True)

    result = checker.audit(tmp_path)

    assert any("has invalid integration status: draft" in failure for failure in result.failures)
    assert any("recorded_date must be YYYY-MM-DD" in failure for failure in result.failures)
