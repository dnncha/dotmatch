import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_publication_ready.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_publication_ready", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_minimal_repo(root: Path) -> None:
    files = {
        "README.md": "# DotMatch\n\n`v0.1.0` includes stable launch artifacts.\n\nSee [Scientific claim ledger](docs/scientific-claims.md).\n",
        "CHANGELOG.md": "# Changelog\n\n## 0.1.0\n\n- Initial launch.\n",
        "LICENSE": "Apache-2.0\n",
        "CITATION.cff": "cff-version: 1.2.0\ntitle: DotMatch\nversion: \"0.1.0\"\n",
        "codemeta.json": (
            '{"name": "DotMatch", "license": "https://spdx.org/licenses/Apache-2.0", '
            '"version": "0.1.0", "keywords": ["bioinformatics"]}\n'
        ),
        "CONTRIBUTING.md": "# Contributing\n",
        "CODE_OF_CONDUCT.md": "# Code of Conduct\n",
        "SECURITY.md": "# Security\n",
        "SUPPORT.md": "# Support\n",
        "pyproject.toml": "[project]\nname = \"dotmatch\"\nversion = \"0.1.0\"\nlicense = \"Apache-2.0\"\n",
        "package.json": '{"version": "0.1.0", "license": "Apache-2.0"}\n',
        "MANIFEST.in": "include src/qdalign.c\ninclude include/qdalign.h\n",
        "setup.py": "from setuptools import setup\nsetup()\n",
        ".gitignore": "build/\n.DS_Store\n*.tsbuildinfo\n",
        ".gitattributes": "* text=auto eol=lf\n",
        ".editorconfig": "root = true\n",
        ".github/workflows/ci.yml": "name: ci\n",
        ".github/workflows/codeql.yml": "name: codeql\n",
        ".github/workflows/release.yml": "name: release\n",
        ".github/PULL_REQUEST_TEMPLATE.md": "## Evidence\n",
        ".github/ISSUE_TEMPLATE/bug_report.yml": "name: Bug report\n",
        ".github/ISSUE_TEMPLATE/feature_request.yml": "name: Feature request\n",
        ".github/ISSUE_TEMPLATE/benchmark_claim.yml": "name: Benchmark\n",
        "docs/scientific-claims.md": (
            "# Claims\n\n"
            "`make barcode-sota-gate` currently fails.\n"
            "`make bcl-sota-gate` currently fails.\n"
            "General aligner replacement is blocked.\n"
        ),
        "docs/github-launch-checklist.md": "# Launch\n",
        "docs/release-process.md": "# Release Process\n",
        "docs/methods-and-citation.md": "# Methods\n",
        "docs/open-core-strategy.md": "# Open Core\n",
        "docs/packaging.md": "# Packaging\n",
        "src/qdalign.c": "int x;\n",
        "include/qdalign.h": "int x;\n",
    }
    for path, text in files.items():
        full = root / path
        full.parent.mkdir(parents=True, exist_ok=True)
        full.write_text(text, encoding="utf-8")


def test_publication_ready_accepts_minimal_valid_repo(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)

    result = checker.audit(tmp_path)

    assert result.failures == []
    assert any("required files" in item for item in result.passed)


def test_publication_ready_reports_missing_claim_ledger(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "docs" / "scientific-claims.md").unlink()

    result = checker.audit(tmp_path)

    assert any("docs/scientific-claims.md" in failure for failure in result.failures)


def test_publication_ready_reports_missing_release_workflow(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / ".github" / "workflows" / "release.yml").unlink()

    result = checker.audit(tmp_path)

    assert any(".github/workflows/release.yml" in failure for failure in result.failures)


def test_publication_ready_reports_missing_codeql_workflow(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / ".github" / "workflows" / "codeql.yml").unlink()

    result = checker.audit(tmp_path)

    assert any(".github/workflows/codeql.yml" in failure for failure in result.failures)


def test_publication_ready_reports_missing_changelog(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "CHANGELOG.md").unlink()

    result = checker.audit(tmp_path)

    assert any("CHANGELOG.md" in failure for failure in result.failures)


def test_publication_ready_rejects_unignored_large_or_macos_files(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / ".DS_Store").write_text("metadata", encoding="utf-8")
    large = tmp_path / "benchmarks" / "raw" / "too_large.csv"
    large.parent.mkdir(parents=True, exist_ok=True)
    large.write_bytes(b"x" * 6)
    monkeypatch.setattr(checker, "MAX_FILE_BYTES", 5)

    result = checker.audit(tmp_path)

    assert any(".DS_Store" in failure for failure in result.failures)
    assert any("too_large.csv" in failure for failure in result.failures)


def test_publication_ready_requires_blocked_sota_status(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "docs" / "scientific-claims.md").write_text("# Claims\n", encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("barcode-sota-gate" in failure for failure in result.failures)
    assert any("bcl-sota-gate" in failure for failure in result.failures)


def test_publication_ready_rejects_mismatched_release_versions(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "CITATION.cff").write_text(
        "cff-version: 1.2.0\ntitle: DotMatch\nversion: \"0.1.0-dev\"\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("version mismatch" in failure for failure in result.failures)


def test_publication_ready_rejects_dev_readme_status(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    readme = (tmp_path / "README.md").read_text(encoding="utf-8")
    (tmp_path / "README.md").write_text(readme.replace("v0.1.0", "v0.1.0-dev"), encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("README.md must not describe the launch version as dev" in failure for failure in result.failures)


def test_publication_ready_rejects_local_absolute_paths(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    raw = tmp_path / "benchmarks" / "raw" / "example.csv"
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_text("tool,command\nx,/" + "Users/alice/projects/dotmatch/dotmatch count\n", encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("local absolute path" in failure and "example.csv" in failure for failure in result.failures)
