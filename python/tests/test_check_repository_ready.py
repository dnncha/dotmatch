import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_repository_ready.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_repository_ready", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_minimal_repo(root: Path) -> None:
    files = {
        "README.md": (
            "# DotMatch\n\n"
            "`v0.1.0` includes stable release artifacts.\n\n"
            "See [Evidence Notes](docs/scientific-claims.md).\n\n"
            "Installation and distribution status is tracked in "
            "[Packaging Notes](docs/packaging.md) and "
            "[distribution release record](docs/distribution-release.json). "
            "PyPI, Bioconda, GHCR/BioContainers, and Zenodo release artifacts are not yet published. "
            "Until those channels are listed here as available, install from source.\n\n"
            "See [Release process](docs/release-process.md).\n"
        ),
        "CHANGELOG.md": "# Changelog\n\n## 0.1.0\n\n- Initial release.\n",
        "LICENSE": "Apache-2.0\n",
        "CITATION.cff": "cff-version: 1.2.0\ntitle: DotMatch\nversion: \"0.1.0\"\n",
        "codemeta.json": (
            '{"name": "DotMatch", "license": "https://spdx.org/licenses/Apache-2.0", '
            '"version": "0.1.0", "softwareVersion": "0.1.0", "keywords": ["bioinformatics"]}\n'
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
        ".zenodo.json": (
            '{"title": "DotMatch", "upload_type": "software", "version": "0.1.0", '
            '"license": "Apache-2.0", "keywords": ["bioinformatics"]}\n'
        ),
        ".github/workflows/ci.yml": "name: ci\n",
        ".github/workflows/codeql.yml": "name: codeql\n",
        ".github/workflows/release.yml": "name: release\n",
        ".github/PULL_REQUEST_TEMPLATE.md": (
            "## Evidence\n\n"
            "- [ ] `make test`\n"
            "- [ ] `make cli-test`\n"
            "- [ ] `make python-test`\n"
            "- [ ] `make pretag-ready` if release surfaces changed.\n\n"
            "## Claim Boundary\n\n"
            "- [ ] This PR does not broaden README/docs claims beyond checked evidence.\n"
        ),
        ".github/ISSUE_TEMPLATE/bug_report.yml": "name: Bug report\n",
        ".github/ISSUE_TEMPLATE/feature_request.yml": "name: Feature request\n",
        ".github/ISSUE_TEMPLATE/benchmark_evidence.yml": "name: Benchmark\n",
        "docs/scientific-claims.md": (
            "# Evidence Notes\n\n"
            "`docs/assay-evidence.json` tracks assay lanes.\n"
            "`make barcode-comparison-gate` requires real-data rows.\n"
            "`make bcl-comparison-gate` requires real run folders.\n"
            "General aligner replacement is blocked.\n"
        ),
        "docs/assay-evidence.json": '{"schema_version": 1, "assays": []}\n',
        "docs/distribution-release.json": '{"schema_version": 1, "status": "not_released", "channels": []}\n',
        "docs/workflow-adoption.json": '{"schema_version": 1, "status": "not_ready", "integrations": []}\n',
        "docs/release-process.md": "# Release Process\n",
        "docs/methods-and-citation.md": "# Methods\n",
        "docs/packaging.md": "# Packaging\n",
        "docs/native-comparator-scope.md": (
            "# Native Comparator Scope\n\n"
            "Current native comparator: Edlib exhaustive global edit-distance assignment with zero mismatches.\n\n"
            "SeqAn and Parasail are not part of the checked README, website, or release-note performance comparison set yet. "
            "They require equivalent global edit-distance or documented semi-global scoring semantics, fixed threshold k, "
            "fixed threshold `k`, native dependency/version capture, raw CSV rows, and zero assignment mismatches before claims change.\n"
            "Supported wording is limited to Edlib exhaustive global edit-distance assignment scans.\n"
        ),
        "docs/schemas.md": "# Schemas\n",
        "examples/workflows/galaxy/dotmatch_crispr_count.xml": "<tool id=\"dotmatch_crispr_count\" />\n",
        "examples/workflows/multiqc/multiqc_config.yaml": "custom_data:\n  dotmatch_sample_qc:\n",
        "examples/workflows/nf-core/README.md": "# nf-core-style Module Candidate\n",
        "examples/workflows/nf-core/modules/local/dotmatch/crispr_count/main.nf": "process DOTMATCH_CRISPR_COUNT\n",
        "examples/workflows/nf-core/modules/local/dotmatch/crispr_count/meta.yml": "name: dotmatch_crispr_count\n",
        "examples/workflows/nextflow/main.nf": "nextflow.enable.dsl=2\n",
        "examples/workflows/snakemake/Snakefile": "rule all:\n",
        "packaging/bioconda/meta.yaml": "package:\n  name: dotmatch\n",
        "packaging/bioconda/build.sh": "#!/usr/bin/env bash\n",
        "scripts/check_assay_evidence.py": "#!/usr/bin/env python3\n",
        "scripts/check_alphabet_policy.py": "#!/usr/bin/env python3\n",
        "scripts/check_citation_metadata.py": "#!/usr/bin/env python3\n",
        "scripts/check_distribution_channels.py": "#!/usr/bin/env python3\n",
        "scripts/check_distribution_record.py": "#!/usr/bin/env python3\n",
        "scripts/check_bioconda_recipe.py": "#!/usr/bin/env python3\n",
        "scripts/check_native_comparator_scope.py": "#!/usr/bin/env python3\n",
        "scripts/check_workflow_adoption.py": "#!/usr/bin/env python3\n",
        "src/qdalign.c": "int x;\n",
        "include/qdalign.h": "int x;\n",
    }
    for path, text in files.items():
        full = root / path
        full.parent.mkdir(parents=True, exist_ok=True)
        full.write_text(text, encoding="utf-8")


def test_repository_ready_accepts_minimal_valid_repo(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)

    result = checker.audit(tmp_path)

    assert result.failures == []
    assert any("required files" in item for item in result.passed)


def test_repository_ready_reports_missing_evidence_notes(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "docs" / "scientific-claims.md").unlink()

    result = checker.audit(tmp_path)

    assert any("docs/scientific-claims.md" in failure for failure in result.failures)


def test_repository_ready_reports_missing_release_process_link(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    readme = (tmp_path / "README.md").read_text(encoding="utf-8")
    (tmp_path / "README.md").write_text(
        readme.replace("See [Release process](docs/release-process.md).\n", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("docs/release-process.md" in failure for failure in result.failures)


def test_repository_ready_reports_missing_source_install_direction(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    readme = (tmp_path / "README.md").read_text(encoding="utf-8")
    (tmp_path / "README.md").write_text(
        readme.replace("Until those channels are listed here as available, install from source.", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("install from source" in failure for failure in result.failures)


def test_repository_ready_requires_distribution_status_in_readme(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "README.md").write_text(
        "# DotMatch\n\n"
        "`v0.1.0` includes stable release artifacts.\n\n"
        "See [Evidence Notes](docs/scientific-claims.md).\n\n"
        "See [Release process](docs/release-process.md).\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("docs/packaging.md" in failure for failure in result.failures)
    assert any("docs/distribution-release.json" in failure for failure in result.failures)
    assert any("public package availability" in failure for failure in result.failures)


def test_repository_ready_reports_missing_release_workflow(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / ".github" / "workflows" / "release.yml").unlink()

    result = checker.audit(tmp_path)

    assert any(".github/workflows/release.yml" in failure for failure in result.failures)


def test_repository_ready_reports_missing_codeql_workflow(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / ".github" / "workflows" / "codeql.yml").unlink()

    result = checker.audit(tmp_path)

    assert any(".github/workflows/codeql.yml" in failure for failure in result.failures)


def test_repository_ready_reports_incomplete_pull_request_template(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / ".github" / "PULL_REQUEST_TEMPLATE.md").write_text(
        "## Evidence\n\n- [ ] `make test`\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("make cli-test" in failure for failure in result.failures)
    assert any("make pretag-ready" in failure for failure in result.failures)
    assert any("claim-boundary" in failure for failure in result.failures)


def test_repository_ready_reports_missing_changelog(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "CHANGELOG.md").unlink()

    result = checker.audit(tmp_path)

    assert any("CHANGELOG.md" in failure for failure in result.failures)


def test_repository_ready_reports_missing_zenodo_metadata(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / ".zenodo.json").unlink()

    result = checker.audit(tmp_path)

    assert any(".zenodo.json" in failure for failure in result.failures)


def test_repository_ready_rejects_codemeta_without_software_version(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "codemeta.json").write_text(
        '{"name": "DotMatch", "license": "https://spdx.org/licenses/Apache-2.0", '
        '"version": "0.1.0", "keywords": ["bioinformatics"]}\n',
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("softwareVersion" in failure for failure in result.failures)


def test_repository_ready_reports_missing_workflow_and_bioconda_artifacts(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "examples" / "workflows" / "nf-core" / "modules" / "local" / "dotmatch" / "crispr_count" / "main.nf").unlink()
    (tmp_path / "examples" / "workflows" / "nextflow" / "main.nf").unlink()
    (tmp_path / "packaging" / "bioconda" / "meta.yaml").unlink()

    result = checker.audit(tmp_path)

    assert any("examples/workflows/nf-core/modules/local/dotmatch/crispr_count/main.nf" in failure for failure in result.failures)
    assert any("examples/workflows/nextflow/main.nf" in failure for failure in result.failures)
    assert any("packaging/bioconda/meta.yaml" in failure for failure in result.failures)


def test_repository_ready_reports_missing_distribution_channel_verifier(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "scripts" / "check_distribution_channels.py").unlink()

    result = checker.audit(tmp_path)

    assert any("scripts/check_distribution_channels.py" in failure for failure in result.failures)


def test_repository_ready_reports_missing_alphabet_policy_verifier(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "scripts" / "check_alphabet_policy.py").unlink()

    result = checker.audit(tmp_path)

    assert any("scripts/check_alphabet_policy.py" in failure for failure in result.failures)


def test_repository_ready_reports_missing_citation_metadata_verifier(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "scripts" / "check_citation_metadata.py").unlink()

    result = checker.audit(tmp_path)

    assert any("scripts/check_citation_metadata.py" in failure for failure in result.failures)


def test_repository_ready_reports_missing_distribution_release_record(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "docs" / "distribution-release.json").unlink()
    (tmp_path / "scripts" / "check_distribution_record.py").unlink()

    result = checker.audit(tmp_path)

    assert any("docs/distribution-release.json" in failure for failure in result.failures)
    assert any("scripts/check_distribution_record.py" in failure for failure in result.failures)


def test_repository_ready_reports_missing_bioconda_recipe_verifier(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "scripts" / "check_bioconda_recipe.py").unlink()

    result = checker.audit(tmp_path)

    assert any("scripts/check_bioconda_recipe.py" in failure for failure in result.failures)


def test_repository_ready_reports_missing_native_comparator_scope(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "docs" / "native-comparator-scope.md").unlink()
    (tmp_path / "scripts" / "check_native_comparator_scope.py").unlink()

    result = checker.audit(tmp_path)

    assert any("docs/native-comparator-scope.md" in failure for failure in result.failures)
    assert any("scripts/check_native_comparator_scope.py" in failure for failure in result.failures)


def test_repository_ready_requires_native_comparator_boundary(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "docs" / "native-comparator-scope.md").write_text(
        "# Native Comparator Scope\n\nSeqAn and Parasail are fast.\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("native comparator scope" in failure for failure in result.failures)


def test_repository_ready_reports_missing_workflow_adoption_artifacts(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "docs" / "workflow-adoption.json").unlink()
    (tmp_path / "scripts" / "check_workflow_adoption.py").unlink()

    result = checker.audit(tmp_path)

    assert any("docs/workflow-adoption.json" in failure for failure in result.failures)
    assert any("scripts/check_workflow_adoption.py" in failure for failure in result.failures)


def test_repository_ready_rejects_unignored_large_or_macos_files(tmp_path, monkeypatch):
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


def test_repository_ready_requires_evidence_boundaries_status(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "docs" / "scientific-claims.md").write_text("# Evidence Notes\n", encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("barcode-comparison-gate" in failure for failure in result.failures)
    assert any("bcl-comparison-gate" in failure for failure in result.failures)


def test_repository_ready_rejects_mismatched_release_versions(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    (tmp_path / "CITATION.cff").write_text(
        "cff-version: 1.2.0\ntitle: DotMatch\nversion: \"0.1.0-dev\"\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("version mismatch" in failure for failure in result.failures)


def test_repository_ready_rejects_dev_readme_status(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    readme = (tmp_path / "README.md").read_text(encoding="utf-8")
    (tmp_path / "README.md").write_text(readme.replace("v0.1.0", "v0.1.0-dev"), encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("README.md must not describe the release version as dev" in failure for failure in result.failures)


def test_repository_ready_rejects_local_absolute_paths(tmp_path):
    checker = _load_checker()
    _write_minimal_repo(tmp_path)
    raw = tmp_path / "benchmarks" / "raw" / "example.csv"
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_text("tool,command\nx,/" + "Users/alice/projects/dotmatch/dotmatch count\n", encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("local absolute path" in failure and "example.csv" in failure for failure in result.failures)
