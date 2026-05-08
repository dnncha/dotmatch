from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_distribution_submission.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_distribution_submission", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_repo(root: Path, dossier: str | None = None) -> None:
    files = {
        "docs/distribution-submission.md": dossier
        or (
            "# Distribution Submission Dossier\n\n"
            "This is a pre-release handoff, not a public distribution claim.\n\n"
            "## Required Checks\n\n"
            "```bash\n"
            "make release-ready\n"
            "make python-package-test\n"
            "make bioconda-recipe-ready\n"
            "make distribution-record-ready\n"
            "make distribution-submission-ready\n"
            "make distribution-channels\n"
            "```\n\n"
            "## Submission Targets\n\n"
            "- PyPI trusted publishing: push the immutable v0.1.0 tag and publish the source distribution "
            "verified by scripts/check_python_wheel.py --sdist-only --out-dir dist plus repaired manylinux/musllinux wheels "
            "from the dotmatch-linux-repaired-wheels artifact.\n"
            "- GHCR: verify ghcr.io/dnncha/dotmatch:v0.1.0, OCI labels, `dotmatch --version`, and `dotmatch dist ACGT AGGT`.\n"
            "- Bioconda: replace REPLACE_WITH_RELEASE_TARBALL_SHA256 in packaging/bioconda/meta.yaml and submit to bioconda-recipes.\n"
            "- BioContainers: verify quay.io/biocontainers/dotmatch after Bioconda acceptance.\n"
            "- Zenodo: archive the immutable GitHub release, mint a DOI, and add it to CITATION.cff.\n\n"
            "## Record Update\n\n"
            "Update docs/distribution-release.json with `public_url`, `evidence_url`, and `verified_date` for every channel.\n"
            "Do not set status to `released` until stable public HTTPS URLs exist and `make distribution-channels` passes.\n"
        ),
    }
    for path, text in files.items():
        full = root / path
        full.parent.mkdir(parents=True, exist_ok=True)
        full.write_text(text, encoding="utf-8")


def test_distribution_submission_accepts_complete_dossier(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)

    result = checker.audit(tmp_path)

    assert result.failures == []
    assert any("distribution submission dossier present" in item for item in result.passed)


def test_distribution_submission_rejects_missing_dossier(tmp_path):
    checker = _load_checker()

    result = checker.audit(tmp_path)

    assert any("docs/distribution-submission.md could not be read" in failure for failure in result.failures)


def test_distribution_submission_rejects_premature_public_claim(tmp_path):
    checker = _load_checker()
    _write_repo(
        tmp_path,
        dossier="# Distribution Submission Dossier\n\nDotMatch is now available on PyPI and Bioconda.\n",
    )

    result = checker.audit(tmp_path)

    assert any("must avoid public distribution claims" in failure for failure in result.failures)


def test_distribution_submission_requires_all_release_channels(tmp_path):
    checker = _load_checker()
    _write_repo(
        tmp_path,
        dossier=(
            "# Distribution Submission Dossier\n\n"
            "This is a pre-release handoff, not a public distribution claim.\n\n"
            "make release-ready\n"
            "make python-package-test\n"
            "make bioconda-recipe-ready\n"
            "make distribution-record-ready\n"
            "make distribution-submission-ready\n"
            "make distribution-channels\n"
            "PyPI trusted publishing only.\n"
            "Do not set status to `released` until stable public HTTPS URLs exist and `make distribution-channels` passes.\n"
        ),
    )

    result = checker.audit(tmp_path)

    assert any("must include GHCR" in failure for failure in result.failures)
    assert any("must include Bioconda" in failure for failure in result.failures)
    assert any("must include Zenodo" in failure for failure in result.failures)


def test_distribution_submission_requires_manifest_update_fields(tmp_path):
    checker = _load_checker()
    _write_repo(
        tmp_path,
        dossier=(
            "# Distribution Submission Dossier\n\n"
            "This is a pre-release handoff, not a public distribution claim.\n\n"
            "make release-ready\n"
            "make python-package-test\n"
            "make bioconda-recipe-ready\n"
            "make distribution-record-ready\n"
            "make distribution-submission-ready\n"
            "make distribution-channels\n"
            "PyPI trusted publishing, source distribution, repaired manylinux/musllinux wheels, "
            "dotmatch-linux-repaired-wheels, scripts/check_python_wheel.py --sdist-only --out-dir dist, "
            "GHCR, Bioconda, BioContainers, Zenodo, CITATION.cff, "
            "packaging/bioconda/meta.yaml, REPLACE_WITH_RELEASE_TARBALL_SHA256, bioconda-recipes, "
            "ghcr.io/dnncha/dotmatch, quay.io/biocontainers/dotmatch.\n"
            "Do not set status to `released` yet.\n"
        ),
    )

    result = checker.audit(tmp_path)

    assert any("must include public_url" in failure for failure in result.failures)
    assert any("must include evidence_url" in failure for failure in result.failures)
    assert any("must include verified_date" in failure for failure in result.failures)


def test_distribution_submission_requires_bioconda_recipe_gate(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    dossier = (tmp_path / "docs" / "distribution-submission.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "distribution-submission.md").write_text(
        dossier.replace("make bioconda-recipe-ready\n", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("make bioconda-recipe-ready" in failure for failure in result.failures)


def test_distribution_submission_requires_python_package_gate(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    dossier = (tmp_path / "docs" / "distribution-submission.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "distribution-submission.md").write_text(
        dossier.replace("make python-package-test\n", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("make python-package-test" in failure for failure in result.failures)


def test_distribution_submission_requires_verified_pypi_sdist_handoff(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    dossier = (tmp_path / "docs" / "distribution-submission.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "distribution-submission.md").write_text(
        dossier.replace(" verified by scripts/check_python_wheel.py --sdist-only --out-dir dist", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("PyPI: scripts/check_python_wheel.py --sdist-only --out-dir dist" in failure for failure in result.failures)


def test_distribution_submission_requires_repaired_pypi_wheel_handoff(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    dossier = (tmp_path / "docs" / "distribution-submission.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "distribution-submission.md").write_text(
        dossier.replace(" plus repaired manylinux/musllinux wheels from the dotmatch-linux-repaired-wheels artifact", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("PyPI: repaired manylinux/musllinux wheels" in failure for failure in result.failures)
    assert any("PyPI: dotmatch-linux-repaired-wheels" in failure for failure in result.failures)
