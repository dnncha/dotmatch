from __future__ import annotations

import re
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def _pyproject_version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    assert match is not None
    return match.group(1)


def test_bioconda_recipe_tracks_release_metadata() -> None:
    recipe = ROOT / "packaging" / "bioconda" / "meta.yaml"
    text = recipe.read_text(encoding="utf-8")

    assert '{% set name = "dotmatch" %}' in text
    assert f'{{% set version = "{_pyproject_version()}" %}}' in text
    assert "https://github.com/dnncha/dotmatch/archive/refs/tags/v{{ version }}.tar.gz" in text
    assert "license: Apache-2.0" in text
    assert "license_file: LICENSE" in text
    assert "summary: Fast exact short-DNA known-target assignment" in text
    assert "recipe-maintainers:" in text


def test_bioconda_recipe_builds_native_cli_and_smoke_tests() -> None:
    recipe = (ROOT / "packaging" / "bioconda" / "meta.yaml").read_text(encoding="utf-8")
    build = (ROOT / "packaging" / "bioconda" / "build.sh").read_text(encoding="utf-8")

    assert "- {{ compiler('c') }}" in recipe
    assert "- {{ stdlib('c') }}" in recipe
    assert "{{ pin_subpackage(\"dotmatch\", max_pin=\"x.x\") }}" in recipe
    assert "- make" in recipe
    assert "- zlib" in recipe
    assert "dotmatch dist ACGT AGGT | grep '^1$'" in recipe
    assert "dotmatch leq 1 ACGT AGGT | grep '^true$'" in recipe
    assert "dotmatch count --help" not in recipe
    assert 'CC="${CC}"' in build
    assert "dotmatch libdotmatch.a shared" in build
    assert 'install -m 755 dotmatch "${PREFIX}/bin/dotmatch"' in build
    assert 'install -m 644 include/qdalign.h "${PREFIX}/include/qdalign.h"' in build


def test_bioconda_recipe_gate_is_wired_into_release_ready() -> None:
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")

    assert "bioconda-recipe-ready:" in makefile
    assert re.search(r"^release-ready: .*bioconda-recipe-ready", makefile, flags=re.MULTILINE)
    assert "python3 scripts/check_bioconda_recipe.py" in makefile


def test_zenodo_metadata_is_release_ready_without_claiming_doi() -> None:
    metadata = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))

    assert metadata["title"] == "DotMatch: Streaming Exact One-Edit Barcode and Guide Assignment Without Exhaustive Scanning"
    assert metadata["upload_type"] == "software"
    assert metadata["version"] == _pyproject_version()
    assert metadata["license"] == "Apache-2.0"
    assert metadata["access_right"] == "open"
    assert metadata["creators"] == [{"name": "O'Toole, Donncha"}]
    assert "known-target assignment" in metadata["keywords"]
    assert all("doi.org" not in str(value).lower() for value in metadata.values())


def test_codemeta_tracks_package_citation_and_no_doi_claim() -> None:
    codemeta = json.loads((ROOT / "codemeta.json").read_text(encoding="utf-8"))
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    zenodo = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))

    assert codemeta["@type"] == "SoftwareSourceCode"
    assert codemeta["name"] == "DotMatch"
    assert codemeta["codeRepository"] == "https://github.com/dnncha/dotmatch"
    assert codemeta["url"] == "https://github.com/dnncha/dotmatch"
    assert codemeta["version"] == _pyproject_version()
    assert codemeta["softwareVersion"] == _pyproject_version()
    assert codemeta["license"] == "https://spdx.org/licenses/Apache-2.0"
    assert codemeta["citation"].endswith("/CITATION.cff")
    assert codemeta["author"] == [{"@type": "Person", "givenName": "Donncha", "familyName": "O'Toole"}]
    assert "version: \"0.1.0\"" in citation
    assert codemeta["softwareVersion"] == zenodo["version"]
    assert "known-target assignment" in codemeta["keywords"]
    assert "CRISPR" in codemeta["keywords"]
    doi_claim_fields = {key: value for key, value in codemeta.items() if key != "@context"}
    assert all("doi.org" not in json.dumps(value).lower() for value in doi_claim_fields.values())


def test_codemeta_is_included_in_source_distribution_manifest() -> None:
    manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    verifier = (ROOT / "scripts" / "check_python_wheel.py").read_text(encoding="utf-8")

    assert "include CITATION.cff" in manifest
    assert "include codemeta.json" in manifest
    assert "/CITATION.cff" in verifier
    assert "/codemeta.json" in verifier


def test_python_package_verifier_checks_installed_cli_version() -> None:
    verifier = (ROOT / "scripts" / "check_python_wheel.py").read_text(encoding="utf-8")

    assert "project_version()" in verifier
    assert '"dotmatch.cli", "--version"' in verifier
    assert 'venv_script(env_dir, "dotmatch")' in verifier
    assert '"--version"' in verifier
    assert "dotmatch-native" in verifier
    assert '"assay", "check"' in verifier
    assert '"infer"' in verifier
    assert '"autopsy"' in verifier


def test_python_package_build_bundles_native_cli() -> None:
    setup = (ROOT / "setup.py").read_text(encoding="utf-8")
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")

    assert "dotmatch-native" in setup
    assert "src/qda.c" in setup
    assert "DOTMATCH_VERSION" in setup
    assert 'tomli; python_version < \\"3.11\\"' in pyproject


def test_release_workflow_builds_and_smoke_tests_container() -> None:
    workflow = (ROOT / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")

    assert "packages: write" in workflow
    assert "docker/setup-buildx-action" in workflow
    assert "docker/metadata-action" in workflow
    assert "docker/build-push-action" in workflow
    assert "ghcr.io/dnncha/dotmatch" in workflow
    assert f"docker run --rm dotmatch:ci --version | grep '^dotmatch {_pyproject_version()}$'" in workflow
    assert "docker run --rm dotmatch:ci dist ACGT AGGT | grep '^1$'" in workflow
    assert "docker image inspect dotmatch:ci" in workflow
    assert "org.opencontainers.image.version" in workflow


def test_dockerfile_has_release_aligned_oci_metadata() -> None:
    dockerfile = (ROOT / "Dockerfile").read_text(encoding="utf-8")

    assert 'LABEL org.opencontainers.image.title="DotMatch"' in dockerfile
    assert 'org.opencontainers.image.source="https://github.com/dnncha/dotmatch"' in dockerfile
    assert 'org.opencontainers.image.url="https://github.com/dnncha/dotmatch"' in dockerfile
    assert f'org.opencontainers.image.version="{_pyproject_version()}"' in dockerfile
    assert 'org.opencontainers.image.licenses="Apache-2.0"' in dockerfile
    assert 'org.opencontainers.image.description=' in dockerfile
    assert 'org.opencontainers.image.documentation="https://github.com/dnncha/dotmatch#readme"' in dockerfile
    assert 'org.opencontainers.image.authors=' in dockerfile


def test_release_workflow_publishes_pypi_sdist_and_repaired_linux_wheels() -> None:
    workflow = (ROOT / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")
    packaging = (ROOT / "docs" / "packaging.md").read_text(encoding="utf-8")

    assert "id-token: write" in workflow
    assert "python scripts/check_python_wheel.py --sdist-only --out-dir dist" in workflow
    assert "Publish PyPI sdist and repaired Linux wheels" in workflow
    assert "pypa/gh-action-pypi-publish@release/v1" in workflow
    assert "name: dotmatch-sdist" in workflow
    assert "name: dotmatch-linux-repaired-wheels" in workflow
    assert "needs: [preflight, sdist, linux-repaired-wheels]" in workflow
    assert "path: dist-pypi" in workflow
    assert "packages-dir: dist-pypi" in workflow
    assert "dotmatch-wheel-Linux" not in workflow
    assert "trusted publishing" in packaging
    assert "publishes the source distribution and repaired manylinux/musllinux wheels" in packaging


def test_cibuildwheel_linux_repaired_wheel_path_is_configured() -> None:
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    workflow = (ROOT / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")
    packaging = (ROOT / "docs" / "packaging.md").read_text(encoding="utf-8")
    readme = (ROOT / "README.md").read_text(encoding="utf-8")

    assert "[tool.cibuildwheel]" in pyproject
    assert "cp39-manylinux_x86_64" in pyproject
    assert "cp312-musllinux_x86_64" in pyproject
    assert "dotmatch dist ACGT AGGT" in pyproject
    assert "pypa/cibuildwheel" in workflow
    assert "dotmatch-linux-repaired-wheels" in workflow
    assert "dist-linux/*.whl" in workflow
    assert "manylinux/musllinux" in packaging
    assert "repaired Linux wheel artifacts" in packaging
    assert "[x] GitHub release manylinux/musllinux x86_64 wheel artifact build" in readme
    assert "[x] PyPI trusted-publishing path for repaired manylinux/musllinux Linux wheels" in readme
    assert "[ ] Public PyPI availability of repaired manylinux/musllinux Linux wheels" in readme


def test_release_workflow_publishing_jobs_depend_on_preflight() -> None:
    workflow = (ROOT / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")

    assert "\n  preflight:" in workflow
    assert "Release preflight gates" in workflow
    assert "python -m pip install build pytest" in workflow
    assert "make test" in workflow
    assert "make cli-test" in workflow
    assert "make python-test" in workflow
    assert "make repository-ready" in workflow
    assert "make release-ready" in workflow
    assert "make python-package-test" in workflow
    assert "needs: [preflight]" in workflow
    assert "needs: [preflight, sdist, linux-repaired-wheels]" in workflow
    assert "needs: [preflight, wheel, sdist, linux-repaired-wheels]" in workflow


def test_distribution_docs_include_clean_pypi_install_verification() -> None:
    packaging = (ROOT / "docs" / "packaging.md").read_text(encoding="utf-8")
    checker = (ROOT / "scripts" / "check_distribution_channels.py").read_text(encoding="utf-8")

    assert "pip install dotmatch==" in checker
    assert "must include repaired manylinux and musllinux wheels" in checker
    assert "must not include raw linux_x86_64 wheels" in checker
    assert "source distribution plus repaired manylinux/musllinux wheels" in packaging
    assert "rejects raw" in packaging and "linux_x86_64" in packaging
    assert "clean virtual environment" in packaging
    assert "pip install dotmatch==<version>" in packaging


def test_distribution_docs_include_ghcr_runtime_verification() -> None:
    packaging = (ROOT / "docs" / "packaging.md").read_text(encoding="utf-8")
    checker = (ROOT / "scripts" / "check_distribution_channels.py").read_text(encoding="utf-8")

    assert '"docker", "run", "--rm", image, "--version"' in checker
    assert "docker run --rm ghcr.io/dnncha/dotmatch:v<version>" in packaging


def test_distribution_docs_include_bioconda_install_verification() -> None:
    packaging = (ROOT / "docs" / "packaging.md").read_text(encoding="utf-8")
    checker = (ROOT / "scripts" / "check_distribution_channels.py").read_text(encoding="utf-8")

    assert '"micromamba"' in checker
    assert "conda create -p <env> -c conda-forge -c bioconda dotmatch=<version>" in packaging


def test_distribution_docs_include_biocontainers_runtime_verification() -> None:
    packaging = (ROOT / "docs" / "packaging.md").read_text(encoding="utf-8")
    checker = (ROOT / "scripts" / "check_distribution_channels.py").read_text(encoding="utf-8")

    assert "quay.io/api/v1/repository/biocontainers/dotmatch/tag/" in checker
    assert '"docker", "run", "--rm", image, "dotmatch", "leq", "1", "ACGT", "AGGT"' in checker
    assert "quay.io/biocontainers/dotmatch:<version>--<build>" in packaging
