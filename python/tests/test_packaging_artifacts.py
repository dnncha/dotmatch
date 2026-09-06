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


def _squash_ws(text: str) -> str:
    return " ".join(text.split())


def test_bioconda_recipe_tracks_release_metadata() -> None:
    recipe = ROOT / "packaging" / "bioconda" / "meta.yaml"
    text = recipe.read_text(encoding="utf-8")

    assert '{% set name = "dotmatch" %}' in text
    assert f'{{% set version = "{_pyproject_version()}" %}}' in text
    assert "https://github.com/dnncha/dotmatch/archive/refs/tags/v{{ version }}.tar.gz" in text
    assert "license: Apache-2.0" in text
    assert "license_file: LICENSE" in text
    assert "summary: Known-target short-DNA assignment from FASTQ" in text
    assert "recipe-maintainers:" in text


def test_bioconda_recipe_builds_python_console_script_and_smoke_tests() -> None:
    recipe = (ROOT / "packaging" / "bioconda" / "meta.yaml").read_text(encoding="utf-8")
    build = (ROOT / "packaging" / "bioconda" / "build.sh").read_text(encoding="utf-8")

    assert "- {{ compiler('c') }}" in recipe
    assert "- {{ stdlib('c') }}" in recipe
    assert "{{ pin_subpackage(\"dotmatch\", max_pin=\"x.x\") }}" in recipe
    assert "- make" in recipe
    assert "- python >=3.9" in recipe
    assert "- pip" in recipe
    assert "- setuptools >=77" in recipe
    assert recipe.count("- zlib") == 1
    assert "dotmatch dist ACGT AGGT | grep '^1$'" in recipe
    assert "dotmatch leq 1 ACGT AGGT | grep '^true$'" in recipe
    assert "dotmatch assay --help" in recipe
    assert "dotmatch barcode --help" in recipe
    assert "dotmatch feature --help" in recipe
    assert "dotmatch feature matrix" in recipe
    assert "dotmatch panel --help" in recipe
    assert "dotmatch assay init" in recipe
    assert "dotmatch barcode infer" in recipe
    assert "dotmatch panel design" in recipe
    assert "dotmatch count --help | grep 'Hamming supports k=0..3'" in recipe
    assert "dotmatch crispr-count --help | grep 'MAGeCK-ready'" in recipe
    assert "dotmatch audit --help | grep 'safe_at_hamming_k3'" in recipe
    assert "dotmatch audit --targets audit_targets.tsv --k 3 --audit-mode exact" in recipe
    assert "dotmatch crispr-count --library crispr_guides.csv" in recipe
    assert 'CC="${CC}"' in build
    assert "libdotmatch.a shared" in build
    assert "${PYTHON} -m pip install . -vv --no-deps --no-build-isolation" in build
    assert 'install -m 755 dotmatch "${PREFIX}/bin/dotmatch"' not in build
    assert 'install -m 644 include/qdalign.h "${PREFIX}/include/qdalign.h"' in build


def test_bioconda_recipe_gate_is_wired_into_release_ready() -> None:
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")

    assert "bioconda-recipe-ready:" in makefile
    assert re.search(r"^release-ready: .*bioconda-recipe-ready", makefile, flags=re.MULTILINE)
    assert "python3 scripts/check_bioconda_recipe.py" in makefile


def test_zenodo_metadata_tracks_the_concept_doi_before_release_minting() -> None:
    metadata = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))

    assert metadata["title"] == "DotMatch: deterministic known-target short-DNA assignment for sequencing workflows"
    assert metadata["upload_type"] == "software"
    assert metadata["version"] == _pyproject_version()
    assert metadata["license"] == "Apache-2.0"
    assert metadata["access_right"] == "open"
    assert metadata["creators"] == [
        {
            "name": "O'Toole, Donncha",
            "orcid": "0009-0003-5012-7229",
            "affiliation": "Independent researcher",
        }
    ]
    assert "known-target assignment" in metadata["keywords"]
    assert "doi" not in metadata
    assert metadata["conceptdoi"] == "10.5281/zenodo.20541628"


def test_codemeta_tracks_release_citation_and_concept_doi_after_minting() -> None:
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
    assert codemeta["identifier"] == "https://doi.org/10.5281/zenodo.20541628"
    assert codemeta["author"] == [
        {
            "@type": "Person",
            "givenName": "Donncha",
            "familyName": "O'Toole",
            "@id": "https://orcid.org/0009-0003-5012-7229",
        }
    ]
    assert f"version: \"{_pyproject_version()}\"" in citation
    record = json.loads((ROOT / "docs/distribution-release.json").read_text(encoding="utf-8"))
    channel = next(item for item in record["channels"] if item["id"] == "zenodo")
    doi_lines = [line for line in citation.splitlines() if line.startswith("doi:")]
    if doi_lines:
        assert record["release_version"] == _pyproject_version()
        assert channel["status"] == "verified"
        assert doi_lines[0].split(":", 1)[1].strip().strip('"') == channel["version_doi"]
    else:
        assert channel["status"] in {"prepared", "blocked"}
        assert "doi" not in zenodo
    assert 'doi: "10.5281/zenodo.22214073"' not in citation
    assert codemeta["softwareVersion"] == zenodo["version"]
    assert "known-target assignment" in codemeta["keywords"]
    assert "CRISPR" in codemeta["keywords"]


def test_codemeta_is_included_in_source_distribution_manifest() -> None:
    manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    verifier = (ROOT / "scripts" / "check_python_wheel.py").read_text(encoding="utf-8")

    assert "include CITATION.cff" in manifest
    assert "include codemeta.json" in manifest
    assert "include docs/assay-evidence.json" in manifest
    assert "include include/qdmetal.h" in manifest
    assert "include src/qdmetal_stub.c" in manifest
    assert "/CITATION.cff" in verifier
    assert "/codemeta.json" in verifier
    assert "/docs/assay-evidence.json" in verifier
    assert "dotmatch/data/assay-evidence.json" in verifier
    assert "evidence_boundary" in verifier


def test_python_package_verifier_checks_installed_cli_version() -> None:
    verifier = (ROOT / "scripts" / "check_python_wheel.py").read_text(encoding="utf-8")

    assert "project_version()" in verifier
    assert '"dotmatch.cli", "--version"' in verifier
    assert 'venv_script(env_dir, "dotmatch")' in verifier
    assert '"--version"' in verifier
    assert "dotmatch-native" in verifier
    assert '"assay", "check"' in verifier
    assert '"crispr"' in verifier
    assert '"qc"' in verifier
    assert '"crispr-qc"' in verifier
    assert '"infer"' in verifier
    assert '"autopsy"' in verifier


def test_python_package_verifier_smokes_feature_matrix_and_paired_fastq() -> None:
    verifier = (ROOT / "scripts" / "check_python_wheel.py").read_text(encoding="utf-8")

    assert '"feature",' in verifier
    assert '"matrix",' in verifier
    assert "feature_matrix" in verifier
    assert "matrix.mtx" in verifier
    assert '"pair-count",' in verifier
    assert '"--left-reads",' in verifier
    assert '"--right-reads",' in verifier
    assert "paired FASTQ pair-count" in verifier


def test_python_package_build_bundles_native_cli() -> None:
    setup = (ROOT / "setup.py").read_text(encoding="utf-8")
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")

    assert "dotmatch-native" in setup
    assert "assay-evidence.json" in setup
    assert "src/qda.c" in setup
    assert "src/qdmetal_stub.c" in setup
    assert "DOTMATCH_VERSION" in setup
    assert 'tomli; python_version < \\"3.11\\"' in pyproject


def test_release_workflow_builds_and_smoke_tests_container() -> None:
    workflow = (ROOT / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")

    assert "packages: write" in workflow
    assert "docker/setup-buildx-action" in workflow
    assert "docker/metadata-action" in workflow
    assert "docker/build-push-action" in workflow
    assert "ghcr.io/dnncha/dotmatch" in workflow
    assert "VERSION=$(python -c" in workflow
    assert 'docker run --rm dotmatch:ci --version | grep "^dotmatch ${VERSION}$"' in workflow
    assert "docker run --rm dotmatch:ci dist ACGT AGGT | grep '^1$'" in workflow
    assert "docker image inspect dotmatch:ci" in workflow
    assert "org.opencontainers.image.version" in workflow
    assert "org.opencontainers.image.title=DotMatch" in workflow
    assert "org.opencontainers.image.licenses=Apache-2.0" in workflow
    assert "org.opencontainers.image.documentation=https://dotmatch.readthedocs.io/" in workflow
    assert "org.opencontainers.image.authors=Donncha O'Toole" in workflow


def test_dockerfile_has_release_aligned_oci_metadata() -> None:
    dockerfile = (ROOT / "Dockerfile").read_text(encoding="utf-8")

    assert 'LABEL org.opencontainers.image.title="DotMatch"' in dockerfile
    assert 'org.opencontainers.image.source="https://github.com/dnncha/dotmatch"' in dockerfile
    assert 'org.opencontainers.image.url="https://dotmatch.readthedocs.io/"' in dockerfile
    assert f'org.opencontainers.image.version="{_pyproject_version()}"' in dockerfile
    assert 'org.opencontainers.image.licenses="Apache-2.0"' in dockerfile
    assert 'org.opencontainers.image.description=' in dockerfile
    assert 'org.opencontainers.image.documentation="https://dotmatch.readthedocs.io/"' in dockerfile
    assert 'org.opencontainers.image.authors=' in dockerfile


def test_release_workflow_publishes_pypi_sdist_and_repaired_linux_wheels() -> None:
    workflow = (ROOT / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")
    packaging = (ROOT / "docs" / "packaging.md").read_text(encoding="utf-8")

    assert "id-token: write" in workflow
    assert "python scripts/check_python_wheel.py --sdist-only --out-dir dist" in workflow
    assert "Publish PyPI sdist, macOS wheel, and repaired Linux wheels" in workflow
    assert "pypa/gh-action-pypi-publish@release/v1" in workflow
    assert "name: dotmatch-sdist" in workflow
    assert "name: dotmatch-linux-repaired-wheels" in workflow
    assert "needs: [preflight, sdist, wheel, linux-repaired-wheels]" in workflow
    assert "path: dist-pypi" in workflow
    assert "packages-dir: dist-pypi" in workflow
    assert "dotmatch-wheel-Linux" not in workflow
    assert "trusted publishing" in packaging
    assert "source distribution plus a macOS wheel and" in packaging
    assert "repaired manylinux/musllinux wheels" in packaging


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
    assert "repaired manylinux/musllinux wheels" in packaging
    assert "https://pypi.org/project/dotmatch/" in readme
    assert "packaging.html" in readme


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
    assert "needs: [preflight, sdist, wheel, linux-repaired-wheels]" in workflow
    assert "needs: [preflight, wheel, sdist, linux-repaired-wheels]" in workflow


def test_distribution_docs_include_clean_pypi_install_verification() -> None:
    packaging = (ROOT / "docs" / "packaging.md").read_text(encoding="utf-8")
    checker = (ROOT / "scripts" / "check_distribution_channels.py").read_text(encoding="utf-8")

    assert "pip install dotmatch==" in checker
    assert "must include repaired manylinux and musllinux wheels" in checker
    assert "must not include raw" in checker
    assert "source distribution plus a macOS wheel" in packaging
    assert "repaired manylinux/musllinux wheels" in packaging
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
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    checker = (ROOT / "scripts" / "check_distribution_channels.py").read_text(encoding="utf-8")
    version = _pyproject_version()

    assert "quay.io/api/v1/repository/biocontainers/dotmatch/tag/" in checker
    assert '"docker", "run", "--rm", image, "dotmatch", "leq", "1", "ACGT", "AGGT"' in checker
    assert "quay.io/biocontainers/dotmatch:<version>--<build>" in packaging
    assert "BioContainers images for DotMatch are generated from the accepted Bioconda" in packaging
    assert "newly tagged version has not reached Bioconda yet" in readme
