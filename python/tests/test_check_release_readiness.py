import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_release_readiness.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_release_readiness", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_release_repo(root: Path) -> None:
    files = {
        "pyproject.toml": '[project]\nname = "dotmatch"\nversion = "0.1.0"\nlicense = "Apache-2.0"\n',
        "package.json": '{"version": "0.1.0", "license": "Apache-2.0"}\n',
        "CITATION.cff": (
            'cff-version: 1.2.0\n'
            'title: "DotMatch"\n'
            'version: "0.1.0"\n'
            'repository-code: "https://github.com/dnncha/dotmatch"\n'
            "license: Apache-2.0\n"
        ),
        "codemeta.json": (
            '{"name": "DotMatch", "version": "0.1.0", "softwareVersion": "0.1.0", '
            '"license": "https://spdx.org/licenses/Apache-2.0", '
            '"codeRepository": "https://github.com/dnncha/dotmatch", '
            '"keywords": ["bioinformatics", "known-target assignment"]}\n'
        ),
        ".zenodo.json": (
            '{"title": "DotMatch", "upload_type": "software", "version": "0.1.0", '
            '"license": "Apache-2.0", "access_right": "open", '
            '"keywords": ["known-target assignment"]}\n'
        ),
        "MANIFEST.in": "include CITATION.cff\ninclude codemeta.json\ninclude src/qdalign.c\ninclude include/qdalign.h\n",
        "scripts/check_python_wheel.py": (
            'required_suffixes = ["/CITATION.cff", "/codemeta.json", '
            '"/src/qdalign.c", "/include/qdalign.h"]\n'
        ),
        "Dockerfile": (
            'FROM debian:bookworm-slim\n'
            'LABEL org.opencontainers.image.title="DotMatch" \\\n'
            '      org.opencontainers.image.source="https://github.com/dnncha/dotmatch" \\\n'
            '      org.opencontainers.image.url="https://github.com/dnncha/dotmatch" \\\n'
            '      org.opencontainers.image.version="0.1.0" \\\n'
            '      org.opencontainers.image.licenses="Apache-2.0"\n'
        ),
        "packaging/bioconda/meta.yaml": (
            '{% set name = "dotmatch" %}\n'
            '{% set version = "0.1.0" %}\n'
            "source:\n"
            "  url: https://github.com/dnncha/dotmatch/archive/refs/tags/v{{ version }}.tar.gz\n"
            "  sha256: REPLACE_WITH_RELEASE_TARBALL_SHA256\n"
            "test:\n"
            "  commands:\n"
            "    - dotmatch dist ACGT AGGT | grep '^1$'\n"
        ),
        "packaging/bioconda/build.sh": "#!/usr/bin/env bash\nmake dotmatch libdotmatch.a shared\n",
        ".github/workflows/release.yml": (
            "permissions:\n"
            "  id-token: write\n"
            "  packages: write\n"
            "jobs:\n"
            "  preflight:\n"
            "    name: Release preflight gates\n"
            "    steps:\n"
            "      - run: python -m pip install build pytest\n"
            "      - run: make test\n"
            "      - run: make cli-test\n"
            "      - run: make python-test\n"
            "      - run: make repository-ready\n"
            "      - run: make release-ready\n"
            "      - run: make python-package-test\n"
            "  container:\n"
            "    needs: [preflight]\n"
            "    steps:\n"
            "      - uses: docker/metadata-action@v5\n"
            "        with:\n"
            "          images: ghcr.io/dnncha/dotmatch\n"
            "      - uses: docker/build-push-action@v6\n"
            "      - run: docker image inspect dotmatch:ci --format '{{ index .Config.Labels \"org.opencontainers.image.version\" }}'\n"
            "  sdist:\n"
            "    steps:\n"
            "      - run: python scripts/check_python_wheel.py --sdist-only --out-dir dist\n"
            "      - uses: actions/upload-artifact@v7\n"
            "        with:\n"
            "          name: dotmatch-sdist\n"
            "          path: dist/*.tar.gz\n"
            "  wheel:\n"
            "    needs: [preflight]\n"
            "    steps:\n"
            "      - uses: actions/upload-artifact@v7\n"
            "        with:\n"
            "          name: dotmatch-wheel-macos\n"
            "          path: dist/*.whl\n"
            "  linux-repaired-wheels:\n"
            "    needs: [preflight]\n"
            "    steps:\n"
            "      - uses: pypa/cibuildwheel@v3.3.0\n"
            "      - uses: actions/upload-artifact@v7\n"
            "        with:\n"
            "          name: dotmatch-linux-repaired-wheels\n"
            "          path: dist-linux/*.whl\n"
            "  pypi-sdist:\n"
            "    name: Publish PyPI sdist, macOS wheel, and repaired Linux wheels\n"
            "    needs: [preflight, sdist, wheel, linux-repaired-wheels]\n"
            "    steps:\n"
            "      - uses: actions/download-artifact@v8\n"
            "        with:\n"
            "          name: dotmatch-sdist\n"
            "          path: dist-pypi\n"
            "      - uses: actions/download-artifact@v8\n"
            "        with:\n"
            "          name: dotmatch-wheel-macos\n"
            "          path: dist-pypi\n"
            "      - uses: actions/download-artifact@v8\n"
            "        with:\n"
            "          name: dotmatch-linux-repaired-wheels\n"
            "          path: dist-pypi\n"
            "      - uses: pypa/gh-action-pypi-publish@release/v1\n"
            "        with:\n"
            "          packages-dir: dist-pypi\n"
            "  github-release:\n"
            "    needs: [preflight, wheel, sdist, linux-repaired-wheels]\n"
            "    steps:\n"
            "      - run: sha256sum * > SHA256SUMS.txt\n"
        ),
        "docs/packaging.md": (
            "# Packaging\n\n"
            "The release workflow uses trusted publishing and publishes the source distribution, the native macOS wheel, and repaired manylinux/musllinux Linux wheels.\n"
            "Raw `linux_x86_64` wheels are not uploaded to PyPI.\n"
            "Images are pushed to ghcr.io/dnncha/dotmatch with OCI labels.\n"
            "BioContainers images are checked at quay.io/biocontainers/dotmatch after Bioconda publication.\n"
            "Run make bioconda-recipe-ready before copying the recipe to bioconda-recipes.\n"
            "docs/distribution-release.json records the prepared public package channels before publication.\n"
        ),
        "docs/release-process.md": (
            "# Release Process\n\n"
            "Run `make pretag-ready`, `make release-ready`, `make assay-evidence-ready`, "
            "`make distribution-record-ready`, `make alphabet-policy-ready`, "
            "`make bioconda-recipe-ready`, "
            "`make citation-metadata-ready`, `make native-comparator-scope-ready`, "
            "`make public-crispr-evidence-gate`, `make crispr-comparison-gate`, and "
            "`make barcode-comparison-gate`, `make feature-barcode-public-gate`, and "
            "`make perturb-seq-public-gate`, `make amplicon-panel-public-gate`, and "
            "`make bcl-tiny-public-gate`, `make oligo-adapter-public-gate`, and "
            "`make workflow-examples-ready` before tagging.\n"
            "Keep `make distribution-channels`, `make workflow-adoption-status`, and "
            "`make bcl-comparison-gate` separate because they require external evidence.\n"
            "Publish the PyPI source distribution and repaired manylinux/musllinux wheels through trusted publishing.\n"
        ),
        "Makefile": (
            "pretag-ready:\n"
            "\t$(MAKE) test\n"
            "\t$(MAKE) cli-test\n"
            "\t$(MAKE) python-test\n"
            "\t$(MAKE) python-package-test\n"
            "\t$(MAKE) repository-ready\n"
            "\t$(MAKE) release-ready\n"
            "\t$(MAKE) coverage\n"
            "\tnpm run lint\n"
            "\tnpm audit --audit-level=moderate\n"
            "\tnpm run build\n"
            "\tNEXT_OUTPUT=export NEXT_PUBLIC_BASE_PATH=/dotmatch NEXT_PUBLIC_SITE_URL=https://dnncha.github.io/dotmatch npm run build\n"
        ),
    }
    for path, text in files.items():
        full = root / path
        full.parent.mkdir(parents=True, exist_ok=True)
        full.write_text(text, encoding="utf-8")


def test_release_readiness_accepts_minimal_release_repo(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)

    result = checker.audit(tmp_path)

    assert result.failures == []
    assert any("versions aligned" in item for item in result.passed)
    assert any("distribution surfaces" in item for item in result.passed)


def test_release_readiness_rejects_doi_claim_before_minted_release(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    (tmp_path / "CITATION.cff").write_text(
        'cff-version: 1.2.0\ntitle: "DotMatch"\nversion: "0.1.0"\ndoi: 10.5281/zenodo.123\n',
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("DOI" in failure for failure in result.failures)


def test_release_readiness_rejects_unaligned_container_label(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    dockerfile = (tmp_path / "Dockerfile").read_text(encoding="utf-8").replace(
        'org.opencontainers.image.version="0.1.0"',
        'org.opencontainers.image.version="0.2.0"',
    )
    (tmp_path / "Dockerfile").write_text(dockerfile, encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("Dockerfile" in failure and "version" in failure for failure in result.failures)


def test_release_readiness_rejects_missing_sdist_metadata(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    (tmp_path / "MANIFEST.in").write_text("include src/qdalign.c\n", encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("MANIFEST.in" in failure and "codemeta.json" in failure for failure in result.failures)


def test_release_readiness_requires_assay_evidence_gate_in_release_process(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    (tmp_path / "docs" / "release-process.md").write_text(
        "# Release Process\n\nRun `make release-ready` before tagging.\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("assay-evidence-ready" in failure for failure in result.failures)


def test_release_readiness_requires_pretag_ready_in_release_process(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    release_process = (tmp_path / "docs" / "release-process.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "release-process.md").write_text(
        release_process.replace("`make pretag-ready`, ", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("pretag-ready" in failure for failure in result.failures)


def test_release_readiness_requires_pretag_ready_target(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    (tmp_path / "Makefile").write_text(
        "pretag-ready:\n"
        "\t$(MAKE) test\n"
        "\t$(MAKE) python-test\n"
        "\tnpm run build\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("pretag-ready target must include $(MAKE) cli-test" in failure for failure in result.failures)
    assert any("pretag-ready target must include npm run lint" in failure for failure in result.failures)
    assert any("pretag-ready target must include NEXT_OUTPUT=export" in failure for failure in result.failures)


def test_release_readiness_rejects_post_release_gates_in_pretag_ready(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    makefile = (tmp_path / "Makefile").read_text(encoding="utf-8")
    (tmp_path / "Makefile").write_text(
        makefile.replace("pretag-ready:\n", "pretag-ready:\n\t$(MAKE) distribution-channels\n"),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("pretag-ready target must not include distribution-channels" in failure for failure in result.failures)


def test_release_readiness_requires_post_release_gate_boundaries_in_release_process(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    release_process = (tmp_path / "docs" / "release-process.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "release-process.md").write_text(
        release_process
        .replace("`make distribution-channels`, ", "")
        .replace("`make workflow-adoption-status`, and ", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("separate distribution-channels gate" in failure for failure in result.failures)
    assert any("separate workflow-adoption-status gate" in failure for failure in result.failures)


def test_release_readiness_requires_distribution_record_gate_in_release_process(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    release_process = (tmp_path / "docs" / "release-process.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "release-process.md").write_text(
        release_process.replace("`make distribution-record-ready`, ", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("distribution-record-ready" in failure for failure in result.failures)


def test_release_readiness_requires_repaired_linux_wheel_release_process_note(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    release_process = (tmp_path / "docs" / "release-process.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "release-process.md").write_text(
        release_process.replace("repaired manylinux/musllinux wheels", "source distribution only"),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("repaired manylinux/musllinux PyPI wheel publishing" in failure for failure in result.failures)


def test_release_readiness_requires_bioconda_recipe_gate_in_release_process(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    release_process = (tmp_path / "docs" / "release-process.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "release-process.md").write_text(
        release_process.replace("`make bioconda-recipe-ready`, ", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("bioconda-recipe-ready" in failure for failure in result.failures)


def test_release_readiness_requires_alphabet_policy_gate_in_release_process(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    release_process = (tmp_path / "docs" / "release-process.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "release-process.md").write_text(
        release_process.replace("`make alphabet-policy-ready`, ", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("alphabet-policy-ready" in failure for failure in result.failures)


def test_release_readiness_requires_citation_metadata_gate_in_release_process(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    release_process = (tmp_path / "docs" / "release-process.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "release-process.md").write_text(
        release_process.replace("`make citation-metadata-ready`, ", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("citation-metadata-ready" in failure for failure in result.failures)


def test_release_readiness_requires_native_comparator_scope_gate_in_release_process(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    release_process = (tmp_path / "docs" / "release-process.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "release-process.md").write_text(
        release_process.replace("`make native-comparator-scope-ready`, ", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("native-comparator-scope-ready" in failure for failure in result.failures)


def test_release_readiness_requires_public_evidence_gates_in_release_process(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    (tmp_path / "docs" / "release-process.md").write_text(
        "# Release Process\n\n"
        "Run `make release-ready` and `make assay-evidence-ready` before tagging.\n",
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("public-crispr-evidence-gate" in failure for failure in result.failures)
    assert any("crispr-comparison-gate" in failure for failure in result.failures)
    assert any("barcode-comparison-gate" in failure for failure in result.failures)
    assert any("feature-barcode-public-gate" in failure for failure in result.failures)
    assert any("perturb-seq-public-gate" in failure for failure in result.failures)
    assert any("amplicon-panel-public-gate" in failure for failure in result.failures)
    assert any("bcl-tiny-public-gate" in failure for failure in result.failures)
    assert any("oligo-adapter-public-gate" in failure for failure in result.failures)


def test_release_readiness_requires_workflow_examples_gate_in_release_process(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    release_process = (tmp_path / "docs" / "release-process.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "release-process.md").write_text(
        release_process.replace("`make workflow-examples-ready`", "`make workflows-missing`"),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("workflow-examples-ready" in failure for failure in result.failures)


def test_release_readiness_requires_biocontainers_packaging_docs(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    packaging = (tmp_path / "docs" / "packaging.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "packaging.md").write_text(
        packaging.replace("quay.io/biocontainers/dotmatch", "quay.io/missing/dotmatch"),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("BioContainers" in failure for failure in result.failures)


def test_release_readiness_requires_distribution_record_packaging_docs(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    packaging = (tmp_path / "docs" / "packaging.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "packaging.md").write_text(
        packaging.replace("docs/distribution-release.json", "docs/missing-release.json"),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("distribution-release.json" in failure for failure in result.failures)


def test_release_readiness_requires_verified_sdist_job(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    workflow = (tmp_path / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")
    (tmp_path / ".github" / "workflows" / "release.yml").write_text(
        workflow.replace(
            "      - run: python scripts/check_python_wheel.py --sdist-only --out-dir dist\n",
            "      - run: python -m build --sdist --outdir dist\n",
        ),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("sdist job must verify the PyPI source distribution artifact" in failure for failure in result.failures)


def test_release_readiness_requires_bioconda_recipe_packaging_docs(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    packaging = (tmp_path / "docs" / "packaging.md").read_text(encoding="utf-8")
    (tmp_path / "docs" / "packaging.md").write_text(
        packaging.replace("make bioconda-recipe-ready", "make missing-bioconda-gate"),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("bioconda-recipe-ready" in failure for failure in result.failures)


def test_release_readiness_requires_preflight_before_publish_jobs(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    workflow = (tmp_path / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")
    (tmp_path / ".github" / "workflows" / "release.yml").write_text(
        workflow
        .replace("  preflight:\n", "  missing-preflight:\n")
        .replace("    needs: [preflight]\n", "")
        .replace("    needs: [preflight, sdist, wheel, linux-repaired-wheels]\n", "    needs: [sdist, wheel, linux-repaired-wheels]\n")
        .replace("    needs: [preflight, wheel, sdist, linux-repaired-wheels]\n", "    needs: [wheel, sdist, linux-repaired-wheels]\n"),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("release workflow missing preflight job" in failure for failure in result.failures)
    assert any("container publish job must depend on preflight" in failure for failure in result.failures)
    assert any("PyPI publish job must depend on preflight, sdist, macOS wheel, and repaired Linux wheels" in failure for failure in result.failures)
    assert any("GitHub release job must depend on preflight, wheels, sdist, and repaired Linux wheels" in failure for failure in result.failures)


def test_release_readiness_requires_repository_ready_in_preflight(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    workflow = (tmp_path / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")
    (tmp_path / ".github" / "workflows" / "release.yml").write_text(
        workflow.replace("      - run: make repository-ready\n", ""),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("preflight job must run make repository-ready" in failure for failure in result.failures)


def test_release_readiness_requires_tests_in_preflight(tmp_path):
    checker = _load_checker()
    _write_release_repo(tmp_path)
    workflow = (tmp_path / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")
    (tmp_path / ".github" / "workflows" / "release.yml").write_text(
        workflow
        .replace("      - run: make test\n", "")
        .replace("      - run: make cli-test\n", "")
        .replace("      - run: make python-test\n", "")
        .replace("      - run: python -m pip install build pytest\n", "      - run: python -m pip install build\n"),
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("preflight job must install pytest" in failure for failure in result.failures)
    assert any("preflight job must run make test" in failure for failure in result.failures)
    assert any("preflight job must run make cli-test" in failure for failure in result.failures)
    assert any("preflight job must run make python-test" in failure for failure in result.failures)
