from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_bioconda_recipe.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_bioconda_recipe", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _meta(version: str = "0.1.0") -> str:
    return (
        '{% set name = "dotmatch" %}\n'
        f'{{% set version = "{version}" %}}\n'
        '{% set sha256 = "REPLACE_WITH_RELEASE_TARBALL_SHA256" %}\n\n'
        "package:\n"
        "  name: {{ name|lower }}\n"
        "  version: {{ version }}\n\n"
        "source:\n"
        "  url: https://github.com/dnncha/dotmatch/archive/refs/tags/v{{ version }}.tar.gz\n"
        "  sha256: {{ sha256 }}\n\n"
        "build:\n"
        "  number: 0\n"
        "  run_exports:\n"
        "    - {{ pin_subpackage(\"dotmatch\", max_pin=\"x.x\") }}\n"
        "  skip: true  # [win]\n\n"
        "requirements:\n"
        "  build:\n"
        "    - {{ compiler('c') }}\n"
        "    - {{ stdlib('c') }}\n"
        "    - make\n"
        "  host:\n"
        "    - zlib\n\n"
        "test:\n"
        "  commands:\n"
        "    - dotmatch dist ACGT AGGT | grep '^1$'\n"
        "    - dotmatch leq 1 ACGT AGGT | grep '^true$'\n\n"
        "about:\n"
        "  home: https://github.com/dnncha/dotmatch\n"
        "  license: Apache-2.0\n"
        "  license_file: LICENSE\n"
        "  summary: Fast exact short-DNA known-target assignment\n"
        "  description: |\n"
        "    DotMatch is a deterministic known-target short-DNA assignment engine for\n"
        "    CRISPR guides, barcodes, primers, panels, and whitelist-style target sets.\n"
        "    It is not a genome aligner and does not emit SAM/BAM or CIGAR output.\n"
        "  dev_url: https://github.com/dnncha/dotmatch\n"
        "  doc_url: https://github.com/dnncha/dotmatch#readme\n\n"
        "extra:\n"
        "  recipe-maintainers:\n"
        "    - dnncha\n"
    )


def _build() -> str:
    return (
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        "make \\\n"
        '  CC="${CC}" \\\n'
        '  CFLAGS="${CFLAGS:-} ${CPPFLAGS:-} -std=c11 -Wall -Wextra -Wpedantic -Iinclude" \\\n'
        '  LDFLAGS="${LDFLAGS:-}" \\\n'
        "  dotmatch libdotmatch.a shared\n\n"
        'mkdir -p "${PREFIX}/bin" \\\n'
        '         "${PREFIX}/include" \\\n'
        '         "${PREFIX}/lib" \\\n'
        '         "${PREFIX}/share/${PKG_NAME}"\n\n'
        'install -m 755 dotmatch "${PREFIX}/bin/dotmatch"\n'
        'install -m 644 include/qdalign.h "${PREFIX}/include/qdalign.h"\n'
        'install -m 644 libdotmatch.a "${PREFIX}/lib/libdotmatch.a"\n'
        'install -m 644 LICENSE "${PREFIX}/share/${PKG_NAME}/LICENSE"\n\n'
        'if [[ "$(uname -s)" == "Darwin" ]]; then\n'
        '    install -m 755 libdotmatch.dylib "${PREFIX}/lib/libdotmatch.dylib"\n'
        "else\n"
        '    install -m 755 libdotmatch.so "${PREFIX}/lib/libdotmatch.so"\n'
        "fi\n"
    )


def _write_repo(
    root: Path,
    *,
    pyproject_version: str = "0.1.0",
    meta: str | None = None,
    build: str | None = None,
) -> None:
    files = {
        "pyproject.toml": f'[project]\nname = "dotmatch"\nversion = "{pyproject_version}"\n',
        "packaging/bioconda/meta.yaml": meta or _meta(),
        "packaging/bioconda/build.sh": build or _build(),
        "Makefile": "dotmatch:\n\ttrue\nshared:\n\ttrue\nlibdotmatch.a:\n\ttrue\n",
    }
    for path, text in files.items():
        full = root / path
        full.parent.mkdir(parents=True, exist_ok=True)
        full.write_text(text, encoding="utf-8")


def test_bioconda_recipe_checker_exists() -> None:
    assert CHECKER.is_file()


def test_bioconda_recipe_accepts_valid_release_template(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)

    result = checker.audit(tmp_path)

    assert result.failures == []
    assert any("Bioconda recipe" in item for item in result.passed)


def test_bioconda_recipe_rejects_version_mismatch(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path, meta=_meta(version="0.2.0"))

    result = checker.audit(tmp_path)

    assert any("version mismatch" in failure for failure in result.failures)


def test_bioconda_recipe_rejects_resolved_sha_before_release_tarball(tmp_path):
    checker = _load_checker()
    _write_repo(
        tmp_path,
        meta=_meta().replace("REPLACE_WITH_RELEASE_TARBALL_SHA256", "0" * 64),
    )

    result = checker.audit(tmp_path)

    assert any("SHA256 placeholder" in failure for failure in result.failures)


def test_bioconda_recipe_requires_cli_smoke_commands(tmp_path):
    checker = _load_checker()
    meta = (
        _meta()
        .replace("    - dotmatch dist ACGT AGGT | grep '^1$'\n", "")
        .replace("    - dotmatch leq 1 ACGT AGGT | grep '^true$'\n", "")
    )
    _write_repo(tmp_path, meta=meta)

    result = checker.audit(tmp_path)

    assert any("dotmatch dist ACGT AGGT" in failure for failure in result.failures)
    assert any("dotmatch leq 1 ACGT AGGT" in failure for failure in result.failures)


def test_bioconda_recipe_rejects_broad_genome_aligner_claim(tmp_path):
    checker = _load_checker()
    meta = _meta().replace(
        "    It is not a genome aligner and does not emit SAM/BAM or CIGAR output.\n",
        "    It is a genome aligner for short sequencing reads.\n",
    )
    _write_repo(tmp_path, meta=meta)

    result = checker.audit(tmp_path)

    assert any("not a genome aligner" in failure for failure in result.failures)


def test_bioconda_recipe_requires_native_install_steps(tmp_path):
    checker = _load_checker()
    build = (
        _build()
        .replace('install -m 755 dotmatch "${PREFIX}/bin/dotmatch"\n', "")
        .replace('install -m 644 include/qdalign.h "${PREFIX}/include/qdalign.h"\n', "")
        .replace('install -m 644 libdotmatch.a "${PREFIX}/lib/libdotmatch.a"\n', "")
        .replace('install -m 644 LICENSE "${PREFIX}/share/${PKG_NAME}/LICENSE"\n', "")
        .replace('install -m 755 libdotmatch.dylib "${PREFIX}/lib/libdotmatch.dylib"\n', "")
        .replace('install -m 755 libdotmatch.so "${PREFIX}/lib/libdotmatch.so"\n', "")
    )
    _write_repo(tmp_path, build=build)

    result = checker.audit(tmp_path)

    assert any("dotmatch" in failure for failure in result.failures)
    assert any("qdalign.h" in failure for failure in result.failures)
    assert any("libdotmatch.a" in failure for failure in result.failures)
    assert any("LICENSE" in failure for failure in result.failures)
    assert any("libdotmatch.dylib" in failure for failure in result.failures)
    assert any("libdotmatch.so" in failure for failure in result.failures)
