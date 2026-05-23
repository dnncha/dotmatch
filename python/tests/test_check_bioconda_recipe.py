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
        "  skip: true  # [win or py<39]\n\n"
        "requirements:\n"
        "  build:\n"
        "    - {{ compiler('c') }}\n"
        "    - {{ stdlib('c') }}\n"
        "    - make\n"
        "  host:\n"
        "    - python >=3.9\n"
        "    - pip\n"
        "    - setuptools >=77\n"
        "    - wheel\n"
        "    - zlib\n"
        "  run:\n"
        "    - python >=3.9\n"
        "    - tomli  # [py<311]\n\n"
        "test:\n"
        "  commands:\n"
        "    - python -c \"import dotmatch; assert dotmatch.distance('ACGT', 'AGGT') == 1\"\n"
        "    - python -c \"from dotmatch.native import find_native_cli; p=find_native_cli(); assert p.name == 'dotmatch-native' and p.exists()\"\n"
        "    - dotmatch --version | grep '^dotmatch {{ version }}$'\n"
        "    - dotmatch dist ACGT AGGT | grep '^1$'\n"
        "    - dotmatch leq 1 ACGT AGGT | grep '^true$'\n"
        "    - dotmatch --help | grep 'Workflow namespaces:'\n"
        "    - dotmatch assay --help | grep 'dotmatch assay'\n"
        "    - dotmatch barcode --help | grep 'dotmatch barcode'\n"
        "    - dotmatch panel --help | grep 'dotmatch panel'\n"
        "    - test -f \"${PREFIX}/include/qdalign.h\"\n"
        "    - test -f \"${PREFIX}/lib/libdotmatch.a\"\n"
        "    - test -f \"${PREFIX}/lib/libdotmatch.so\" || test -f \"${PREFIX}/lib/libdotmatch.dylib\"\n"
        "    - dotmatch assay init --template crispr --out assay.toml\n"
        "    - grep 'assay_type = \"crispr\"' assay.toml\n"
        "    - printf 'target_id\\ttarget_seq\\nbc0\\tACGT\\n' > targets.tsv\n"
        "    - printf '@r0\\nACGT\\n+\\nIIII\\n' > reads.fastq\n"
        "    - dotmatch count --targets targets.tsv --reads reads.fastq --sample-label sample --target-start 0 --target-length 4 --k 0 --metric hamming --out counts.tsv\n"
        "    - awk -F '\\t' 'NR==2 { exit !($1==\"bc0\" && $2==\"ACGT\" && $3==\"\" && $4==\"0\" && $5==\"1\" && $10==\"1\") }' counts.tsv\n\n"
        "    - printf 'barcode_id\\tbarcode_seq\\ns1\\tACGT\\ns2\\tTTTT\\n' > barcodes.tsv\n"
        "    - printf '@r1\\nNACGTAAAA\\n+\\nIIIIIIIII\\n@r2\\nNTTTTAAAA\\n+\\nIIIIIIIII\\n' > barcode_reads.fastq\n"
        "    - dotmatch barcode infer --barcodes barcodes.tsv --reads barcode_reads.fastq --scan-starts 0:2 --barcode-length 4 --sample-reads 10 --out offset_scan.tsv --summary barcode_summary.json\n"
        "    - grep '\"recommended_start\": 1' barcode_summary.json\n"
        "    - dotmatch panel design --n 2 --length 4 --candidate-pool-size 100 --restarts 1 --min-hamming-distance 2 --min-levenshtein-distance 2 --out-dir panel_out\n"
        "    - test -f panel_out/barcodes.tsv\n"
        "    - test -f panel_out/design_report.json\n\n"
        "about:\n"
        "  home: https://github.com/dnncha/dotmatch\n"
        "  license: Apache-2.0\n"
        "  license_file: LICENSE\n"
        "  summary: Fast exact short-DNA known-target assignment\n"
        "  description: |\n"
        "    DotMatch provides deterministic known-target short-DNA assignment for\n"
        "    CRISPR guides, barcodes, primers, panels, and whitelist-style target sets.\n"
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
        '  DOTMATCH_VERSION="${PKG_VERSION}" \\\n'
        '  CC="${CC}" \\\n'
        '  CFLAGS="${CFLAGS:-} ${CPPFLAGS:-} -std=c11 -Wall -Wextra -Wpedantic -Iinclude" \\\n'
        '  LDFLAGS="${LDFLAGS:-}" \\\n'
        "  libdotmatch.a shared\n\n"
        'mkdir -p "${PREFIX}/bin" \\\n'
        '         "${PREFIX}/include" \\\n'
        '         "${PREFIX}/lib" \\\n'
        '         "${PREFIX}/share/${PKG_NAME}"\n\n'
        '${PYTHON} -m pip install . -vv --no-deps --no-build-isolation\n\n'
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


def test_bioconda_recipe_requires_native_install_steps(tmp_path):
    checker = _load_checker()
    build = (
        _build()
        .replace('${PYTHON} -m pip install . -vv --no-deps --no-build-isolation\n\n', "")
        .replace('install -m 644 include/qdalign.h "${PREFIX}/include/qdalign.h"\n', "")
        .replace('install -m 644 libdotmatch.a "${PREFIX}/lib/libdotmatch.a"\n', "")
        .replace('install -m 644 LICENSE "${PREFIX}/share/${PKG_NAME}/LICENSE"\n', "")
        .replace('install -m 755 libdotmatch.dylib "${PREFIX}/lib/libdotmatch.dylib"\n', "")
        .replace('install -m 755 libdotmatch.so "${PREFIX}/lib/libdotmatch.so"\n', "")
    )
    _write_repo(tmp_path, build=build)

    result = checker.audit(tmp_path)

    assert any("Python console script" in failure for failure in result.failures)
    assert any("qdalign.h" in failure for failure in result.failures)
    assert any("libdotmatch.a" in failure for failure in result.failures)
    assert any("LICENSE" in failure for failure in result.failures)
    assert any("libdotmatch.dylib" in failure for failure in result.failures)
    assert any("libdotmatch.so" in failure for failure in result.failures)
