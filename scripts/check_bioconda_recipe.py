#!/usr/bin/env python3
"""Audit the Bioconda pre-submission recipe template."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


META = Path("packaging") / "bioconda" / "meta.yaml"
BUILD = Path("packaging") / "bioconda" / "build.sh"
PYPROJECT = Path("pyproject.toml")
SHA_PLACEHOLDER = "REPLACE_WITH_RELEASE_TARBALL_SHA256"


class AuditResult:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def _read(path: Path, result: AuditResult) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except Exception as exc:
        result.failures.append(f"{path.as_posix()} could not be read: {exc}")
        return ""


def _pyproject_version(text: str) -> str | None:
    in_project = False
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            in_project = line == "[project]"
            continue
        if not in_project or not line.startswith("version"):
            continue
        match = re.match(r'version\s*=\s*["\']([^"\']+)["\']', line)
        if match:
            return match.group(1)
    return None


def _recipe_version(text: str) -> str | None:
    match = re.search(r'{%\s*set\s+version\s*=\s*"([^"]+)"\s*%}', text)
    return match.group(1) if match else None


def _check_versions(pyproject: str, meta: str, result: AuditResult) -> None:
    project_version = _pyproject_version(pyproject)
    recipe_version = _recipe_version(meta)
    if not project_version:
        result.failures.append("pyproject.toml must declare [project] version")
    if not recipe_version:
        result.failures.append("Bioconda recipe must declare Jinja version")
    if project_version and recipe_version and project_version != recipe_version:
        result.failures.append(
            f"Bioconda recipe version mismatch: pyproject.toml={project_version}, "
            f"packaging/bioconda/meta.yaml={recipe_version}"
        )


def _require(text: str, fragment: str, message: str, result: AuditResult) -> None:
    if fragment not in text:
        result.failures.append(message)


def _check_meta(meta: str, result: AuditResult) -> None:
    required_fragments = [
        ('{% set name = "dotmatch" %}', "Bioconda recipe must set package name to dotmatch"),
        (
            "https://github.com/dnncha/dotmatch/archive/refs/tags/v{{ version }}.tar.gz",
            "Bioconda recipe must use the immutable GitHub release tarball URL",
        ),
        ("sha256: {{ sha256 }}", "Bioconda recipe must wire source sha256 through the Jinja sha256 variable"),
        ("skip: true  # [win]", "Bioconda recipe must skip unsupported Windows builds"),
        ("{{ compiler('c') }}", "Bioconda recipe must request the C compiler"),
        ("- make", "Bioconda recipe must include make in build requirements"),
        ("- zlib", "Bioconda recipe must include zlib requirements"),
        ("license: Apache-2.0", "Bioconda recipe must declare Apache-2.0 license"),
        ("license_file: LICENSE", "Bioconda recipe must install and declare LICENSE"),
        ("recipe-maintainers:", "Bioconda recipe must declare recipe maintainers"),
        ("- dnncha", "Bioconda recipe must list dnncha as recipe maintainer"),
    ]
    for fragment, message in required_fragments:
        _require(meta, fragment, message, result)

    if SHA_PLACEHOLDER not in meta:
        result.failures.append("Bioconda recipe must retain SHA256 placeholder before submission")

    for command in [
        "dotmatch --version",
        "dotmatch dist ACGT AGGT",
        "dotmatch leq 1 ACGT AGGT",
    ]:
        if command not in meta:
            result.failures.append(f"Bioconda recipe test commands must include {command}")

    lower_meta = meta.lower()
    if "not a genome aligner" not in lower_meta:
        result.failures.append("Bioconda recipe description must state DotMatch is not a genome aligner")
    if re.search(r"\b(accepted|published|released|available)\s+(on|in|from)\s+bioconda\b", meta, flags=re.I):
        result.failures.append("Bioconda recipe must not claim public Bioconda availability before submission")


def _check_build(build: str, result: AuditResult) -> None:
    required_fragments = [
        ("set -euo pipefail", "Bioconda build.sh must fail fast with set -euo pipefail"),
        ("make dotmatch libdotmatch.a shared", "Bioconda build.sh must build CLI, static library, and shared library"),
        ('"${PREFIX}/bin"', "Bioconda build.sh must create the bin install directory"),
        ('"${PREFIX}/include"', "Bioconda build.sh must create the include install directory"),
        ('"${PREFIX}/lib"', "Bioconda build.sh must create the lib install directory"),
        ('"${PREFIX}/share/${PKG_NAME}"', "Bioconda build.sh must create the package share directory"),
        ('install -m 755 dotmatch "${PREFIX}/bin/dotmatch"', "Bioconda build.sh must install dotmatch"),
        ('install -m 644 include/qdalign.h "${PREFIX}/include/qdalign.h"', "Bioconda build.sh must install qdalign.h"),
        ('install -m 644 libdotmatch.a "${PREFIX}/lib/libdotmatch.a"', "Bioconda build.sh must install libdotmatch.a"),
        ('install -m 644 LICENSE "${PREFIX}/share/${PKG_NAME}/LICENSE"', "Bioconda build.sh must install LICENSE"),
        ('libdotmatch.dylib "${PREFIX}/lib/libdotmatch.dylib"', "Bioconda build.sh must install libdotmatch.dylib on Darwin"),
        ('libdotmatch.so "${PREFIX}/lib/libdotmatch.so"', "Bioconda build.sh must install libdotmatch.so on Linux"),
    ]
    for fragment, message in required_fragments:
        _require(build, fragment, message, result)
    if "uname -s" not in build or "Darwin" not in build:
        result.failures.append("Bioconda build.sh must branch shared-library install by platform")


def audit(root: Path) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    pyproject = _read(root / PYPROJECT, result)
    meta = _read(root / META, result)
    build = _read(root / BUILD, result)
    if not pyproject or not meta or not build:
        return result

    _check_versions(pyproject, meta, result)
    _check_meta(meta, result)
    _check_build(build, result)

    if not result.failures:
        result.passed.append("Bioconda recipe pre-submission template is ready")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="repository root")
    args = parser.parse_args()

    result = audit(Path(args.root))
    for item in result.passed:
        print(f"PASS: {item}")
    for item in result.failures:
        print(f"FAIL: {item}")
    if result.ok:
        print("BIOCONDA RECIPE: PASS")
        return 0
    print("BIOCONDA RECIPE: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
