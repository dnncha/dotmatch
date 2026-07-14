#!/usr/bin/env python3
"""Fail closed when the AssayCode Bioconda metapackage drifts."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "packaging" / "bioconda" / "assaycode-meta.yaml"
PYPROJECT = ROOT / "pyproject.toml"


def _version(text: str) -> str | None:
    match = re.search(r'^version\\s*=\\s*"([^"]+)"', text, flags=re.MULTILINE)
    return match.group(1) if match else None


def audit() -> list[str]:
    errors: list[str] = []
    try:
        meta = META.read_text(encoding="utf-8")
        project = PYPROJECT.read_text(encoding="utf-8")
    except OSError as exc:
        return [str(exc)]

    project_version = _version(project)
    recipe_match = re.search(r'{%\\s*set\\s+version\\s*=\\s*"([^"]+)"\\s*%}', meta)
    recipe_version = recipe_match.group(1) if recipe_match else None
    if not project_version or recipe_version != project_version:
        errors.append(
            f"version mismatch: pyproject={project_version!r}, assaycode recipe={recipe_version!r}"
        )

    required = [
        "name: assaycode",
        "noarch: generic",
        "number: 0",
        "- dotmatch =={{ version }}",
        "assaycode --version",
        "dotmatch --version",
        "import assaycode, dotmatch; assert assaycode.engine is dotmatch",
        "license: Apache-2.0",
        "recipe-maintainers:",
        "- dnncha",
    ]
    for fragment in required:
        if fragment not in meta:
            errors.append(f"missing required text: {fragment}")

    forbidden = [
        "REPLACE_WITH_RELEASE_TARBALL_SHA256",
        "source:",
        "build.sh",
        "pip install",
    ]
    for fragment in forbidden:
        if fragment in meta:
            errors.append(f"metapackage must not contain: {fragment}")

    run_match = re.search(r"requirements:\\n  run:\\n(?P<body>(?:    - .+\\n)+)", meta)
    dependencies = [] if run_match is None else [
        line.strip()[2:] for line in run_match.group("body").splitlines()
    ]
    if dependencies != ["dotmatch =={{ version }}"]:
        errors.append(f"metapackage must have exactly one pinned dependency; saw {dependencies}")
    return errors


def main() -> int:
    errors = audit()
    if errors:
        for error in errors:
            print(f"ASSAYCODE BIOCONDA: FAIL: {error}")
        return 1
    print("ASSAYCODE BIOCONDA: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
