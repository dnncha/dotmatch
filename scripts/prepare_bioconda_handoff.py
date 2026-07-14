#!/usr/bin/env python3
"""Render release-ready DotMatch and AssayCode Bioconda recipe directories."""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PLACEHOLDER = "REPLACE_WITH_RELEASE_TARBALL_SHA256"


def _project_version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    if match is None:
        raise ValueError("pyproject.toml does not declare a version")
    return match.group(1)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def render(release_tarball: Path, output: Path) -> tuple[Path, Path, str]:
    if not release_tarball.is_file():
        raise FileNotFoundError(release_tarball)
    version = _project_version()
    accepted_names = {f"v{version}.tar.gz", f"dotmatch-{version}.tar.gz", f"{version}.tar.gz"}
    if release_tarball.name not in accepted_names:
        raise ValueError(
            f"release archive name must identify {version}; expected one of {sorted(accepted_names)}"
        )

    digest = _sha256(release_tarball)
    dotmatch_dir = output / "recipes" / "dotmatch"
    assaycode_dir = output / "recipes" / "assaycode"
    dotmatch_dir.mkdir(parents=True, exist_ok=True)
    assaycode_dir.mkdir(parents=True, exist_ok=True)

    dotmatch_template = (ROOT / "packaging" / "bioconda" / "meta.yaml").read_text(
        encoding="utf-8"
    )
    if PLACEHOLDER not in dotmatch_template:
        raise ValueError("DotMatch recipe template has no release checksum placeholder")
    (dotmatch_dir / "meta.yaml").write_text(
        dotmatch_template.replace(PLACEHOLDER, digest), encoding="utf-8"
    )
    shutil.copy2(ROOT / "packaging" / "bioconda" / "build.sh", dotmatch_dir / "build.sh")
    shutil.copy2(
        ROOT / "packaging" / "bioconda" / "assaycode-meta.yaml",
        assaycode_dir / "meta.yaml",
    )
    return dotmatch_dir, assaycode_dir, digest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-tarball", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    try:
        dotmatch_dir, assaycode_dir, digest = render(args.release_tarball, args.out)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(f"sha256={digest}")
    print(f"dotmatch_recipe={dotmatch_dir}")
    print(f"assaycode_recipe={assaycode_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
