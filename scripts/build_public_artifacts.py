"""Build an integrity-checked public download bundle from an already tested source.

This script writes local artifacts only. It does not create GitHub repositories,
release tags, registry publications or scientific validation claims.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

from release_manifest import release_manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--run-url", required=True)
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    import tomllib
    version = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
    inventory = json.loads((root / "release-files.json").read_text(encoding="utf-8"))
    if inventory != release_manifest(root):
        raise RuntimeError("Source inventory mismatch; refusing to distribute unchecked files")
    destination = args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=False)  # Never replace an already published version.
    subprocess.run([sys.executable, "-m", "build"], cwd=root, check=True)
    for source in sorted((root / "dist").glob(f"editwitness-{version}*")):
        if source.suffix == ".whl" or source.name.endswith(".tar.gz"):
            shutil.copy2(source, destination / source.name)
    archive = destination / f"editwitness-{version}-source.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as handle:
        for name in sorted([*inventory["files"], "release-files.json"]):
            info = zipfile.ZipInfo(f"editwitness/{name}", date_time=(2026, 9, 6, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            handle.writestr(info, (root / name).read_bytes())
    def cli(*arguments: str) -> None:
        subprocess.run([sys.executable, "-m", "editwitness", *arguments], cwd=root, check=True)
    cli("demo", "-o", str(destination / "demo.json"))
    cli("analyze", str(destination / "demo.json"), "-o", str(destination / "evidence.json"),
        "--html", str(destination / "report.html"))
    cli("verify", str(destination / "evidence.json"), "--manifest", str(destination / "demo.json"))
    hashes = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
              for path in sorted(destination.iterdir()) if path.is_file()}
    receipt = {
        "package": "editwitness", "version": version, "source_commit": args.source_commit,
        "verification_run": args.run_url,
        "artifact_policy": "Version directory is create-only; SHA256 hashes are integrity checks, not authentication.",
        "verification_policy": "Publication job depends on successful OS/Python matrix, strict typing and coverage jobs.",
        "artifacts": hashes,
        "scientific_status": "Software-tested research model; no empirical biological validation.",
        "standalone_repository_created": False, "pypi_published": False,
        "github_prerelease": {"status": "not_attempted"},
    }
    (destination / "publication.json").write_text(json.dumps(receipt, indent=2)+"\n", encoding="utf-8")
    (destination / "SHA256SUMS").write_text("".join(f"{sha}  {name}\n" for name, sha in hashes.items()), encoding="utf-8")
    (destination / "README.md").write_text(
        f"# EditWitness {version}\n\nResearch-alpha downloads built from `{args.source_commit}`.\n\n"
        f"[Verification run]({args.run_url}). Read `publication.json` for the actual publication status.\n\n"
        f"Install the wheel with `python -m pip install ./editwitness-{version}-py3-none-any.whl`. "
        "Verify its SHA-256 against `SHA256SUMS`. The source ZIP has independent package sources, "
        "not DotMatch history or temporary transport code.\n\n"
        "The report and JSON example are explicitly synthetic. No experimental validation, clinical "
        "suitability, standalone repository or PyPI publication is implied.\n", encoding="utf-8")
    print(json.dumps({"public_artifacts": str(destination), "source_commit": args.source_commit,
                      "version": version, "files": sorted(p.name for p in destination.iterdir())}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
