"""List and hash distributable source files; omit local data, caches and build products."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT_FILES = {
    "README.md", "LICENSE", "NOTICE", "CITATION.cff", "CHANGELOG.md", "CONTRIBUTING.md",
    "SECURITY.md", "AGENTS.md", "BUILD_STATUS.md", "llms.txt", "roadmap.json", "pyproject.toml",
    ".gitignore", "MANIFEST.in",
}
DIRECTORIES = {"src", "tests", "docs", "examples", "scripts", "benchmarks", "skills", ".github"}
ALLOWED_SUFFIXES = {".py", ".json", ".md", ".yml", ".yaml", ".toml", ".fasta", ".txt"}


def source_files(root: Path) -> list[Path]:
    files = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        relative = path.relative_to(root)
        if relative.as_posix() in {".github/workflows/bootstrap-editwitness.yml",
                                   ".github/workflows/public-editwitness.yml"}:
            continue  # Temporary transport-only workflow, never part of the standalone package.
        if any(part in {"__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache"}
               or part.endswith(".egg-info") for part in relative.parts):
            continue
        included = (len(relative.parts) == 1 and relative.name in ROOT_FILES) or (
            relative.parts[0] in DIRECTORIES and (path.suffix in ALLOWED_SUFFIXES or path.name == "py.typed")
        )
        if included:
            if path.is_symlink():
                raise ValueError(f"refusing source symlink: {relative}")
            files.append(path)
    return sorted(files)


def release_manifest(root: Path) -> dict[str, object]:
    return {
        "format": "editwitness.source-manifest.v1",
        "files": {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in source_files(root)},
        "note": "Integrity inventory, not an authenticated signature. Regenerate only after reviewed changes.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    target = root / "release-files.json"
    expected = release_manifest(root)
    if args.check:
        if not target.is_file() or json.loads(target.read_text(encoding="utf-8")) != expected:
            parser.exit(1, "Source inventory differs; inspect changes before regenerating.\n")
        print("Source inventory verified.")
    else:
        target.write_text(json.dumps(expected, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("Wrote release-files.json.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
