"""Create a NEW standalone GitHub repository from the checked source inventory.

Requires an already authenticated local GitHub CLI. Does not publish to PyPI,
modify an existing repository, push DotMatch history, or ask for credentials.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from release_manifest import release_manifest


def command(*args: str, cwd: Path | None = None, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(list(args), cwd=cwd, text=True, capture_output=True, check=check)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    visibility = parser.add_mutually_exclusive_group(required=True)
    visibility.add_argument("--public", action="store_true")
    visibility.add_argument("--private", action="store_true")
    parser.add_argument("--owner", default="dnncha")
    parser.add_argument("--dry-run", action="store_true", help="Check sources and show the plan; no network or writes")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9-]{0,38}", args.owner):
        parser.error("invalid GitHub user name")
    root = Path(__file__).resolve().parents[1]
    inventory_path = root / "release-files.json"
    if not inventory_path.is_file():
        parser.error("release-files.json is missing; review sources before generating it")
    inventory = json.loads(inventory_path.read_text())
    if inventory != release_manifest(root):
        parser.error("source inventory mismatch; review changes before regenerating release-files.json")
    repo = f"{args.owner}/editwitness"
    plan = {"repository": repo, "visibility": "public" if args.public else "private",
            "source_files": len(inventory["files"]), "fresh_history": True,
            "existing_repository_policy": "refuse", "publishes_to_pypi": False}
    if args.dry_run:
        print(json.dumps(plan, indent=2))
        return 0
    for program in ("git", "gh"):
        if shutil.which(program) is None:
            parser.error(f"{program} is required locally; this script cannot create GitHub repositories without it")
    work: Path | None = None
    try:
        command("gh", "auth", "status")
        user = json.loads(command("gh", "api", "user").stdout)
        if user["login"].lower() != args.owner.lower():
            parser.error("authenticated GitHub user does not match --owner; no repository was created")
        exists = command("gh", "api", f"repos/{repo}", "--silent", check=False)
        if exists.returncode == 0:
            parser.error(f"{repo} already exists; refusing to modify or overwrite it")
        if "404" not in exists.stderr:
            parser.error("repository preflight failed for a reason other than Not Found; no writes attempted")
        work = Path(tempfile.mkdtemp(prefix="editwitness-publish-"))
        for relative in inventory["files"]:
            destination = work / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(root / relative, destination)
        shutil.copy2(inventory_path, work / inventory_path.name)
        command("git", "init", "-b", "main", cwd=work)
        command("git", "add", ".", cwd=work)
        command("git", "-c", f"user.name={user.get('name') or user['login']}", "-c",
                f"user.email={user['id']}+{user['login']}@users.noreply.github.com", "commit", "-m",
                "Release EditWitness 0.1.0a1 research alpha", cwd=work)
        created = command("gh", "repo", "create", repo, "--public" if args.public else "--private",
                          "--description", "CRISPR assay observability: explicit counterexamples, candidate panels, local agent-friendly evidence.",
                          "--source", str(work), "--remote", "origin", "--push", cwd=work)
        print(created.stdout.strip())
        print(f"Published a fresh repository: https://github.com/{repo}")
        print("No PyPI publication, clinical validation, or passing remote CI is implied.")
        shutil.rmtree(work)
        return 0
    except subprocess.CalledProcessError as error:
        print(error.stderr or str(error))
        if work is not None:
            print(f"Recovery source retained at {work}. Check whether the repository was created before retrying.")
            print("The script never deletes a remote repository or force-pushes an existing branch.")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
