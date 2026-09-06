"""Publish a standalone repository, optionally a CI-gated GitHub research prerelease.

Requires an already authenticated GitHub CLI. Never asks for credentials,
force-pushes, changes visibility of an existing repository, or publishes to PyPI.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
import tomllib
from pathlib import Path

from release_manifest import release_manifest


def command(*args: str, cwd: Path | None = None, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(list(args), cwd=cwd, text=True, encoding="utf-8", capture_output=True,
                          check=check, timeout=1800, env=dict(os.environ, GIT_TERMINAL_PROMPT="0"))


def checked_version(root: Path) -> str:
    version = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
    if not re.fullmatch(r"\d+\.\d+\.\d+a\d+", version):
        raise ValueError("this publisher releases explicitly alpha versions only")
    return str(version)


def verify_artifacts(directory: Path, *, expected_version: str | None = None) -> list[Path]:
    """Validate every CI-produced asset; disallow paths and extra unlisted files."""
    checksums = directory / "SHA256SUMS"
    expected = {}
    for line in checksums.read_text(encoding="ascii").splitlines():
        match = re.fullmatch(r"([a-f0-9]{64})  ([A-Za-z0-9_.-]+)", line)
        if match is None:
            raise ValueError("invalid distribution checksum entry")
        sha, name = match.groups()
        if name in expected or name == "SHA256SUMS":
            raise ValueError("duplicate or recursive distribution checksum entry")
        expected[name] = sha
    wheels = [name for name in expected if name.endswith(".whl")]
    sources = [name for name in expected if name.endswith(".tar.gz")]
    if len(expected) != 2 or len(wheels) != 1 or len(sources) != 1:
        raise ValueError("CI artifacts must contain exactly one checked wheel and source distribution")
    if expected_version is not None and (
        wheels[0] != f"editwitness-{expected_version}-py3-none-any.whl"
        or sources[0] != f"editwitness-{expected_version}.tar.gz"
    ):
        raise ValueError("CI artifact names do not match the reviewed package version")
    actual = {p.name for p in directory.iterdir()}
    if actual != set(expected) | {"SHA256SUMS"}:
        raise ValueError("unlisted or missing CI distribution assets")
    for name, sha in expected.items():
        path = directory / name
        if path.is_symlink() or not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != sha:
            raise ValueError(f"distribution checksum mismatch: {name}")
    return [directory / name for name in sorted(expected)] + [checksums]


def wait_for_ci(repo: str, sha: str) -> int:
    deadline = time.monotonic() + 120
    while time.monotonic() < deadline:
        response = command("gh", "run", "list", "--repo", repo, "--commit", sha,
            "--workflow", "ci.yml", "--limit", "10", "--json", "databaseId,headSha,event")
        runs = [r for r in json.loads(response.stdout) if r["headSha"] == sha and r["event"] == "push"]
        if runs:
            run_id = int(runs[0]["databaseId"])
            command("gh", "run", "watch", str(run_id), "--repo", repo, "--exit-status", "--interval", "5")
            return run_id
        time.sleep(3)
    raise ValueError("No matching main-branch CI run appeared. Source may be published, but no release was created.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    visibility = parser.add_mutually_exclusive_group(required=True)
    visibility.add_argument("--public", action="store_true")
    visibility.add_argument("--private", action="store_true")
    parser.add_argument("--owner", default="dnncha")
    parser.add_argument("--release", action="store_true", help="Wait for matching CI; publish its verified assets as a prerelease")
    parser.add_argument("--resume", action="store_true", help="Resume only when the existing repo has exactly this reviewed source")
    parser.add_argument("--dry-run", action="store_true", help="Verify inventory and print a plan; no network or changes")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9-]{0,38}", args.owner):
        parser.error("invalid GitHub user name")
    root = Path(__file__).resolve().parents[1]
    inventory_path = root / "release-files.json"
    try:
        inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
        if inventory != release_manifest(root):
            raise ValueError("source inventory mismatch; review edits before regenerating release-files.json")
        version = checked_version(root)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    repo = f"{args.owner}/editwitness"
    if args.dry_run:
        print(json.dumps({"repository": repo, "visibility": "public" if args.public else "private",
            "version": version, "source_files": len(inventory["files"]), "fresh_history": not args.resume,
            "release_after_matching_ci": args.release, "publishes_to_pypi": False,
            "resume_policy": "exact source and visibility match; never overwrite"}, indent=2))
        return 0
    for program in ("git", "gh"):
        if shutil.which(program) is None:
            parser.error(f"{program} is required and unavailable here; no remote writes were attempted")
    work: Path | None = None
    try:
        command("gh", "auth", "status")
        user = json.loads(command("gh", "api", "user").stdout)
        if user["login"].lower() != args.owner.lower():
            raise ValueError("authenticated account does not match --owner; no writes attempted")
        check = command("gh", "api", f"repos/{repo}", check=False)
        exists = check.returncode == 0
        if not exists and "404" not in check.stderr:
            raise ValueError("repository lookup failed for a reason other than Not Found; no writes attempted")
        if exists and not args.resume:
            raise ValueError("repository exists; use --resume only for an exact reviewed-source match")
        work = Path(tempfile.mkdtemp(prefix="editwitness-publish-"))
        identity = ("-c", f"user.name={user.get('name') or user['login']}", "-c",
                    f"user.email={user['id']}+{user['login']}@users.noreply.github.com")
        if exists:
            remote = json.loads(check.stdout)
            if remote["private"] != args.private:
                raise ValueError("existing repository visibility differs; refusing to change it")
            clone = work / "source"
            command("gh", "repo", "clone", repo, str(clone), "--", "--branch", "main", "--single-branch")
            work = clone
            tracked = set(command("git", "ls-files", "-z", cwd=work).stdout.split("\0")) - {""}
            expected_files = set(inventory["files"]) | {"release-files.json"}
            if tracked != expected_files:
                raise ValueError("remote tracked files differ from the reviewed source set")
            remote_inventory = json.loads((work / "release-files.json").read_text(encoding="utf-8"))
            if remote_inventory != inventory or release_manifest(work) != inventory:
                raise ValueError("remote source differs from reviewed inventory; refusing to overwrite it")
        else:
            for relative in inventory["files"]:
                destination = work / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(root / relative, destination)
            shutil.copy2(inventory_path, work / inventory_path.name)
            command("git", "init", "-b", "main", cwd=work)
            command("git", "add", ".", cwd=work)
            command("git", *identity, "commit", "-m", f"EditWitness {version}: audited research alpha", cwd=work)
            command("gh", "repo", "create", repo, "--public" if args.public else "--private",
                "--description", "CRISPR assay blind spots: sequence-level counterexamples, primer rematching and validation design.",
                "--source", str(work), "--remote", "origin", "--push", cwd=work)
        sha = command("git", "rev-parse", "HEAD", cwd=work).stdout.strip()
        print(f"Source confirmed at https://github.com/{repo}/commit/{sha}", flush=True)
        if args.release:
            run_id = wait_for_ci(repo, sha)
            dist = work / "ci-distributions"
            command("gh", "run", "download", str(run_id), "--repo", repo,
                "--name", "editwitness-distributions", "--dir", str(dist))
            assets = verify_artifacts(dist, expected_version=version)
            tag = f"v{version}"
            existing_tag = command("git", "ls-remote", "--tags", "origin", f"refs/tags/{tag}", cwd=work).stdout.strip()
            if existing_tag:
                raise ValueError("release tag already exists; no tag or release was overwritten")
            command("git", *identity, "tag", "-a", tag, "-m", f"EditWitness {version} research alpha", cwd=work)
            command("git", "push", "origin", f"refs/tags/{tag}", cwd=work)
            command("gh", "release", "create", tag, "--repo", repo, "--verify-tag", "--prerelease",
                "--title", f"EditWitness {version}", "--notes-file", str(work / "RELEASE_NOTES.md"),
                *map(str, assets), cwd=work)
            print(f"GitHub prerelease confirmed: https://github.com/{repo}/releases/tag/{tag}", flush=True)
        print("No PyPI publication or empirical biological validation is implied.")
        return 0
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, ValueError, OSError) as error:
        print(getattr(error, "stderr", None) or str(error))
        if work is not None:
            print(f"Recovery source retained at {work}; inspect the remote before retrying.")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
