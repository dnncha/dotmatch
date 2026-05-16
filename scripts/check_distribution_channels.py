#!/usr/bin/env python3
"""Verify that a DotMatch release is publicly available on distribution channels.

This is a post-release verifier. It is expected to fail before a tag has been
published to PyPI, Bioconda, GHCR, and Zenodo.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import tempfile
import urllib.error
import urllib.request
import venv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


PYPI_URL = "https://pypi.org/pypi/dotmatch/{version}/json"
BIOCONDA_URL = "https://api.anaconda.org/package/bioconda/dotmatch"
GHCR_IMAGE = "ghcr.io/dnncha/dotmatch:v{version}"
BIOCONTAINERS_TAGS_URL = (
    "https://quay.io/api/v1/repository/biocontainers/dotmatch/tag/?onlyActiveTags=true&page={page}&limit=100"
)
BIOCONTAINERS_IMAGE = "quay.io/biocontainers/dotmatch:{tag}"


@dataclass(frozen=True)
class ChannelMessage:
    channel: str
    message: str


@dataclass
class AuditResult:
    passed: list[ChannelMessage] = field(default_factory=list)
    failures: list[ChannelMessage] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.failures


def fetch_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "DotMatch distribution verifier"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.loads(response.read().decode("utf-8"))


def url_ok(url: str) -> bool:
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "DotMatch distribution verifier"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return 200 <= int(response.status) < 400
    except urllib.error.HTTPError as exc:
        if exc.code not in {405, 501}:
            return False
    request = urllib.request.Request(url, headers={"User-Agent": "DotMatch distribution verifier"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return 200 <= int(response.status) < 400
    except Exception:
        return False


def project_version(root: Path) -> str:
    text = (root / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    if match is None:
        raise ValueError("pyproject.toml does not declare project version")
    return match.group(1)


def citation_doi(root: Path) -> str:
    text = (root / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r'^\s*doi\s*:\s*["\']?([^"\'\s]+)', text, flags=re.MULTILINE)
    return match.group(1).strip() if match else ""


def clean_install_env() -> dict[str, str]:
    env = os.environ.copy()
    env.pop("DOTMATCH_LIB", None)
    env.pop("QUICKDNA_LIB", None)
    env.pop("PYTHONPATH", None)
    return env


def venv_python(env_dir: Path) -> Path:
    if os.name == "nt":
        return env_dir / "Scripts" / "python.exe"
    return env_dir / "bin" / "python"


def venv_script(env_dir: Path, name: str) -> Path:
    if os.name == "nt":
        return env_dir / "Scripts" / f"{name}.exe"
    return env_dir / "bin" / name


def run_checked(cmd: list[str], *, cwd: Path, env: dict[str, str]) -> str:
    proc = subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if proc.returncode != 0:
        detail = (proc.stderr or "").strip() or (proc.stdout or "").strip()
        if detail:
            raise RuntimeError(detail)
        raise RuntimeError(f"command failed with exit code {proc.returncode}: {' '.join(cmd)}")
    return (proc.stdout or "").strip()


def verify_pypi_install(version: str) -> None:
    with tempfile.TemporaryDirectory(prefix="dotmatch-pypi-install-") as tmp:
        root = Path(tmp)
        env_dir = root / "venv"
        venv.EnvBuilder(with_pip=True).create(env_dir)
        py = venv_python(env_dir)
        env = clean_install_env()

        run_checked([str(py), "-m", "pip", "install", "--quiet", f"dotmatch=={version}"], cwd=root, env=env)
        run_checked(
            [
                str(py),
                "-c",
                "import dotmatch; assert dotmatch.distance('ACGT', 'AGGT') == 1; print('import ok')",
            ],
            cwd=root,
            env=env,
        )
        observed_version = run_checked([str(venv_script(env_dir, "dotmatch")), "--version"], cwd=root, env=env)
        expected_version = f"dotmatch {version}"
        if observed_version != expected_version:
            raise RuntimeError(f"dotmatch --version reported {observed_version!r}, expected {expected_version!r}")
        observed_distance = run_checked([str(venv_script(env_dir, "dotmatch")), "dist", "ACGT", "AGGT"], cwd=root, env=env)
        if observed_distance != "1":
            raise RuntimeError(f"dotmatch dist smoke test reported {observed_distance!r}, expected '1'")


def verify_ghcr_run(image: str, version: str) -> None:
    env = os.environ.copy()
    cwd = Path.cwd()
    observed_version = run_checked(["docker", "run", "--rm", image, "--version"], cwd=cwd, env=env)
    expected_version = f"dotmatch {version}"
    if observed_version != expected_version:
        raise RuntimeError(f"docker image --version reported {observed_version!r}, expected {expected_version!r}")
    observed_distance = run_checked(["docker", "run", "--rm", image, "dist", "ACGT", "AGGT"], cwd=cwd, env=env)
    if observed_distance != "1":
        raise RuntimeError(f"docker image dist smoke test reported {observed_distance!r}, expected '1'")


def verify_bioconda_install(version: str) -> None:
    conda = shutil.which("micromamba") or shutil.which("conda")
    if conda is None:
        raise RuntimeError("micromamba or conda is required to verify the Bioconda install")

    with tempfile.TemporaryDirectory(prefix="dotmatch-bioconda-install-") as tmp:
        root = Path(tmp)
        env = os.environ.copy()
        if Path(conda).name == "micromamba":
            env.setdefault("MAMBA_ROOT_PREFIX", str(root / "mamba-root"))
        prefix = root / "env"
        channels = ["-c", "conda-forge", "-c", "bioconda"]
        run_checked([conda, "create", "-y", "-p", str(prefix), *channels, f"dotmatch={version}"], cwd=root, env=env)
        observed_version = run_checked([conda, "run", "-p", str(prefix), "dotmatch", "--version"], cwd=root, env=env)
        expected_version = f"dotmatch {version}"
        if observed_version != expected_version:
            raise RuntimeError(f"Bioconda dotmatch --version reported {observed_version!r}, expected {expected_version!r}")
        observed_distance = run_checked([conda, "run", "-p", str(prefix), "dotmatch", "dist", "ACGT", "AGGT"], cwd=root, env=env)
        if observed_distance != "1":
            raise RuntimeError(f"Bioconda dotmatch dist smoke test reported {observed_distance!r}, expected '1'")
        observed_threshold = run_checked([conda, "run", "-p", str(prefix), "dotmatch", "leq", "1", "ACGT", "AGGT"], cwd=root, env=env)
        if observed_threshold != "true":
            raise RuntimeError(f"Bioconda dotmatch leq smoke test reported {observed_threshold!r}, expected 'true'")


def verify_biocontainers_run(image: str, version: str) -> None:
    env = os.environ.copy()
    cwd = Path.cwd()
    observed_distance = run_checked(["docker", "run", "--rm", image, "dotmatch", "dist", "ACGT", "AGGT"], cwd=cwd, env=env)
    if observed_distance != "1":
        raise RuntimeError(f"BioContainers dotmatch dist smoke test reported {observed_distance!r}, expected '1'")
    observed_threshold = run_checked(["docker", "run", "--rm", image, "dotmatch", "leq", "1", "ACGT", "AGGT"], cwd=cwd, env=env)
    if observed_threshold != "true":
        raise RuntimeError(f"BioContainers dotmatch leq smoke test reported {observed_threshold!r}, expected 'true'")


def check_pypi(version: str, result: AuditResult) -> None:
    channel = "pypi"
    try:
        data = fetch_json(PYPI_URL.format(version=version))
    except Exception as exc:
        result.failures.append(ChannelMessage(channel, f"PyPI version {version} is not reachable: {exc}"))
        return
    urls = data.get("urls") or []
    has_sdist = any(item.get("packagetype") == "sdist" for item in urls if isinstance(item, dict))
    wheels = [item for item in urls if isinstance(item, dict) and item.get("packagetype") == "bdist_wheel"]
    has_macos_wheel = any("macosx_" in str(item.get("filename") or "") for item in wheels)
    has_manylinux_wheel = any("manylinux" in str(item.get("filename") or "") for item in wheels)
    has_musllinux_wheel = any("musllinux" in str(item.get("filename") or "") for item in wheels)
    has_raw_linux_wheel = any(
        "linux_x86_64" in str(item.get("filename") or "")
        and "manylinux" not in str(item.get("filename") or "")
        and "musllinux" not in str(item.get("filename") or "")
        for item in wheels
    )
    if data.get("info", {}).get("version") != version or not has_sdist:
        result.failures.append(ChannelMessage(channel, f"PyPI version {version} is not available as an sdist"))
        return
    if not has_macos_wheel:
        result.failures.append(ChannelMessage(channel, f"PyPI version {version} must include a macOS wheel"))
        return
    if not has_manylinux_wheel or not has_musllinux_wheel:
        result.failures.append(
            ChannelMessage(channel, f"PyPI version {version} must include repaired manylinux and musllinux wheels")
        )
        return
    if has_raw_linux_wheel:
        result.failures.append(ChannelMessage(channel, f"PyPI version {version} must not include raw linux_x86_64 wheels"))
        return
    result.passed.append(
        ChannelMessage(channel, f"PyPI sdist, macOS wheel, and repaired Linux wheels are available for {version}")
    )
    try:
        verify_pypi_install(version)
    except Exception as exc:
        result.failures.append(ChannelMessage("pypi-install", f"PyPI one-command install failed for {version}: {exc}"))
        return
    result.passed.append(ChannelMessage("pypi-install", f"pip install dotmatch=={version} works in a clean environment"))


def check_bioconda(version: str, result: AuditResult) -> None:
    channel = "bioconda"
    try:
        data = fetch_json(BIOCONDA_URL)
    except Exception as exc:
        result.failures.append(ChannelMessage(channel, f"Bioconda package metadata is not reachable: {exc}"))
        return
    files = data.get("files") or []
    if not any(item.get("version") == version for item in files if isinstance(item, dict)):
        result.failures.append(ChannelMessage(channel, f"Bioconda version {version} is not available"))
        return
    result.passed.append(ChannelMessage(channel, f"Bioconda package is available for {version}"))
    try:
        verify_bioconda_install(version)
    except Exception as exc:
        result.failures.append(ChannelMessage("bioconda-install", f"Bioconda one-command install failed for {version}: {exc}"))
        return
    result.passed.append(ChannelMessage("bioconda-install", f"Bioconda install and CLI smoke tests pass for {version}"))


def biocontainers_tags_for_version(version: str) -> list[str]:
    tags: list[str] = []
    page = 1
    while page <= 20:
        data = fetch_json(BIOCONTAINERS_TAGS_URL.format(page=page))
        for item in data.get("tags") or []:
            if not isinstance(item, dict):
                continue
            name = str(item.get("name") or "")
            if name.startswith(f"{version}--"):
                tags.append(name)
        if not data.get("has_additional"):
            break
        page += 1
    return sorted(tags)


def check_biocontainers(version: str, result: AuditResult) -> None:
    channel = "biocontainers"
    try:
        tags = biocontainers_tags_for_version(version)
    except Exception as exc:
        result.failures.append(ChannelMessage(channel, f"BioContainers tags are not reachable: {exc}"))
        return
    if not tags:
        result.failures.append(ChannelMessage(channel, f"BioContainers image tag for version {version} is not available"))
        return
    tag = tags[0]
    image = BIOCONTAINERS_IMAGE.format(tag=tag)
    try:
        proc = subprocess.run(
            ["docker", "manifest", "inspect", image],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        result.failures.append(ChannelMessage(channel, "docker is required to verify the BioContainers image manifest"))
        return
    if proc.returncode != 0:
        detail = (proc.stderr or "").strip() or (proc.stdout or "").strip()
        suffix = f": {detail}" if detail else ""
        result.failures.append(ChannelMessage(channel, f"BioContainers image tag {image} is not available{suffix}"))
        return
    result.passed.append(ChannelMessage(channel, f"BioContainers image tag is available: {image}"))
    try:
        verify_biocontainers_run(image, version)
    except Exception as exc:
        result.failures.append(ChannelMessage("biocontainers-run", f"BioContainers image runtime smoke test failed for {image}: {exc}"))
        return
    result.passed.append(ChannelMessage("biocontainers-run", f"BioContainers docker run smoke tests pass for {image}"))


def check_ghcr(version: str, result: AuditResult) -> None:
    channel = "ghcr"
    image = GHCR_IMAGE.format(version=version)
    try:
        proc = subprocess.run(
            ["docker", "manifest", "inspect", image],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        result.failures.append(ChannelMessage(channel, "docker is required to verify the GHCR image manifest"))
        return
    if proc.returncode != 0:
        detail = (proc.stderr or "").strip() or (proc.stdout or "").strip()
        suffix = f": {detail}" if detail else ""
        result.failures.append(ChannelMessage(channel, f"GHCR image tag {image} is not available{suffix}"))
        return
    result.passed.append(ChannelMessage(channel, f"GHCR image tag is available: {image}"))
    try:
        verify_ghcr_run(image, version)
    except Exception as exc:
        result.failures.append(ChannelMessage("ghcr-run", f"GHCR image runtime smoke test failed for {image}: {exc}"))
        return
    result.passed.append(ChannelMessage("ghcr-run", f"docker run smoke tests pass for {image}"))


def check_zenodo(root: Path, result: AuditResult) -> None:
    channel = "zenodo"
    doi = citation_doi(root)
    if not doi:
        result.failures.append(ChannelMessage(channel, "CITATION.cff must include a DOI after Zenodo release"))
        return
    url = f"https://doi.org/{doi}"
    if not url_ok(url):
        result.failures.append(ChannelMessage(channel, f"Zenodo DOI does not resolve: {doi}"))
        return
    result.passed.append(ChannelMessage(channel, f"Zenodo DOI resolves: {doi}"))


def audit(root: Path, version: Optional[str] = None) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    try:
        release_version = version or project_version(root)
    except Exception as exc:
        result.failures.append(ChannelMessage("metadata", str(exc)))
        return result
    check_pypi(release_version, result)
    check_bioconda(release_version, result)
    check_biocontainers(release_version, result)
    check_ghcr(release_version, result)
    check_zenodo(root, result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="repository root")
    parser.add_argument("--version", default="", help="release version; defaults to pyproject.toml")
    args = parser.parse_args()

    result = audit(Path(args.root), args.version or None)
    for item in result.passed:
        print(f"PASS [{item.channel}]: {item.message}")
    for item in result.failures:
        print(f"FAIL [{item.channel}]: {item.message}")
    if result.ok:
        print("DISTRIBUTION CHANNELS: PASS")
        return 0
    print("DISTRIBUTION CHANNELS: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
