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
GHCR_TOKEN_URL = "https://ghcr.io/token?service=ghcr.io&scope=repository:{repository}:pull"
BIOCONTAINERS_TAGS_URL = (
    "https://quay.io/api/v1/repository/biocontainers/dotmatch/tag/?onlyActiveTags=true&page={page}&limit=100"
)
BIOCONTAINERS_IMAGE = "quay.io/biocontainers/dotmatch:{tag}"
ZENODO_RECORD_URL = "https://zenodo.org/api/records/{record_id}"


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


def fetch_registry_manifest(image: str) -> tuple[dict, str]:
    registry_repo, reference = image.rsplit(":", 1)
    registry, repository = registry_repo.split("/", 1)
    headers = {
        "Accept": (
            "application/vnd.oci.image.index.v1+json, "
            "application/vnd.docker.distribution.manifest.list.v2+json, "
            "application/vnd.oci.image.manifest.v1+json, "
            "application/vnd.docker.distribution.manifest.v2+json"
        ),
        "User-Agent": "DotMatch distribution verifier",
    }
    if registry == "ghcr.io":
        token_data = fetch_json(GHCR_TOKEN_URL.format(repository=repository))
        token = str(token_data.get("token") or "")
        if not token:
            raise RuntimeError(f"GHCR did not return a pull token for {repository}")
        headers["Authorization"] = f"Bearer {token}"

    request = urllib.request.Request(f"https://{registry}/v2/{repository}/manifests/{reference}", headers=headers)
    with urllib.request.urlopen(request, timeout=30) as response:
        digest = str(response.headers.get("Docker-Content-Digest") or "")
        data = json.loads(response.read().decode("utf-8"))
    return data, digest


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


def verify_ghcr_manifest(image: str) -> str:
    data, digest = fetch_registry_manifest(image)
    if int(data.get("schemaVersion") or 0) != 2:
        raise RuntimeError("GHCR manifest must use schemaVersion 2")
    manifests = data.get("manifests") or []
    if manifests:
        linux_amd64 = [
            item
            for item in manifests
            if isinstance(item, dict)
            and isinstance(item.get("platform"), dict)
            and item["platform"].get("os") == "linux"
            and item["platform"].get("architecture") == "amd64"
        ]
        if not linux_amd64:
            raise RuntimeError("GHCR manifest list must include linux/amd64")
    if not digest:
        digest = str(data.get("config", {}).get("digest") or "")
    if not digest.startswith("sha256:"):
        raise RuntimeError("GHCR manifest did not include a sha256 digest")
    return digest


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
        (root / "gc_library.tsv").write_text("guide\tbases\tgene\ng0\tACGT\tGENE0\n", encoding="utf-8")
        (root / "gc_reads.fastq").write_text("@r0\nACGT\n+\nIIII\n", encoding="utf-8")
        run_checked(
            [
                conda,
                "run",
                "-p",
                str(prefix),
                "dotmatch",
                "guide-counter",
                "count",
                "--input",
                "gc_reads.fastq",
                "--samples",
                "sample",
                "--library",
                "gc_library.tsv",
                "--offset-sample-size",
                "1",
                "--offset-min-fraction",
                "0.1",
                "--output",
                "gc_out",
            ],
            cwd=root,
            env=env,
        )
        counts = (root / "gc_out.counts.txt").read_text(encoding="utf-8").splitlines()
        extended = (root / "gc_out.extended-counts.txt").read_text(encoding="utf-8").splitlines()
        stats = (root / "gc_out.stats.txt").read_text(encoding="utf-8").splitlines()
        if counts != ["guide\tgene\tsample", "g0\tGENE0\t1"]:
            raise RuntimeError("Bioconda guide-counter counts smoke test produced unexpected counts")
        if not extended or extended[0] != "guide\tgene\tguide_type\tsample":
            raise RuntimeError("Bioconda guide-counter extended-counts smoke test produced unexpected header")
        if len(stats) < 2 or not stats[1].startswith("gc_reads.fastq\tsample\t1\t1\t1\t1.0000"):
            raise RuntimeError("Bioconda guide-counter stats smoke test produced unexpected stats")


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
    result.passed.append(
        ChannelMessage("bioconda-install", f"Bioconda install, CLI, and GuideCounter-compatible smoke tests pass for {version}")
    )


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
        digest = verify_ghcr_manifest(image)
    except Exception as exc:
        result.failures.append(ChannelMessage(channel, f"GHCR image tag {image} is not available: {exc}"))
        return
    result.passed.append(ChannelMessage(channel, f"GHCR image tag is available: {image} ({digest})"))
    try:
        verify_ghcr_run(image, version)
    except FileNotFoundError:
        result.failures.append(ChannelMessage("ghcr-run", "docker is required to run GHCR image smoke tests"))
        return
    except Exception as exc:
        result.failures.append(ChannelMessage("ghcr-run", f"GHCR image runtime smoke test failed for {image}: {exc}"))
        return
    result.passed.append(ChannelMessage("ghcr-run", f"docker run smoke tests pass for {image}"))


def check_zenodo(root: Path, version: str, result: AuditResult) -> None:
    channel = "zenodo"
    doi = citation_doi(root)
    if not doi:
        result.failures.append(ChannelMessage(channel, "CITATION.cff must include a DOI after Zenodo release"))
        return
    record_id = doi.rsplit(".", 1)[-1] if doi.startswith("10.5281/zenodo.") else ""
    if not record_id:
        result.failures.append(ChannelMessage(channel, f"Zenodo DOI is not a Zenodo record DOI: {doi}"))
        return
    try:
        data = fetch_json(ZENODO_RECORD_URL.format(record_id=record_id))
    except Exception as exc:
        result.failures.append(ChannelMessage(channel, f"Zenodo record metadata is not reachable for {doi}: {exc}"))
        return
    metadata = data.get("metadata") if isinstance(data, dict) else {}
    record_version = str((metadata or {}).get("version") or "")
    if record_version != version:
        result.failures.append(
            ChannelMessage(channel, f"Zenodo record {doi} reports version {record_version or '<missing>'}, expected {version}")
        )
        return
    url = f"https://doi.org/{doi}"
    if not url_ok(url):
        result.failures.append(ChannelMessage(channel, f"Zenodo DOI does not resolve: {doi}"))
        return
    result.passed.append(ChannelMessage(channel, f"Zenodo DOI resolves and reports version {version}: {doi}"))


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
    check_zenodo(root, release_version, result)
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
