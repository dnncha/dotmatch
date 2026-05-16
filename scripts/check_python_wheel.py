#!/usr/bin/env python3
from __future__ import annotations

import argparse
from email.parser import Parser
from email.message import Message
import os
import platform
import re
import shutil
import subprocess
import sys
import sysconfig
import tarfile
import tempfile
import venv
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(cmd: list[str], *, cwd: Path | None = None, env: dict[str, str] | None = None) -> None:
    subprocess.run(cmd, cwd=cwd, env=env, check=True)


def run_text(cmd: list[str], *, cwd: Path | None = None, env: dict[str, str] | None = None) -> str:
    result = subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return result.stdout.strip()


def project_version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    if match is None:
        raise SystemExit("pyproject.toml does not declare project version")
    return match.group(1)


def wheel_native_members(wheel: Path) -> list[str]:
    with zipfile.ZipFile(wheel) as archive:
        return [
            name
            for name in archive.namelist()
            if name.startswith("dotmatch/") and (name.endswith(".so") or name.endswith(".dylib"))
        ]


def wheel_native_cli_members(wheel: Path) -> list[str]:
    with zipfile.ZipFile(wheel) as archive:
        return [name for name in archive.namelist() if name == "dotmatch/dotmatch-native"]


def check_sdist_members(sdist: Path) -> None:
    required_suffixes = [
        "/src/qdalign.c",
        "/src/qda.c",
        "/include/qdalign.h",
        "/setup.py",
        "/pyproject.toml",
        "/README.md",
        "/CITATION.cff",
        "/codemeta.json",
        "/LICENSE",
    ]
    with tarfile.open(sdist, "r:gz") as archive:
        names = archive.getnames()
    missing = [
        suffix
        for suffix in required_suffixes
        if not any(name.endswith(suffix) for name in names)
    ]
    if missing:
        raise SystemExit(f"{sdist.name} is missing required source files: {', '.join(missing)}")


def read_distribution_metadata(artifact: Path) -> Message:
    if artifact.suffix == ".whl":
        with zipfile.ZipFile(artifact) as archive:
            metadata_members = [
                name
                for name in archive.namelist()
                if name.endswith(".dist-info/METADATA")
            ]
            if len(metadata_members) != 1:
                raise SystemExit(
                    f"{artifact.name} must contain exactly one .dist-info/METADATA file, found {len(metadata_members)}"
                )
            text = archive.read(metadata_members[0]).decode("utf-8")
    elif artifact.name.endswith(".tar.gz"):
        with tarfile.open(artifact, "r:gz") as archive:
            all_metadata_members = [
                member
                for member in archive.getmembers()
                if member.name.endswith("/PKG-INFO") or member.name == "PKG-INFO"
            ]
            metadata_members = [
                member
                for member in all_metadata_members
                if len(Path(member.name).parts) == 2 and Path(member.name).name == "PKG-INFO"
            ]
            if len(metadata_members) != 1:
                raise SystemExit(
                    f"{artifact.name} must contain exactly one top-level PKG-INFO file, found {len(metadata_members)}"
                )
            extracted = archive.extractfile(metadata_members[0])
            if extracted is None:
                raise SystemExit(f"{artifact.name} PKG-INFO could not be read")
            text = extracted.read().decode("utf-8")
    else:
        raise SystemExit(f"{artifact.name} is not a supported Python distribution artifact")
    return Parser().parsestr(text)


def _metadata_values(metadata: Message, key: str) -> list[str]:
    return [str(value) for value in metadata.get_all(key, [])]


def _metadata_contains(values: list[str], fragment: str) -> bool:
    return fragment.lower() in "\n".join(values).lower()


def _project_url_labels(metadata: Message) -> set[str]:
    labels: set[str] = set()
    for value in _metadata_values(metadata, "Project-URL"):
        label, _sep, _url = value.partition(",")
        labels.add(label.strip())
    return labels


def check_distribution_metadata(artifact: Path, expected_version: str) -> None:
    metadata = read_distribution_metadata(artifact)
    failures: list[str] = []

    if metadata.get("Name") != "dotmatch":
        failures.append("Name must be dotmatch")
    if metadata.get("Version") != expected_version:
        failures.append(f"Version must be {expected_version}")
    if "known-target short-DNA assignment" not in str(metadata.get("Summary", "")):
        failures.append("Summary must mention known-target short-DNA assignment")

    license_text = "\n".join(_metadata_values(metadata, "License-Expression") + _metadata_values(metadata, "License"))
    if "Apache-2.0" not in license_text:
        failures.append("License metadata must include Apache-2.0")

    keywords = _metadata_values(metadata, "Keywords")
    for keyword in ["bioinformatics", "computational biology", "CRISPR", "FASTQ", "known-target assignment"]:
        if not _metadata_contains(keywords, keyword):
            failures.append(f"Keywords must include {keyword}")

    classifiers = _metadata_values(metadata, "Classifier")
    for classifier in [
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ]:
        if classifier not in classifiers:
            failures.append(f"Classifier must include {classifier}")

    project_urls = _project_url_labels(metadata)
    for label in ["Homepage", "Repository", "Issues", "Documentation"]:
        if label not in project_urls:
            failures.append(f"Project-URL must include {label}")

    if failures:
        raise SystemExit(f"{artifact.name} has invalid PyPI metadata: {'; '.join(failures)}")


def venv_python(env_dir: Path) -> Path:
    if os.name == "nt":
        return env_dir / "Scripts" / "python.exe"
    return env_dir / "bin" / "python"


def venv_script(env_dir: Path, name: str) -> Path:
    if os.name == "nt":
        return env_dir / "Scripts" / f"{name}.exe"
    return env_dir / "bin" / name


def clean_import_env() -> dict[str, str]:
    env = os.environ.copy()
    env.pop("DOTMATCH_LIB", None)
    env.pop("QUICKDNA_LIB", None)
    env.pop("PYTHONPATH", None)
    return env


def verify_clean_install(artifact: Path, install_root: Path, expected_version: str) -> None:
    env_dir = install_root / "venv"
    venv.EnvBuilder(with_pip=True).create(env_dir)
    py = venv_python(env_dir)
    run([str(py), "-m", "pip", "install", "--quiet", str(artifact)])

    probe_dir = install_root / "probe"
    probe_dir.mkdir()
    probe = (
        "import dotmatch, quickdna; "
        "assert dotmatch.distance('ACGT', 'AGGT') == 1; "
        "assert quickdna.distance_leq('ACGT', 'AGGT', 1); "
        "print('dotmatch package import ok')"
    )
    env = clean_import_env()
    run([str(py), "-c", probe], cwd=probe_dir, env=env)
    for label, cmd in [
        ("module CLI", [str(py), "-m", "dotmatch.cli", "--version"]),
        ("console CLI", [str(venv_script(env_dir, "dotmatch")), "--version"]),
    ]:
        observed = run_text(cmd, cwd=probe_dir, env=env)
        expected = f"dotmatch {expected_version}"
        if observed != expected:
            raise SystemExit(f"{artifact.name} {label} reported {observed!r}, expected {expected!r}")

    dist_observed = run_text([str(venv_script(env_dir, "dotmatch")), "dist", "ACGT", "AGGT"], cwd=probe_dir, env=env)
    if dist_observed != "1":
        raise SystemExit(f"{artifact.name} console CLI distance smoke test returned {dist_observed!r}")

    targets = probe_dir / "targets.tsv"
    reads = probe_dir / "reads.fastq"
    spec = probe_dir / "assay.toml"
    targets.write_text("guide_a\tACGT\tGENEA\n", encoding="utf-8")
    reads.write_text("@r0\nACGT\n+\nIIII\n", encoding="utf-8")
    spec.write_text(
        f"""
schema_version = 1
mode = "count"
assay_type = "crispr"
targets = "{targets}"

[[samples]]
id = "sample"
fastq = "{reads}"

[run]
out_dir = "{probe_dir / 'assay_out'}"

[extract]
start = 0
length = 4

[assignment]
k = 1
metric = "hamming"
""".lstrip(),
        encoding="utf-8",
    )
    run([str(venv_script(env_dir, "dotmatch")), "assay", "check", str(spec)], cwd=probe_dir, env=env)
    inferred = probe_dir / "inferred.toml"
    inferred_report = probe_dir / "inference_report.json"
    run(
        [
            str(venv_script(env_dir, "dotmatch")),
            "assay",
            "infer",
            "--mode",
            "count",
            "--assay-type",
            "crispr",
            "--targets",
            str(targets),
            "--reads",
            str(reads),
            "--sample-id",
            "sample",
            "--out",
            str(inferred),
            "--report",
            str(inferred_report),
        ],
        cwd=probe_dir,
        env=env,
    )
    run([str(venv_script(env_dir, "dotmatch")), "assay", "check", str(inferred)], cwd=probe_dir, env=env)
    run(
        [
            str(venv_script(env_dir, "dotmatch")),
            "assay",
            "autopsy",
            str(spec),
            "--out-dir",
            str(probe_dir / "autopsy"),
        ],
        cwd=probe_dir,
        env=env,
    )


def check_macos_tag(wheel: Path) -> None:
    if platform.system() != "Darwin":
        return
    configured = sysconfig.get_config_var("MACOSX_DEPLOYMENT_TARGET") or "10.9"
    major, _sep, minor = configured.partition(".")
    expected_major = int(major)
    expected_minor = int(minor or 0)
    if "universal2" in wheel.name or platform.machine() == "arm64":
        expected_major = max(expected_major, 11)
        if expected_major == 11:
            expected_minor = 0
    expected_prefix = f"macosx_{expected_major}_{expected_minor}"
    if expected_prefix not in wheel.name:
        raise SystemExit(
            f"{wheel.name} does not use the interpreter deployment target prefix {expected_prefix}"
        )


def check_macos_architecture(wheel: Path, native_member: str) -> None:
    if platform.system() != "Darwin":
        return
    with tempfile.TemporaryDirectory(prefix="dotmatch-wheel-arch-") as tmp:
        tmp_path = Path(tmp)
        with zipfile.ZipFile(wheel) as archive:
            archive.extract(native_member, tmp_path)
        native_path = tmp_path / native_member
        result = subprocess.run(
            ["file", str(native_path)],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        description = result.stdout
        if "universal2" in wheel.name and ("arm64" not in description or "x86_64" not in description):
            raise SystemExit(
                f"{wheel.name} is tagged universal2 but {native_member} is not universal: {description.strip()}"
            )


def build_and_verify_sdist(out_dir: Path, install_root: Path, expected_version: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    run([sys.executable, "-m", "build", "--sdist", "--outdir", str(out_dir)], cwd=ROOT)
    sdists = sorted(out_dir.glob("dotmatch-*.tar.gz"))
    if len(sdists) != 1:
        raise SystemExit(f"expected exactly one dotmatch sdist in {out_dir}, found {len(sdists)}")
    sdist = sdists[0]
    check_sdist_members(sdist)
    check_distribution_metadata(sdist, expected_version)
    verify_clean_install(sdist, install_root, expected_version)
    return sdist


def build_and_verify_wheel(out_dir: Path, install_root: Path, expected_version: str) -> tuple[Path, list[str]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    run([sys.executable, "-m", "build", "--wheel", "--outdir", str(out_dir)], cwd=ROOT)
    wheels = sorted(out_dir.glob("dotmatch-*.whl"))
    if len(wheels) != 1:
        raise SystemExit(f"expected exactly one dotmatch wheel in {out_dir}, found {len(wheels)}")
    wheel = wheels[0]
    if "-py3-none-" not in wheel.name:
        raise SystemExit(f"{wheel.name} should use a py3-none platform tag")
    native_members = wheel_native_members(wheel)
    if not native_members:
        raise SystemExit(f"{wheel.name} does not contain dotmatch/libdotmatch.*")
    native_cli_members = wheel_native_cli_members(wheel)
    if not native_cli_members:
        raise SystemExit(f"{wheel.name} does not contain dotmatch-native")
    check_distribution_metadata(wheel, expected_version)
    check_macos_tag(wheel)
    check_macos_architecture(wheel, native_members[0])
    verify_clean_install(wheel, install_root, expected_version)
    return wheel, native_members + native_cli_members


def main() -> int:
    parser = argparse.ArgumentParser(description="Build and verify the DotMatch Python wheel.")
    parser.add_argument("--out-dir", default="", help="optional wheel output directory")
    parser.add_argument("--sdist-only", action="store_true", help="build and verify only the source distribution")
    args = parser.parse_args()

    if args.out_dir:
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        cleanup_out = False
    else:
        out_dir = Path(tempfile.mkdtemp(prefix="dotmatch-wheel-"))
        cleanup_out = True

    with tempfile.TemporaryDirectory(prefix="dotmatch-wheel-install-") as install_tmp:
        install_root = Path(install_tmp)
        try:
            expected_version = project_version()
            sdist_out_dir = out_dir if args.sdist_only else install_root / "sdist"
            sdist = build_and_verify_sdist(sdist_out_dir, install_root / "sdist-install", expected_version)
            if args.sdist_only:
                print(f"verified {sdist.name} source build")
                return 0

            wheel, native_members = build_and_verify_wheel(out_dir, install_root / "wheel-install", expected_version)
            print(f"verified {wheel.name} with native payload: {', '.join(native_members)}")
            print(f"verified {sdist.name} source build")
        finally:
            if cleanup_out:
                shutil.rmtree(out_dir, ignore_errors=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
