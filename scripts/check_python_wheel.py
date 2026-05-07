#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import platform
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


def wheel_native_members(wheel: Path) -> list[str]:
    with zipfile.ZipFile(wheel) as archive:
        return [
            name
            for name in archive.namelist()
            if name.startswith("dotmatch/") and (name.endswith(".so") or name.endswith(".dylib"))
        ]


def check_sdist_members(sdist: Path) -> None:
    required_suffixes = [
        "/src/qdalign.c",
        "/include/qdalign.h",
        "/setup.py",
        "/pyproject.toml",
        "/README.md",
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


def venv_python(env_dir: Path) -> Path:
    if os.name == "nt":
        return env_dir / "Scripts" / "python.exe"
    return env_dir / "bin" / "python"


def clean_import_env() -> dict[str, str]:
    env = os.environ.copy()
    env.pop("DOTMATCH_LIB", None)
    env.pop("QUICKDNA_LIB", None)
    env.pop("PYTHONPATH", None)
    return env


def verify_clean_install(artifact: Path, install_root: Path) -> None:
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
    run([str(py), "-c", probe], cwd=probe_dir, env=clean_import_env())


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


def main() -> int:
    parser = argparse.ArgumentParser(description="Build and verify the DotMatch Python wheel.")
    parser.add_argument("--out-dir", default="", help="optional wheel output directory")
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
            sdist_dir = install_root / "sdist"
            sdist_dir.mkdir()
            run([sys.executable, "-m", "build", "--sdist", "--outdir", str(sdist_dir)], cwd=ROOT)
            sdists = sorted(sdist_dir.glob("dotmatch-*.tar.gz"))
            if len(sdists) != 1:
                raise SystemExit(f"expected exactly one dotmatch sdist in {sdist_dir}, found {len(sdists)}")
            sdist = sdists[0]
            check_sdist_members(sdist)

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
            check_macos_tag(wheel)
            check_macos_architecture(wheel, native_members[0])

            verify_clean_install(wheel, install_root / "wheel-install")
            verify_clean_install(sdist, install_root / "sdist-install")
            print(f"verified {wheel.name} with native payload: {', '.join(native_members)}")
            print(f"verified {sdist.name} source build")
        finally:
            if cleanup_out:
                shutil.rmtree(out_dir, ignore_errors=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
