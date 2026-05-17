from __future__ import annotations

import os
import platform
import re
import shlex
import subprocess
import sysconfig
from pathlib import Path

from setuptools import Distribution, setup
from setuptools.command.build_py import build_py as _build_py

try:
    from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel
except ImportError:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel


ROOT = Path(__file__).resolve().parent


def _project_version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    if match is None:
        raise RuntimeError("pyproject.toml version not found")
    return match.group(1)


def _macos_arch_flags() -> list[str]:
    flags = shlex.split(sysconfig.get_config_var("CFLAGS") or "")
    arch_flags: list[str] = []
    i = 0
    while i < len(flags):
        if flags[i] == "-arch" and i + 1 < len(flags):
            arch_flags.extend(flags[i : i + 2])
            i += 2
            continue
        i += 1
    return arch_flags


def _macos_deployment_target(arch_flags: list[str]) -> str:
    configured = sysconfig.get_config_var("MACOSX_DEPLOYMENT_TARGET") or "10.9"
    major_text, _sep, minor_text = configured.partition(".")
    major = int(major_text)
    minor = int(minor_text or 0)
    if "arm64" in arch_flags and major < 11:
        return "11.0"
    return f"{major}.{minor}"


if platform.system() == "Darwin":
    _initial_arch_flags = _macos_arch_flags()
    os.environ.setdefault("MACOSX_DEPLOYMENT_TARGET", _macos_deployment_target(_initial_arch_flags))


class BinaryDistribution(Distribution):
    def has_ext_modules(self) -> bool:
        return True


class build_py(_build_py):
    def run(self) -> None:
        super().run()
        self.build_native_library()
        self.build_native_cli()

    def build_native_library(self) -> None:
        system = platform.system()
        arch_flags: list[str] = []
        if system == "Darwin":
            lib_name = "libdotmatch.dylib"
            link_flags = ["-dynamiclib", "-install_name", "@rpath/libdotmatch.dylib"]
            arch_flags = _macos_arch_flags()
            build_env = os.environ.copy()
            build_env.setdefault("MACOSX_DEPLOYMENT_TARGET", _macos_deployment_target(arch_flags))
        elif system == "Linux":
            lib_name = "libdotmatch.so"
            link_flags = ["-shared"]
            build_env = os.environ.copy()
        else:
            raise RuntimeError(f"DotMatch Python wheels are not yet supported on {system}")

        build_cmd = self.get_finalized_command("build")
        build_temp = Path(build_cmd.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        object_path = build_temp / "qdalign.package.o"
        out_dir = Path(self.build_lib) / "dotmatch"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / lib_name

        cc = shlex.split(os.environ.get("CC") or sysconfig.get_config_var("CC") or "cc")
        cflags = shlex.split(os.environ.get("CFLAGS", ""))
        compile_cmd = [
            *cc,
            *arch_flags,
            "-O3",
            "-std=c11",
            "-Wall",
            "-Wextra",
            "-Wpedantic",
            "-fPIC",
            "-I",
            str(ROOT / "include"),
            *cflags,
            "-c",
            str(ROOT / "src" / "qdalign.c"),
            "-o",
            str(object_path),
        ]
        link_cmd = [*cc, *arch_flags, *link_flags, str(object_path), "-o", str(out_path)]
        subprocess.run(compile_cmd, check=True, env=build_env)
        subprocess.run(link_cmd, check=True, env=build_env)

    def build_native_cli(self) -> None:
        system = platform.system()
        arch_flags: list[str] = []
        if system == "Darwin":
            arch_flags = _macos_arch_flags()
            build_env = os.environ.copy()
            build_env.setdefault("MACOSX_DEPLOYMENT_TARGET", _macos_deployment_target(arch_flags))
        elif system == "Linux":
            build_env = os.environ.copy()
        else:
            raise RuntimeError(f"DotMatch Python wheels are not yet supported on {system}")

        build_cmd = self.get_finalized_command("build")
        build_temp = Path(build_cmd.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        out_dir = Path(self.build_lib) / "dotmatch"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "dotmatch-native"

        cc = shlex.split(os.environ.get("CC") or sysconfig.get_config_var("CC") or "cc")
        cflags = shlex.split(os.environ.get("CFLAGS", ""))
        compile_cmd = [
            *cc,
            *arch_flags,
            "-O3",
            "-std=c11",
            "-Wall",
            "-Wextra",
            "-Wpedantic",
            "-I",
            str(ROOT / "include"),
            f'-DDOTMATCH_VERSION="{_project_version()}"',
            *cflags,
            str(ROOT / "src/qda.c"),
            str(ROOT / "src/qdalign.c"),
            "-o",
            str(out_path),
            "-lz",
            "-pthread",
        ]
        subprocess.run(compile_cmd, check=True, env=build_env)
        out_path.chmod(0o755)


class bdist_wheel(_bdist_wheel):
    def finalize_options(self) -> None:
        super().finalize_options()
        self.root_is_pure = False

    def get_tag(self) -> tuple[str, str, str]:
        _python, _abi, platform_tag = super().get_tag()
        return "py3", "none", platform_tag


setup(
    cmdclass={"build_py": build_py, "bdist_wheel": bdist_wheel},
    distclass=BinaryDistribution,
    package_data={"dotmatch": ["libdotmatch.*", "dotmatch-native"]},
)
