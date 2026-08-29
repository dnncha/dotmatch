#!/usr/bin/env python3
"""Run DotMatch workflow-manager integration smoke tests.

This runner keeps the expensive ecosystem checks in one place so local release
work and CI exercise the same artifacts: nf-test modules, the small Nextflow
pipeline, Snakemake, Galaxy wrapper linting and CRISPR-count execution, and
MultiQC custom/plugin reports.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NFCORE_MODULE_ROOT = ROOT / "examples" / "workflows" / "nf-core" / "modules" / "local" / "dotmatch"
NFCORE_MODULES = ["count", "demux", "audit", "panel_check", "crispr_count", "assay_run"]
GALAXY_WRAPPERS = [
    ROOT / "examples" / "workflows" / "galaxy" / "dotmatch_crispr_count.xml",
    ROOT / "examples" / "workflows" / "galaxy" / "dotmatch_demux.xml",
    ROOT / "examples" / "workflows" / "galaxy" / "dotmatch_panel_check.xml",
]
GALAXY_CRISPR_WRAPPER = GALAXY_WRAPPERS[0]


class Failure(Exception):
    """Workflow integration test failure."""


def _tool(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise Failure(f"required workflow test tool is missing from PATH: {name}")
    return path


def _run(cmd: list[str], *, cwd: Path, env: dict[str, str]) -> None:
    print(f"+ ({cwd.relative_to(ROOT) if cwd.is_relative_to(ROOT) else cwd}) {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, cwd=cwd, env=env, check=True)


def _dotmatch_env(tmp: Path) -> dict[str, str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = f"{ROOT / 'python'}{os.pathsep}{env.get('PYTHONPATH', '')}".rstrip(os.pathsep)
    native_cli = ROOT / "dotmatch"
    if native_cli.exists():
        env.setdefault("DOTMATCH_NATIVE_CLI", str(native_cli))

    bin_dir = tmp / "bin"
    bin_dir.mkdir()
    wrapper = bin_dir / "dotmatch"
    python = sys.executable
    wrapper.write_text(
        "#!/usr/bin/env sh\n"
        f"PYTHONPATH='{ROOT / 'python'}' "
        f"DOTMATCH_NATIVE_CLI='{native_cli}' "
        f"exec '{python}' -m dotmatch.cli \"$@\"\n",
        encoding="utf-8",
    )
    wrapper.chmod(0o755)
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    return env


def run_nf_tests(env: dict[str, str]) -> None:
    _tool("nf-test")
    for module in NFCORE_MODULES:
        module_dir = NFCORE_MODULE_ROOT / module
        _run(["nf-test", "test", "tests/main.nf.test"], cwd=module_dir, env=env)


def run_nextflow_pipeline(tmp: Path, env: dict[str, str]) -> None:
    _tool("nextflow")
    outdir = tmp / "nextflow-results"
    workdir = tmp / "nextflow-work"
    _run(
        [
            "nextflow",
            "run",
            "main.nf",
            "--outdir",
            str(outdir),
            "-work-dir",
            str(workdir),
        ],
        cwd=ROOT / "examples" / "workflows" / "nf-core" / "pipeline",
        env=env,
    )


def run_snakemake(tmp: Path, env: dict[str, str]) -> None:
    _tool("snakemake")
    fixtures = ROOT / "examples" / "workflows" / "fixtures"
    with (ROOT / "examples" / "workflows" / "snakemake" / "config.json").open(
        encoding="utf-8"
    ) as handle:
        config = json.load(handle)
    config.update(
        {
            "library": str(fixtures / "crispr_library.csv"),
            "samples": {
                "sample_a": str(fixtures / "sample_a.fastq"),
                "sample_b": str(fixtures / "sample_b.fastq"),
            },
            "guide_start": 0,
            "guide_length": 4,
            "metric": "hamming",
            "indel_window": 0,
            "outdir": str(tmp / "snakemake-output"),
        }
    )
    config_path = tmp / "snakemake-config.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    _run(
        [
            "snakemake",
            "-s",
            str(ROOT / "examples" / "workflows" / "snakemake" / "Snakefile"),
            "--configfile",
            str(config_path),
            "--cores",
            "1",
            "-p",
        ],
        cwd=ROOT,
        env=env,
    )


def run_planemo(env: dict[str, str]) -> None:
    _tool("planemo")
    _run(["planemo", "lint", *[str(path) for path in GALAXY_WRAPPERS]], cwd=ROOT, env=env)
    _run(
        ["planemo", "test", "--install_galaxy", str(GALAXY_CRISPR_WRAPPER)],
        cwd=ROOT,
        env=env,
    )


def run_multiqc(tmp: Path, env: dict[str, str]) -> None:
    _tool("multiqc")
    custom_out = tmp / "multiqc-custom"
    plugin_out = tmp / "multiqc-plugin"
    _run(
        [
            "multiqc",
            str(ROOT / "examples" / "workflows" / "multiqc" / "data"),
            "-c",
            str(ROOT / "examples" / "workflows" / "multiqc" / "multiqc_config.yaml"),
            "-o",
            str(custom_out),
            "--filename",
            "dotmatch_multiqc.html",
            "--force",
        ],
        cwd=ROOT,
        env=env,
    )
    _run(
        [
            "multiqc",
            str(ROOT / "examples" / "workflows" / "multiqc" / "data"),
            "--module",
            "dotmatch",
            "-o",
            str(plugin_out),
            "--filename",
            "dotmatch_plugin_multiqc.html",
            "--force",
        ],
        cwd=ROOT,
        env=env,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-nf-test", action="store_true", help="skip nf-test module checks")
    parser.add_argument("--skip-nextflow", action="store_true", help="skip the Nextflow pipeline check")
    parser.add_argument("--skip-snakemake", action="store_true", help="skip the Snakemake workflow check")
    parser.add_argument(
        "--skip-planemo",
        action="store_true",
        help="skip Galaxy wrapper linting and CRISPR-count execution",
    )
    parser.add_argument("--skip-multiqc", action="store_true", help="skip MultiQC report checks")
    args = parser.parse_args()

    try:
        with tempfile.TemporaryDirectory(prefix="dotmatch-workflow-tests-") as tmp_name:
            tmp = Path(tmp_name)
            env = _dotmatch_env(tmp)
            if not args.skip_nf_test:
                run_nf_tests(env)
            if not args.skip_nextflow:
                run_nextflow_pipeline(tmp, env)
            if not args.skip_snakemake:
                run_snakemake(tmp, env)
            if not args.skip_planemo:
                run_planemo(env)
            if not args.skip_multiqc:
                run_multiqc(tmp, env)
    except subprocess.CalledProcessError as exc:
        print(f"WORKFLOW INTEGRATION: FAIL ({exc})", file=sys.stderr)
        return exc.returncode or 1
    except Failure as exc:
        print(f"WORKFLOW INTEGRATION: FAIL ({exc})", file=sys.stderr)
        return 1

    print("WORKFLOW INTEGRATION: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
