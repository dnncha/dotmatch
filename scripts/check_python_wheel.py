#!/usr/bin/env python3
from __future__ import annotations

import argparse
from email.parser import Parser
from email.message import Message
import hashlib
import json
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
REPAIRED_LINUX_WHEEL_ARCHITECTURES = ("x86_64", "aarch64")
HOST_ARCHITECTURE_ALIASES = {
    "amd64": "x86_64",
    "arm64": "aarch64",
    "x86_64": "x86_64",
    "aarch64": "aarch64",
}


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


def run_expect_exit(
    cmd: list[str],
    expected_returncode: int,
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode != expected_returncode:
        raise SystemExit(
            f"{cmd!r} returned {result.returncode}, expected {expected_returncode}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
    return result


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


def wheel_assay_evidence_members(wheel: Path) -> list[str]:
    with zipfile.ZipFile(wheel) as archive:
        return [name for name in archive.namelist() if name == "dotmatch/data/assay-evidence.json"]


def wheel_agent_discovery_members(wheel: Path) -> list[str]:
    required = {
        "dotmatch/data/agent-capabilities.json",
        "dotmatch/data/agent-capabilities.schema.json",
        "dotmatch/data/agent-tools.json",
        "dotmatch/data/agent-tools.schema.json",
        "dotmatch/data/codex-skill/SKILL.md",
        "dotmatch/data/codex-skill/agents/openai.yaml",
    }
    with zipfile.ZipFile(wheel) as archive:
        return sorted(required.intersection(archive.namelist()))


def check_sdist_members(sdist: Path) -> None:
    required_suffixes = [
        "/src/qdalign.c",
        "/src/qdmetal_stub.c",
        "/src/qda.c",
        "/include/qdalign.h",
        "/include/qdmetal.h",
        "/setup.py",
        "/pyproject.toml",
        "/README.md",
        "/CITATION.cff",
        "/codemeta.json",
        "/agent-capabilities.json",
        "/agent-capabilities.schema.json",
        "/agent-capabilities.v1.json",
        "/agent-capabilities.v1.schema.json",
        "/agent-tools.json",
        "/agent-tools.schema.json",
        "/agent-reference-crispr.json",
        "/extensions/codex/dotmatch-agent/SKILL.md",
        "/docs/assay-evidence.json",
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
    if "known-target short-dna assignment" not in str(metadata.get("Summary", "")).lower():
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
    for label in ["Homepage", "Repository", "Issues", "Documentation", "Agent guide"]:
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

    capability_text = run_text(
        [str(venv_script(env_dir, "dotmatch")), "capabilities", "--json"],
        cwd=probe_dir,
        env=env,
    )
    capability_data = json.loads(capability_text)
    if capability_data.get("generated_for_version") != expected_version:
        raise SystemExit(f"{artifact.name} installed capability version does not match {expected_version}")
    capability_ids = {
        item.get("id")
        for item in capability_data.get("intents", [])
        if isinstance(item, dict)
    }
    for required_id in [
        "crispr-guide-counting",
        "inline-barcode-demultiplexing",
        "feature-barcode-assignment",
        "perturb-seq-guide-capture",
        "barcode-panel-design",
        "known-target-fastq-matching",
    ]:
        if required_id not in capability_ids:
            raise SystemExit(f"{artifact.name} installed capabilities are missing {required_id}")

    tool_contract_text = run_text(
        [str(venv_script(env_dir, "dotmatch")), "agent", "tools", "--json"],
        cwd=probe_dir,
        env=env,
    )
    tool_contract = json.loads(tool_contract_text)
    tool_names = {item.get("name") for item in tool_contract.get("tools", []) if isinstance(item, dict)}
    if tool_names != {"discover", "prepare_assay", "inspect_assay", "run_assay", "review_assay", "handoff_assay"}:
        raise SystemExit(f"{artifact.name} installed agent tool contract is incomplete: {tool_names!r}")
    exported_skill = probe_dir / "exported-dotmatch-skill"
    export_text = run_text(
        [str(venv_script(env_dir, "dotmatch")), "agent", "export-skill", "--target", str(exported_skill)],
        cwd=probe_dir,
        env=env,
    )
    export_result = json.loads(export_text)
    if export_result.get("status") != "passed" or not (exported_skill / "SKILL.md").is_file():
        raise SystemExit(f"{artifact.name} could not export the installed Codex skill")

    dist_observed = run_text([str(venv_script(env_dir, "dotmatch")), "dist", "ACGT", "AGGT"], cwd=probe_dir, env=env)
    if dist_observed != "1":
        raise SystemExit(f"{artifact.name} console CLI distance smoke test returned {dist_observed!r}")

    feature_targets = probe_dir / "feature_targets.tsv"
    feature_observations = probe_dir / "feature_observations.tsv"
    feature_output = probe_dir / "feature_matrix"
    feature_targets.write_text("target_id\ttarget_seq\nfeature_a\tACGT\nfeature_b\tTTTT\n", encoding="utf-8")
    feature_observations.write_text(
        "observation_id\tcell_barcode\tfeature_seq\n"
        "feature_read_1\tcell_a\tACGT\n"
        "feature_read_2\tcell_b\tTTTT\n"
        "feature_read_3\tcell_b\tCCCC\n",
        encoding="utf-8",
    )
    run(
        [
            str(venv_script(env_dir, "dotmatch")),
            "feature",
            "matrix",
            "--observations",
            str(feature_observations),
            "--targets",
            str(feature_targets),
            "--id-column",
            "observation_id",
            "--cell-column",
            "cell_barcode",
            "--sequence-column",
            "feature_seq",
            "--metric",
            "hamming",
            "--k",
            "0",
            "--out-dir",
            str(feature_output),
        ],
        cwd=probe_dir,
        env=env,
    )
    feature_summary = json.loads((feature_output / "summary.json").read_text(encoding="utf-8"))
    if feature_summary.get("assigned_unique") != 2 or feature_summary.get("unmatched") != 1:
        raise SystemExit(f"{artifact.name} installed feature matrix summary is invalid: {feature_summary!r}")
    if not (feature_output / "matrix.mtx").is_file():
        raise SystemExit(f"{artifact.name} installed feature matrix command did not write matrix.mtx")

    left_targets = probe_dir / "left_targets.tsv"
    right_targets = probe_dir / "right_targets.tsv"
    left_reads = probe_dir / "pair_R1.fastq"
    right_reads = probe_dir / "pair_R2.fastq"
    pair_counts = probe_dir / "pair_counts.tsv"
    left_targets.write_text("left_a\tACGT\n", encoding="utf-8")
    right_targets.write_text("right_a\tGGAA\n", encoding="utf-8")
    left_reads.write_text("@pair_1/1\nACGT\n+\nIIII\n", encoding="utf-8")
    right_reads.write_text("@pair_1/2\nGGAA\n+\nIIII\n", encoding="utf-8")
    run(
        [
            str(venv_script(env_dir, "dotmatch")),
            "pair-count",
            "--left-targets",
            str(left_targets),
            "--right-targets",
            str(right_targets),
            "--left-reads",
            str(left_reads),
            "--right-reads",
            str(right_reads),
            "--left-start",
            "0",
            "--left-length",
            "4",
            "--right-start",
            "0",
            "--right-length",
            "4",
            "--k",
            "0",
            "--metric",
            "hamming",
            "--out",
            str(pair_counts),
        ],
        cwd=probe_dir,
        env=env,
    )
    if "left_a\tright_a\t1" not in pair_counts.read_text(encoding="utf-8"):
        raise SystemExit(f"{artifact.name} installed paired FASTQ pair-count output is invalid")

    targets = probe_dir / "targets.tsv"
    reads = probe_dir / "reads.fastq"
    spec = probe_dir / "assay.toml"
    targets.write_text(
        "g1\tACGT\tG1\n"
        "g2\tAGGT\tG2\n"
        "g3\tATGT\tG3\n"
        "g4\tACGG\tG4\n",
        encoding="utf-8",
    )
    reads.write_text(
        "@r0\nACGT\n+\nIIII\n"
        "@r1\nAGGT\n+\nIIII\n"
        "@r2\nATGT\n+\nIIII\n"
        "@r3\nACGG\n+\nIIII\n",
        encoding="utf-8",
    )

    agent_smoke_counts = probe_dir / "agent_smoke_counts.tsv"
    agent_smoke_summary = probe_dir / "agent_smoke_summary.json"
    run(
        [
            str(venv_script(env_dir, "dotmatch")),
            "count",
            "--targets",
            str(targets),
            "--reads",
            str(reads),
            "--sample-label",
            "agent-smoke",
            "--target-start",
            "0",
            "--target-length",
            "4",
            "--k",
            "0",
            "--metric",
            "hamming",
            "--out",
            str(agent_smoke_counts),
            "--summary",
            str(agent_smoke_summary),
        ],
        cwd=probe_dir,
        env=env,
    )
    agent_summary = json.loads(agent_smoke_summary.read_text(encoding="utf-8"))
    samples = agent_summary.get("samples", [])
    if len(samples) != 1 or samples[0].get("total_reads") != 4 or samples[0].get("assigned_unique") != 4:
        raise SystemExit(f"{artifact.name} clean-install FASTQ count returned unexpected summary: {samples!r}")
    if "g1\tACGT\tG1" not in agent_smoke_counts.read_text(encoding="utf-8"):
        raise SystemExit(f"{artifact.name} clean-install FASTQ count did not write the expected target")
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
k = 0
metric = "hamming"
""".lstrip(),
        encoding="utf-8",
    )
    run([str(venv_script(env_dir, "dotmatch")), "assay", "check", str(spec)], cwd=probe_dir, env=env)
    inspect_input = probe_dir / "inspect-agent-input.json"
    inspect_input.write_text(json.dumps({"spec": str(spec)}) + "\n", encoding="utf-8")
    inspect_result = run_expect_exit(
        [str(venv_script(env_dir, "dotmatch")), "agent", "invoke", "inspect_assay", "--input", str(inspect_input)],
        0,
        cwd=probe_dir,
        env=env,
    )
    inspect_envelope = json.loads(inspect_result.stdout)
    if inspect_envelope.get("status") != "passed" or inspect_envelope.get("spec", {}).get("sha256", "") == "":
        raise SystemExit(f"{artifact.name} installed agent inspect returned invalid evidence")
    reliability_summary = json.loads((probe_dir / "assay_out" / "reliability_summary.json").read_text(encoding="utf-8"))
    evidence_boundary = reliability_summary.get("evidence_boundary", {})
    if evidence_boundary.get("status") != "supported" or evidence_boundary.get("id") != "crispr_guide_counting":
        raise SystemExit(f"{artifact.name} installed assay reliability evidence boundary is invalid: {evidence_boundary!r}")
    run_input = probe_dir / "run-agent-input.json"
    run_input.write_text(json.dumps({"spec": str(spec), "minimum_free_bytes": 0}) + "\n", encoding="utf-8")
    run_result = run_expect_exit(
        [str(venv_script(env_dir, "dotmatch")), "agent", "invoke", "run_assay", "--input", str(run_input)],
        0,
        cwd=probe_dir,
        env=env,
    )
    run_envelope = json.loads(run_result.stdout)
    if run_envelope.get("status") != "passed" or not run_envelope.get("result", {}).get("revision_history"):
        raise SystemExit(f"{artifact.name} installed agent run returned invalid evidence")
    reliability_dir = probe_dir / "assay_out"
    postrun_summary = json.loads((reliability_dir / "reliability_summary.json").read_text(encoding="utf-8"))
    postrun_evidence = postrun_summary.get("evidence_boundary", {})
    if postrun_summary.get("stage") != "postrun" or postrun_summary.get("overall_status") not in {"passed", "informational"}:
        raise SystemExit(f"{artifact.name} installed assay run reliability summary is invalid: {postrun_summary!r}")
    if postrun_evidence.get("status") != "supported" or postrun_evidence.get("id") != "crispr_guide_counting":
        raise SystemExit(f"{artifact.name} installed assay run evidence boundary is invalid: {postrun_evidence!r}")
    review_input = probe_dir / "review-agent-input.json"
    review_input.write_text(json.dumps({"spec": str(spec)}) + "\n", encoding="utf-8")
    review_result = run_expect_exit(
        [str(venv_script(env_dir, "dotmatch")), "agent", "invoke", "review_assay", "--input", str(review_input)],
        0,
        cwd=probe_dir,
        env=env,
    )
    review_envelope = json.loads(review_result.stdout)
    if review_envelope.get("status") != "passed" or not review_envelope.get("artifacts"):
        raise SystemExit(f"{artifact.name} installed agent review returned invalid evidence")
    handoff_dir = probe_dir / "agent-handoff"
    handoff_input = probe_dir / "handoff-agent-input.json"
    handoff_input.write_text(
        json.dumps({"spec": str(spec), "output_dir": str(handoff_dir), "minimum_free_bytes": 0}) + "\n",
        encoding="utf-8",
    )
    handoff_result = run_expect_exit(
        [str(venv_script(env_dir, "dotmatch")), "agent", "invoke", "handoff_assay", "--input", str(handoff_input)],
        0,
        cwd=probe_dir,
        env=env,
    )
    handoff_envelope = json.loads(handoff_result.stdout)
    if handoff_envelope.get("status") != "passed" or handoff_envelope.get("result", {}).get("raw_data_included") is not False:
        raise SystemExit(f"{artifact.name} installed agent handoff returned invalid evidence")
    if any(path.suffix in {".fastq", ".fq"} or path.name.endswith((".fastq.gz", ".fq.gz")) for path in handoff_dir.rglob("*")):
        raise SystemExit(f"{artifact.name} installed agent handoff copied raw FASTQ")
    for required in [
        "reliability_summary.json",
        "reliability_findings.tsv",
        "reliability_report.html",
        "reliability_manifest.summary.tsv",
    ]:
        if not (reliability_dir / required).exists():
            raise SystemExit(f"{artifact.name} installed assay run did not write {required}")
    handoff_dir = probe_dir / "handoff"
    run(
        [
            str(venv_script(env_dir, "dotmatch")),
            "assay",
            "handoff",
            str(spec),
            "--out-dir",
            str(handoff_dir),
        ],
        cwd=probe_dir,
        env=env,
    )
    handoff = json.loads((handoff_dir / "handoff_manifest.json").read_text(encoding="utf-8"))
    if handoff.get("bundle_type") != "dotmatch_assay_handoff":
        raise SystemExit(f"{artifact.name} installed assay handoff manifest is invalid: {handoff!r}")
    input_roles = {str(item.get("role")) for item in handoff.get("inputs", []) if isinstance(item, dict)}
    if input_roles != {"targets", "sample:sample"}:
        raise SystemExit(f"{artifact.name} installed assay handoff inputs are invalid: {input_roles!r}")
    review_files = [item for item in handoff.get("review_files", []) if isinstance(item, dict)]
    review_paths = {str(item.get("path")) for item in review_files}
    if "review/reliability_report.html" not in review_paths or "review/counts.mageck.tsv" not in review_paths:
        raise SystemExit(f"{artifact.name} installed assay handoff is missing required review files: {review_paths!r}")
    for item in review_files:
        review_path = handoff_dir / str(item["path"])
        observed_hash = hashlib.sha256(review_path.read_bytes()).hexdigest()
        if observed_hash != item.get("sha256"):
            raise SystemExit(f"{artifact.name} installed assay handoff checksum is invalid for {review_path}")
    counts = probe_dir / "counts.mageck.tsv"
    sample_qc = probe_dir / "sample_qc.tsv"
    crispr_qc = probe_dir / "crispr_qc.json"
    counts.write_text("sgRNA\tGene\tsample\nguide_a\tGENEA\t1\n", encoding="utf-8")
    sample_qc.write_text(
        "sample_id\tassignment_rate\tambiguous_rate\tno_match_rate\tinvalid_rate\n"
        "sample\t1\t0\t0\t0\n",
        encoding="utf-8",
    )
    run(
        [
            str(venv_script(env_dir, "dotmatch")),
            "crispr",
            "qc",
            "--counts",
            str(counts),
            "--sample-qc",
            str(sample_qc),
            "--library",
            str(targets),
            "--out",
            str(crispr_qc),
        ],
        cwd=probe_dir,
        env=env,
    )
    run(
        [
            str(venv_script(env_dir, "dotmatch")),
            "crispr-qc",
            "--counts",
            str(counts),
            "--sample-qc",
            str(sample_qc),
            "--library",
            str(targets),
            "--out",
            str(probe_dir / "crispr_qc_legacy.json"),
        ],
        cwd=probe_dir,
        env=env,
    )
    run(
        [
            str(venv_script(env_dir, "dotmatch")),
            "crispr",
            "infer",
            "--library",
            str(targets),
            "--reads",
            str(reads),
            "--out",
            str(probe_dir / "crispr_inferred.toml"),
            "--report",
            str(probe_dir / "crispr_inference_report.json"),
        ],
        cwd=probe_dir,
        env=env,
    )
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
    check_result = run_expect_exit(
        [str(venv_script(env_dir, "dotmatch")), "assay", "check", str(inferred)],
        2,
        cwd=probe_dir,
        env=env,
    )
    if "draft_assayspec" not in f"{check_result.stdout}\n{check_result.stderr}":
        raise SystemExit("inferred draft assay check did not report the draft_assayspec blocker")
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


def wheel_platform_tags(wheel: Path) -> list[str]:
    if not wheel.name.endswith(".whl"):
        return []
    fields = wheel.name[:-4].rsplit("-", 3)
    if len(fields) != 4:
        return []
    return fields[-1].split(".")


def repaired_linux_wheel_architectures(wheel: Path, family: str) -> set[str]:
    return {
        architecture
        for architecture in REPAIRED_LINUX_WHEEL_ARCHITECTURES
        if any(
            tag.startswith(family) and tag.endswith(f"_{architecture}")
            for tag in wheel_platform_tags(wheel)
        )
    }


def require_repaired_linux_wheel_architectures(wheels: list[Path], required_architectures: list[str]) -> None:
    missing = [
        f"{family}_{architecture}"
        for family in ["manylinux", "musllinux"]
        for architecture in required_architectures
        if not any(architecture in repaired_linux_wheel_architectures(wheel, family) for wheel in wheels)
    ]
    if missing:
        raise SystemExit("missing repaired Linux wheel coverage: " + ", ".join(missing))


def wheel_supported_by_current_platform(wheel: Path) -> bool:
    name = wheel.name
    system = platform.system()
    libc_name = platform.libc_ver()[0].lower()
    if "musllinux" in name:
        architectures = repaired_linux_wheel_architectures(wheel, "musllinux")
        host_architecture = HOST_ARCHITECTURE_ALIASES.get(platform.machine().lower(), platform.machine().lower())
        return system == "Linux" and libc_name == "musl" and host_architecture in architectures
    if "manylinux" in name:
        architectures = repaired_linux_wheel_architectures(wheel, "manylinux")
        host_architecture = HOST_ARCHITECTURE_ALIASES.get(platform.machine().lower(), platform.machine().lower())
        return system == "Linux" and libc_name == "glibc" and host_architecture in architectures
    if "macosx" in name:
        return system == "Darwin"
    if "win_" in name or "win32" in name or "win_amd64" in name:
        return os.name == "nt"
    return True


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
    assay_evidence_members = wheel_assay_evidence_members(wheel)
    if not assay_evidence_members:
        raise SystemExit(f"{wheel.name} does not contain dotmatch/data/assay-evidence.json")
    agent_discovery_members = wheel_agent_discovery_members(wheel)
    if len(agent_discovery_members) != 6:
        raise SystemExit(f"{wheel.name} does not contain all installed agent contract and skill files")
    check_distribution_metadata(wheel, expected_version)
    check_macos_tag(wheel)
    check_macos_architecture(wheel, native_members[0])
    verify_clean_install(wheel, install_root, expected_version)
    return wheel, native_members + native_cli_members + assay_evidence_members + agent_discovery_members


def verify_existing_wheels(wheel_dir: Path, install_root: Path, expected_version: str) -> list[Path]:
    wheels = sorted(wheel_dir.glob("dotmatch-*.whl"))
    if not wheels:
        raise SystemExit(f"expected at least one dotmatch wheel in {wheel_dir}")
    clean_install_index = 0
    for wheel in wheels:
        native_members = wheel_native_members(wheel)
        if not native_members:
            raise SystemExit(f"{wheel.name} does not contain dotmatch/libdotmatch.*")
        if not wheel_native_cli_members(wheel):
            raise SystemExit(f"{wheel.name} does not contain dotmatch-native")
        if not wheel_assay_evidence_members(wheel):
            raise SystemExit(f"{wheel.name} does not contain dotmatch/data/assay-evidence.json")
        if len(wheel_agent_discovery_members(wheel)) != 6:
            raise SystemExit(f"{wheel.name} does not contain all installed agent contract and skill files")
        check_distribution_metadata(wheel, expected_version)
        if not wheel_supported_by_current_platform(wheel):
            print(f"skipping clean install for unsupported wheel tag on this host: {wheel.name}")
            continue
        verify_clean_install(wheel, install_root / f"wheel-{clean_install_index}", expected_version)
        clean_install_index += 1
    return wheels


def main() -> int:
    parser = argparse.ArgumentParser(description="Build and verify the DotMatch Python wheel.")
    parser.add_argument("--out-dir", default="", help="optional wheel output directory")
    parser.add_argument("--sdist-only", action="store_true", help="build and verify only the source distribution")
    parser.add_argument("--wheel-only", action="store_true", help="verify existing wheels in --out-dir without building")
    parser.add_argument(
        "--require-repaired-linux-architectures",
        choices=REPAIRED_LINUX_WHEEL_ARCHITECTURES,
        metavar="ARCH",
        nargs="+",
        help="require repaired manylinux and musllinux wheels for each architecture",
    )
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
            if args.wheel_only:
                wheels = verify_existing_wheels(out_dir, install_root / "existing-wheel-install", expected_version)
                if args.require_repaired_linux_architectures:
                    require_repaired_linux_wheel_architectures(wheels, args.require_repaired_linux_architectures)
                print("verified existing wheels: " + ", ".join(wheel.name for wheel in wheels))
                return 0
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
