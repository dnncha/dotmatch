#!/usr/bin/env python3

import argparse
import json
import re
from pathlib import Path
from typing import Optional


class ReleaseAudit:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _json(path: Path) -> dict[str, object]:
    return json.loads(_read(path))


def _pyproject_version(path: Path) -> Optional[str]:
    in_project = False
    for raw_line in _read(path).splitlines():
        line = raw_line.strip()
        if line.startswith("[") and line.endswith("]"):
            in_project = line == "[project]"
            continue
        if not in_project or not line.startswith("version"):
            continue
        match = re.match(r'version\s*=\s*"([^"]+)"', line)
        if match:
            return match.group(1)
    return None


def _cff_version(path: Path) -> Optional[str]:
    for raw_line in _read(path).splitlines():
        match = re.match(r'\s*version\s*:\s*["\']?([^"\']+)["\']?\s*$', raw_line)
        if match:
            return match.group(1).strip()
    return None


def _docker_label_version(text: str) -> Optional[str]:
    match = re.search(r'org\.opencontainers\.image\.version="([^"]+)"', text)
    return match.group(1) if match else None


def _bioconda_template_version(text: str) -> Optional[str]:
    match = re.search(r'{%\s*set\s+version\s*=\s*"([^"]+)"\s*%}', text)
    return match.group(1) if match else None


def _contains_doi_claim(path: Path, *, allow_codemeta_context: bool = False) -> bool:
    if path.name == "codemeta.json" and allow_codemeta_context:
        data = _json(path)
        data.pop("@context", None)
        return "doi.org" in json.dumps(data).lower() or re.search(r'\bdoi\s*[:=]', json.dumps(data), re.I) is not None
    text = _read(path)
    return "doi.org" in text.lower() or re.search(r'^\s*doi\s*:', text, re.I | re.M) is not None


def _workflow_job_block(workflow: str, job_name: str) -> str:
    match = re.search(rf"^  {re.escape(job_name)}:\n", workflow, flags=re.M)
    if not match:
        return ""
    start = match.start()
    next_job = re.search(r"^  [A-Za-z0-9_-]+:\n", workflow[match.end():], flags=re.M)
    end = match.end() + next_job.start() if next_job else len(workflow)
    return workflow[start:end]


def _make_target_block(makefile: str, target_name: str) -> str:
    match = re.search(rf"^{re.escape(target_name)}:\n", makefile, flags=re.M)
    if not match:
        return ""
    start = match.start()
    next_target = re.search(r"^[A-Za-z0-9_.-]+:(?:\s|$)", makefile[match.end():], flags=re.M)
    end = match.end() + next_target.start() if next_target else len(makefile)
    return makefile[start:end]


def check_versions(root: Path, result: ReleaseAudit) -> None:
    version_files: dict[str, Optional[str]] = {}
    try:
        version_files["pyproject.toml"] = _pyproject_version(root / "pyproject.toml")
        version_files["package.json"] = str(_json(root / "package.json").get("version") or "")
        codemeta = _json(root / "codemeta.json")
        version_files["codemeta.json version"] = str(codemeta.get("version") or "")
        version_files["codemeta.json softwareVersion"] = str(codemeta.get("softwareVersion") or "")
        version_files[".zenodo.json"] = str(_json(root / ".zenodo.json").get("version") or "")
        version_files["CITATION.cff"] = _cff_version(root / "CITATION.cff")
        version_files["Dockerfile OCI label"] = _docker_label_version(_read(root / "Dockerfile"))
        version_files["packaging/bioconda/meta.yaml"] = _bioconda_template_version(
            _read(root / "packaging" / "bioconda" / "meta.yaml")
        )
    except Exception as exc:
        result.failures.append(f"release version metadata could not be read: {exc}")
        return

    missing = [name for name, version in version_files.items() if not version]
    result.failures.extend(f"{name} must declare release version" for name in missing)

    declared = {name: version for name, version in version_files.items() if version}
    unique_versions = sorted(set(declared.values()))
    if len(unique_versions) > 1:
        detail = ", ".join(f"{name}={version}" for name, version in sorted(declared.items()))
        result.failures.append(f"release version mismatch: {detail}")

    if not missing and len(unique_versions) == 1:
        result.passed.append("release versions aligned")


def check_no_unminted_doi_claims(root: Path, result: ReleaseAudit) -> None:
    checked = [
        root / "CITATION.cff",
        root / ".zenodo.json",
        root / "codemeta.json",
    ]
    for path in checked:
        if _contains_doi_claim(path, allow_codemeta_context=True):
            result.failures.append(f"{path.name} must not claim a DOI before an immutable release DOI is minted")
    if not any("DOI" in failure for failure in result.failures):
        result.passed.append("DOI claims deferred until minted release")


def check_sdist_metadata(root: Path, result: ReleaseAudit) -> None:
    manifest = _read(root / "MANIFEST.in")
    verifier = _read(root / "scripts" / "check_python_wheel.py")
    for required in ["CITATION.cff", "codemeta.json", "src/qdalign.c", "include/qdalign.h"]:
        if f"include {required}" not in manifest:
            result.failures.append(f"MANIFEST.in must include {required}")
    for required_suffix in ["/CITATION.cff", "/codemeta.json", "/src/qdalign.c", "/include/qdalign.h"]:
        if required_suffix not in verifier:
            result.failures.append(f"scripts/check_python_wheel.py must verify {required_suffix}")
    if not any("MANIFEST.in" in failure or "check_python_wheel.py" in failure for failure in result.failures):
        result.passed.append("sdist release metadata verified")


def check_distribution_surfaces(root: Path, result: ReleaseAudit) -> None:
    workflow = _read(root / ".github" / "workflows" / "release.yml")
    dockerfile = _read(root / "Dockerfile")
    bioconda = _read(root / "packaging" / "bioconda" / "meta.yaml")
    packaging = _read(root / "docs" / "packaging.md")
    release_process = _read(root / "docs" / "release-process.md")
    makefile = _read(root / "Makefile")

    required_workflow_fragments = [
        "id-token: write",
        "packages: write",
        "pypa/gh-action-pypi-publish@release/v1",
        "packages-dir: dist-pypi",
        "docker/metadata-action",
        "docker/build-push-action",
        "ghcr.io/dnncha/dotmatch",
        "docker image inspect dotmatch:ci",
        "SHA256SUMS.txt",
    ]
    for fragment in required_workflow_fragments:
        if fragment not in workflow:
            result.failures.append(f"release workflow missing {fragment}")
    if "dotmatch-wheel-Linux" in workflow:
        result.failures.append("release workflow must not publish raw Linux wheels to PyPI")
    preflight = _workflow_job_block(workflow, "preflight")
    container_job = _workflow_job_block(workflow, "container")
    sdist_job = _workflow_job_block(workflow, "sdist")
    pypi_job = _workflow_job_block(workflow, "pypi-sdist")
    github_release_job = _workflow_job_block(workflow, "github-release")
    if not preflight:
        result.failures.append("release workflow missing preflight job")
    else:
        if "Release preflight gates" not in preflight:
            result.failures.append("release workflow preflight job must be named Release preflight gates")
        if "python -m pip install build pytest" not in preflight:
            result.failures.append("release workflow preflight job must install pytest")
        if "make test" not in preflight:
            result.failures.append("release workflow preflight job must run make test")
        if "make cli-test" not in preflight:
            result.failures.append("release workflow preflight job must run make cli-test")
        if "make python-test" not in preflight:
            result.failures.append("release workflow preflight job must run make python-test")
        if "make repository-ready" not in preflight:
            result.failures.append("release workflow preflight job must run make repository-ready")
        if "make release-ready" not in preflight:
            result.failures.append("release workflow preflight job must run make release-ready")
        if "make python-package-test" not in preflight:
            result.failures.append("release workflow preflight job must run make python-package-test")
    if "needs: [preflight]" not in container_job:
        result.failures.append("container publish job must depend on preflight")
    if "python scripts/check_python_wheel.py --sdist-only --out-dir dist" not in sdist_job:
        result.failures.append("release workflow sdist job must verify the PyPI source distribution artifact")
    if "Publish PyPI sdist, macOS wheel, and repaired Linux wheels" not in pypi_job:
        result.failures.append("PyPI publish job must publish sdist, macOS wheel, and repaired Linux wheels")
    if "needs: [preflight, sdist, wheel, linux-repaired-wheels]" not in pypi_job:
        result.failures.append("PyPI publish job must depend on preflight, sdist, macOS wheel, and repaired Linux wheels")
    if (
        "name: dotmatch-sdist" not in pypi_job
        or "name: dotmatch-wheel-macos" not in pypi_job
        or "name: dotmatch-linux-repaired-wheels" not in pypi_job
    ):
        result.failures.append("PyPI publish job must download sdist, macOS wheel, and repaired Linux wheel artifacts")
    if "needs: [preflight, wheel, sdist, linux-repaired-wheels]" not in github_release_job:
        result.failures.append("GitHub release job must depend on preflight, wheels, sdist, and repaired Linux wheels")

    required_labels = [
        "org.opencontainers.image.title",
        "org.opencontainers.image.source",
        "org.opencontainers.image.url",
        "org.opencontainers.image.version",
        "org.opencontainers.image.licenses",
    ]
    for label in required_labels:
        if label not in dockerfile:
            result.failures.append(f"Dockerfile missing OCI label {label}")

    if "REPLACE_WITH_RELEASE_TARBALL_SHA256" not in bioconda:
        result.failures.append("Bioconda template must retain release SHA256 placeholder until copying into bioconda-recipes")
    if "dotmatch dist ACGT AGGT" not in bioconda:
        result.failures.append("Bioconda template must include native CLI smoke test")

    if (
        "publishes the source distribution, the native macOS wheel, and repaired manylinux/musllinux Linux wheels"
        not in packaging
    ):
        result.failures.append("docs/packaging.md must document PyPI sdist, macOS wheel, and repaired Linux wheel policy")
    if "Raw `linux_x86_64` wheels" not in packaging:
        result.failures.append("docs/packaging.md must document raw Linux wheels are not uploaded to PyPI")
    if "ghcr.io/dnncha/dotmatch" not in packaging or "OCI labels" not in packaging:
        result.failures.append("docs/packaging.md must document container registry and OCI labels")
    if "quay.io/biocontainers/dotmatch" not in packaging:
        result.failures.append("docs/packaging.md must document BioContainers post-release verification")
    if "make bioconda-recipe-ready" not in packaging:
        result.failures.append("docs/packaging.md must document make bioconda-recipe-ready")
    if "docs/distribution-release.json" not in packaging:
        result.failures.append("docs/packaging.md must document docs/distribution-release.json")
    if "make release-ready" not in release_process:
        result.failures.append("docs/release-process.md must include make release-ready")
    if "make pretag-ready" not in release_process:
        result.failures.append("docs/release-process.md must include make pretag-ready")
    pretag_block = _make_target_block(makefile, "pretag-ready")
    if not pretag_block:
        result.failures.append("Makefile must include pretag-ready target")
    for pretag_fragment in [
        "$(MAKE) test",
        "$(MAKE) cli-test",
        "$(MAKE) python-test",
        "$(MAKE) python-package-test",
        "$(MAKE) repository-ready",
        "$(MAKE) release-ready",
        "$(MAKE) coverage",
        "npm run lint",
        "npm audit --audit-level=moderate",
        "NEXT_OUTPUT=export NEXT_PUBLIC_BASE_PATH=/dotmatch NEXT_PUBLIC_SITE_URL=https://dnncha.github.io/dotmatch npm run build",
    ]:
        if pretag_fragment not in pretag_block:
            result.failures.append(f"Makefile pretag-ready target must include {pretag_fragment}")
    for post_release_gate in [
        "distribution-channels",
        "workflow-adoption-status",
        "bcl-comparison-gate",
    ]:
        if post_release_gate in pretag_block:
            result.failures.append(f"Makefile pretag-ready target must not include {post_release_gate}")
        if f"make {post_release_gate}" not in release_process:
            result.failures.append(f"docs/release-process.md must document separate {post_release_gate} gate")
    if "make assay-evidence-ready" not in release_process:
        result.failures.append("docs/release-process.md must include make assay-evidence-ready")
    if "make distribution-record-ready" not in release_process:
        result.failures.append("docs/release-process.md must include make distribution-record-ready")
    if "repaired manylinux/musllinux wheels" not in release_process:
        result.failures.append("docs/release-process.md must document repaired manylinux/musllinux PyPI wheel publishing")
    if "make bioconda-recipe-ready" not in release_process:
        result.failures.append("docs/release-process.md must include make bioconda-recipe-ready")
    if "make alphabet-policy-ready" not in release_process:
        result.failures.append("docs/release-process.md must include make alphabet-policy-ready")
    if "make citation-metadata-ready" not in release_process:
        result.failures.append("docs/release-process.md must include make citation-metadata-ready")
    if "make native-comparator-scope-ready" not in release_process:
        result.failures.append("docs/release-process.md must include make native-comparator-scope-ready")
    for evidence_gate in [
        "make public-crispr-evidence-gate",
        "make crispr-comparison-gate",
        "make barcode-comparison-gate",
        "make feature-barcode-public-gate",
        "make perturb-seq-public-gate",
        "make amplicon-panel-public-gate",
        "make bcl-tiny-public-gate",
        "make oligo-adapter-public-gate",
        "make workflow-examples-ready",
    ]:
        if evidence_gate not in release_process:
            result.failures.append(f"docs/release-process.md must include {evidence_gate}")

    if not any(
        marker in failure
        for failure in result.failures
        for marker in [
            "release workflow",
            "Dockerfile",
            "Bioconda",
            "docs/packaging.md",
            "release-process",
        ]
    ):
        result.passed.append("distribution surfaces release-ready")


def audit(root: Path) -> ReleaseAudit:
    root = root.resolve()
    result = ReleaseAudit()
    check_versions(root, result)
    check_no_unminted_doi_claims(root, result)
    check_sdist_metadata(root, result)
    check_distribution_surfaces(root, result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit DotMatch release-readiness metadata and distribution surfaces.")
    parser.add_argument("--root", default=".", help="repository root")
    args = parser.parse_args()

    result = audit(Path(args.root))
    for item in result.passed:
        print(f"PASS: {item}")
    for item in result.failures:
        print(f"FAIL: {item}")
    if result.ok:
        print("RELEASE READINESS: PASS")
        return 0
    print("RELEASE READINESS: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
