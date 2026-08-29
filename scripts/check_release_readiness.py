#!/usr/bin/env python3

import argparse
import json
import re
import urllib.request
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


def _c_header_version(text: str) -> Optional[str]:
    match = re.search(r'#define\s+QDALN_VERSION\s+"([^"]+)"', text)
    return match.group(1) if match else None


def _has_release_doi_field(path: Path) -> bool:
    if path.name == "CITATION.cff":
        return re.search(r"^\s*doi\s*:", _read(path), re.I | re.M) is not None
    if path.suffix == ".json":
        data = _json(path)
        return any(key in data for key in {"doi", "conceptdoi"})
    return False


def _doi_values(path: Path) -> list[str]:
    if path.name == "CITATION.cff":
        match = re.search(r'^\s*doi\s*:\s*["\']?([^"\'\s]+)', _read(path), flags=re.I | re.M)
        return [match.group(1)] if match else []
    if path.suffix == ".json":
        data = _json(path)
        values = []
        for key in ["doi", "conceptdoi"]:
            if data.get(key):
                values.append(str(data[key]))
        return values
    return []


def _doi_resolves(doi: str) -> bool:
    request = urllib.request.Request(f"https://doi.org/{doi}", method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return 200 <= int(response.status) < 400
    except Exception:
        return False


def _workflow_job_block(workflow: str, job_name: str) -> str:
    match = re.search(rf"^  {re.escape(job_name)}:\n", workflow, flags=re.M)
    if not match:
        return ""
    start = match.start()
    next_job = re.search(r"^  [A-Za-z0-9_-]+:\n", workflow[match.end():], flags=re.M)
    end = match.end() + next_job.start() if next_job else len(workflow)
    return workflow[start:end]


def _make_target_block(makefile: str, target_name: str) -> str:
    match = re.search(rf"^{re.escape(target_name)}:(?:[^\n]*)\n", makefile, flags=re.M)
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
        version_files["include/qdalign.h"] = _c_header_version(_read(root / "include" / "qdalign.h"))
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


def check_no_unminted_doi_fields(root: Path, result: ReleaseAudit) -> None:
    checked = [
        root / "CITATION.cff",
        root / ".zenodo.json",
        root / "codemeta.json",
    ]
    observed = [doi for path in checked for doi in _doi_values(path)]
    if observed:
        for doi in observed:
            if not _doi_resolves(doi):
                result.failures.append(f"DOI does not resolve through doi.org: {doi}")
        if not any("DOI" in failure for failure in result.failures):
            result.passed.append("DOI fields resolve")
        return
    for path in checked:
        if _has_release_doi_field(path):
            result.failures.append(f"{path.name} has a DOI field before an immutable release DOI is minted")
    if not any("DOI" in failure for failure in result.failures):
        result.passed.append("DOI fields deferred until minted release")


def check_sdist_metadata(root: Path, result: ReleaseAudit) -> None:
    manifest = _read(root / "MANIFEST.in")
    verifier = _read(root / "scripts" / "check_python_wheel.py")
    for required in [
        "CITATION.cff",
        "codemeta.json",
        "docs/assay-evidence.json",
        "src/qdalign.c",
        "src/qdmetal_stub.c",
        "include/qdalign.h",
        "include/qdmetal.h",
    ]:
        if f"include {required}" not in manifest:
            result.failures.append(f"MANIFEST.in must include {required}")
    for required_suffix in [
        "/CITATION.cff",
        "/codemeta.json",
        "/docs/assay-evidence.json",
        "/src/qdalign.c",
        "/src/qdmetal_stub.c",
        "/include/qdalign.h",
        "/include/qdmetal.h",
    ]:
        if required_suffix not in verifier:
            result.failures.append(f"scripts/check_python_wheel.py must verify {required_suffix}")
    for verifier_fragment in ["dotmatch/data/assay-evidence.json", "evidence_boundary"]:
        if verifier_fragment not in verifier:
            result.failures.append(f"scripts/check_python_wheel.py must verify {verifier_fragment}")
    if not any("MANIFEST.in" in failure or "check_python_wheel.py" in failure for failure in result.failures):
        result.passed.append("sdist release metadata verified")


def check_distribution_surfaces(root: Path, result: ReleaseAudit) -> None:
    workflow = _read(root / ".github" / "workflows" / "release.yml")
    pyproject = _read(root / "pyproject.toml")
    dockerfile = _read(root / "Dockerfile")
    bioconda = _read(root / "packaging" / "bioconda" / "meta.yaml")
    packaging = _read(root / "docs" / "packaging.md")
    release_process = _read(root / "docs" / "release-process.md")
    readme = _read(root / "README.md")
    makefile = _read(root / "Makefile")
    project_version = _pyproject_version(root / "pyproject.toml")

    required_workflow_fragments = [
        "id-token: write",
        "packages: write",
        "pypa/gh-action-pypi-publish@release/v1",
        "packages-dir: dist-pypi",
        "docker/metadata-action",
        "docker/login-action",
        "docker/build-push-action",
        "ghcr.io/dnncha/dotmatch",
        "python scripts/check_python_wheel.py --wheel-only --out-dir dist-linux",
        "CIBW_ARCHS_LINUX: \"x86_64 aarch64\"",
        "--require-repaired-linux-architectures x86_64 aarch64",
        "docker/setup-qemu-action@v4",
        "docker buildx build --platform linux/arm64",
        "platforms: linux/amd64,linux/arm64",
        "docker buildx imagetools inspect",
        "scripts/check_oci_manifest.py",
        "docker image inspect dotmatch:ci",
        "SHA256SUMS.txt",
    ]
    for fragment in required_workflow_fragments:
        if fragment not in workflow:
            result.failures.append(f"release workflow missing {fragment}")
    required_cibuildwheel_fragments = [
        "cp39-manylinux_aarch64",
        "cp312-manylinux_aarch64",
        "cp39-musllinux_aarch64",
        "cp312-musllinux_aarch64",
        "test-command",
        "dotmatch leq 1 ACGT AGGT",
    ]
    for fragment in required_cibuildwheel_fragments:
        if fragment not in pyproject:
            result.failures.append(f"pyproject.toml cibuildwheel configuration missing {fragment}")
    if "dotmatch-wheel-Linux" in workflow:
        result.failures.append("release workflow must not publish raw Linux wheels to PyPI")
    container_version_check = re.search(
        r"docker image inspect dotmatch:ci[^\n]*\| grep [\"']\^([^\"']+)\$[\"']",
        workflow,
    )
    if container_version_check and project_version:
        workflow_version = container_version_check.group(1)
        if workflow_version not in {project_version, "${VERSION}"}:
            result.failures.append(
                "release workflow container version smoke test must match "
                f"pyproject.toml ({project_version}); saw {workflow_version}"
            )
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
        if "python -m pip install -r docs/requirements.txt" not in preflight:
            result.failures.append("release workflow preflight job must install documentation tooling")
        if "make test" not in preflight:
            result.failures.append("release workflow preflight job must run make test")
        if "make cli-test" not in preflight:
            result.failures.append("release workflow preflight job must run make cli-test")
        if "make asan" not in preflight:
            result.failures.append("release workflow preflight job must run make asan")
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
    required_published_labels = [
        "org.opencontainers.image.title=DotMatch",
        "org.opencontainers.image.source=https://github.com/dnncha/dotmatch",
        "org.opencontainers.image.url=https://dotmatch.readthedocs.io/",
        "org.opencontainers.image.documentation=https://dotmatch.readthedocs.io/",
        "org.opencontainers.image.licenses=Apache-2.0",
        "org.opencontainers.image.authors=Donncha O'Toole",
    ]
    for label in required_published_labels:
        if label not in workflow:
            result.failures.append(f"release workflow metadata-action missing published OCI label {label}")

    if "REPLACE_WITH_RELEASE_TARBALL_SHA256" not in bioconda:
        result.failures.append("Bioconda template must retain release SHA256 placeholder until copying into bioconda-recipes")
    if "dotmatch dist ACGT AGGT" not in bioconda:
        result.failures.append("Bioconda template must include native CLI smoke test")
    if "additional-platforms:" not in bioconda or "- osx-arm64" not in bioconda:
        result.failures.append("Bioconda template must opt into osx-arm64 / Apple Silicon builds")

    for label, text in [("docs/packaging.md", packaging), ("docs/release-process.md", release_process)]:
        if not text.strip():
            result.failures.append(f"{label} is empty")
    if "osx-arm64" not in packaging or "Apple Silicon" not in packaging:
        result.failures.append("docs/packaging.md must document Bioconda osx-arm64 / Apple Silicon support")
    if "osx-arm64" not in readme or "Apple Silicon" not in readme:
        result.failures.append("README.md must document Bioconda osx-arm64 / Apple Silicon support")
    for release_process_fragment in [
        "make asan",
        "make docs-ready",
        "make scientific-readiness-ready",
        "make native-exact-gate",
        "make pretag-ready",
    ]:
        if release_process_fragment not in release_process:
            result.failures.append(f"docs/release-process.md must mention {release_process_fragment}")
    docs_block = _make_target_block(makefile, "docs-ready")
    if not docs_block:
        result.failures.append("Makefile must include docs-ready target")
    elif "python3 -m sphinx -W -b html docs docs/_build/html" not in docs_block:
        result.failures.append("Makefile docs-ready target must build public docs with Sphinx warnings as errors")
    release_block = _make_target_block(makefile, "release-ready")
    if not release_block:
        result.failures.append("Makefile must include release-ready target")
    elif "docs-ready" not in release_block:
        result.failures.append("Makefile release-ready target must include docs-ready")
    elif "native-exact-gate" not in release_block:
        result.failures.append("Makefile release-ready target must include native-exact-gate")
    crispr_block = _make_target_block(makefile, "crispr-comparison-gate")
    if not crispr_block:
        result.failures.append("Makefile must include crispr-comparison-gate target")
    else:
        for required_k in ["2", "3"]:
            fragment = f"--require-hamming-k23-comparator {required_k}"
            if fragment not in crispr_block:
                result.failures.append(f"Makefile crispr-comparison-gate target must include {fragment}")
    repository_block = _make_target_block(makefile, "repository-ready")
    if not repository_block:
        result.failures.append("Makefile must include repository-ready target")
    elif "$(MAKE) docs-ready" not in repository_block:
        result.failures.append("Makefile repository-ready target must include $(MAKE) docs-ready")
    pretag_block = _make_target_block(makefile, "pretag-ready")
    if not pretag_block:
        result.failures.append("Makefile must include pretag-ready target")
    for pretag_fragment in [
        "$(MAKE) test",
        "$(MAKE) cli-test",
        "$(MAKE) asan",
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
    check_no_unminted_doi_fields(root, result)
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
