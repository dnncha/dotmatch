#!/usr/bin/env python3

import argparse
import json
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


MAX_FILE_BYTES = 25 * 1024 * 1024
MAX_TOTAL_BYTES = 100 * 1024 * 1024

REQUIRED_FILES = [
    "README.md",
    "CHANGELOG.md",
    "LICENSE",
    "NOTICE",
    "CITATION.cff",
    "codemeta.json",
    "CONTRIBUTING.md",
    "CODE_OF_CONDUCT.md",
    "SECURITY.md",
    "SUPPORT.md",
    "TRADEMARKS.md",
    "pyproject.toml",
    "package.json",
    "MANIFEST.in",
    "setup.py",
    ".gitignore",
    ".gitattributes",
    ".editorconfig",
    ".github/workflows/ci.yml",
    ".github/workflows/codeql.yml",
    ".github/workflows/pages.yml",
    ".github/workflows/release.yml",
    ".github/PULL_REQUEST_TEMPLATE.md",
    ".github/ISSUE_TEMPLATE/bug_report.yml",
    ".github/ISSUE_TEMPLATE/feature_request.yml",
    ".github/ISSUE_TEMPLATE/benchmark_evidence.yml",
    ".zenodo.json",
    "docs/scientific-claims.md",
    "docs/scientific-readiness.json",
    "docs/assay-evidence.json",
    "docs/distribution-release.json",
    "docs/workflow-adoption.json",
    "docs/release-process.md",
    "docs/methods-and-citation.md",
    "docs/native-comparator-scope.md",
    "docs/packaging.md",
    "docs/schemas.md",
    "docs/commercial-boundary.md",
    "docs/evidence-packet-v1.md",
    "examples/workflows/galaxy/dotmatch_crispr_count.xml",
    "examples/workflows/multiqc/multiqc_config.yaml",
    "examples/workflows/nf-core/README.md",
    "examples/workflows/nf-core/modules/local/dotmatch/crispr_count/main.nf",
    "examples/workflows/nf-core/modules/local/dotmatch/crispr_count/meta.yml",
    "examples/workflows/nextflow/main.nf",
    "examples/workflows/snakemake/Snakefile",
    "packaging/bioconda/meta.yaml",
    "packaging/bioconda/build.sh",
    "scripts/check_assay_evidence.py",
    "scripts/check_scientific_readiness.py",
    "scripts/check_alphabet_policy.py",
    "scripts/check_citation_metadata.py",
    "scripts/check_distribution_channels.py",
    "scripts/check_distribution_record.py",
    "scripts/check_bioconda_recipe.py",
    "scripts/check_native_comparator_scope.py",
    "scripts/check_workflow_adoption.py",
    "src/qdalign.c",
    "include/qdalign.h",
]

GENERATED_PATH_PARTS = {
    ".next",
    ".pytest_cache",
    "__pycache__",
    "assay_out",
    "node_modules",
    "build",
    "dist",
}

LOCAL_ABSOLUTE_PATH_PREFIXES = [
    "/" + "Users/",
    "/" + "private/tmp/",
    "/" + "var/folders/",
    "/" + "tmp/dotmatch",
]
LOCAL_ABSOLUTE_PATH_PATTERNS = [
    re.compile(re.escape(prefix.encode("utf-8")) + rb"[^,\s\"')<]*")
    for prefix in LOCAL_ABSOLUTE_PATH_PREFIXES
]

RAW_DATA_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".bam", ".bcl")
ALLOWED_RAW_DATA_PREFIXES = (
    "demo-data/",
    "examples/barcode_autopsy/failure_modes/",
    "examples/workflows/fixtures/",
    "examples/workflows/galaxy/test-data/",
    "examples/workflows/nf-core/upstream/modules/nf-core/dotmatch/",
    "benchmarks/real/data/",
)
PRIVATE_DATA_MARKERS = [
    re.compile(rb"\bcustomer\b", re.IGNORECASE),
    re.compile(rb"\bpatient\b", re.IGNORECASE),
    re.compile(rb"\bPHI\b", re.IGNORECASE),
    re.compile(rb"\bproprietary\b", re.IGNORECASE),
]


@dataclass
class AuditResult:
    passed: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    failures: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.failures


def rel(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def in_git_worktree(root: Path) -> bool:
    return (
        subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--is-inside-work-tree"],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=False,
        ).stdout.strip()
        == "true"
    )


def repository_files(root: Path) -> list[Path]:
    if in_git_worktree(root):
        proc = subprocess.run(
            ["git", "-C", str(root), "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
        return [path for item in proc.stdout.split(b"\0") if item and (path := root / item.decode()).is_file()]

    files: list[Path] = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        parts = set(path.relative_to(root).parts)
        if parts & GENERATED_PATH_PARTS:
            continue
        files.append(path)
    return files


def check_required_files(root: Path, result: AuditResult) -> None:
    missing = [path for path in REQUIRED_FILES if not (root / path).is_file()]
    if missing:
        result.failures.extend(f"missing required file: {path}" for path in missing)
    else:
        result.passed.append("required files present")


def check_metadata(root: Path, result: AuditResult) -> None:
    try:
        codemeta = json.loads((root / "codemeta.json").read_text(encoding="utf-8"))
    except Exception as exc:
        result.failures.append(f"codemeta.json is invalid JSON: {exc}")
        return

    if codemeta.get("name") != "DotMatch":
        result.failures.append("codemeta.json name must be DotMatch")
    if "Apache-2.0" not in str(codemeta.get("license", "")):
        result.failures.append("codemeta.json must reference Apache-2.0")
    if not codemeta.get("softwareVersion"):
        result.failures.append("codemeta.json must declare softwareVersion")
    if not codemeta.get("keywords"):
        result.failures.append("codemeta.json must include discovery keywords")
    if not result.failures:
        result.passed.append("metadata files parse")


def _pyproject_version(path: Path) -> Optional[str]:
    in_project = False
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            in_project = line == "[project]"
            continue
        if not in_project or not line.startswith("version"):
            continue
        match = re.match(r'version\s*=\s*["\']([^"\']+)["\']', line)
        if match:
            return match.group(1)
    return None


def _cff_version(path: Path) -> Optional[str]:
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        match = re.match(r'\s*version\s*:\s*["\']?([^"\']+)["\']?\s*$', raw_line)
        if match:
            return match.group(1).strip()
    return None


def check_release_versions(root: Path, result: AuditResult) -> None:
    versions: dict[str, Optional[str]] = {}
    try:
        versions["pyproject.toml"] = _pyproject_version(root / "pyproject.toml")
    except Exception as exc:
        result.failures.append(f"pyproject.toml version could not be read: {exc}")
    try:
        versions["package.json"] = json.loads((root / "package.json").read_text(encoding="utf-8")).get("version")
    except Exception as exc:
        result.failures.append(f"package.json version could not be read: {exc}")
    try:
        versions["codemeta.json"] = json.loads((root / "codemeta.json").read_text(encoding="utf-8")).get("version")
    except Exception as exc:
        result.failures.append(f"codemeta.json version could not be read: {exc}")
    try:
        versions["CITATION.cff"] = _cff_version(root / "CITATION.cff")
    except Exception as exc:
        result.failures.append(f"CITATION.cff version could not be read: {exc}")

    missing = [name for name, version in versions.items() if not version]
    result.failures.extend(f"{name} must declare release version" for name in missing)

    declared = {name: version for name, version in versions.items() if version}
    unique_versions = sorted(set(declared.values()))
    if len(unique_versions) > 1:
        detail = ", ".join(f"{name}={version}" for name, version in sorted(declared.items()))
        result.failures.append(f"release version mismatch: {detail}")

    readme = (root / "README.md").read_text(encoding="utf-8")
    if "v0.1.0-dev" in readme or "0.1.0-dev" in readme:
        result.failures.append("README.md must not describe the release version as dev")

    if not any("version" in failure or "README.md must not describe" in failure for failure in result.failures):
        result.passed.append("release versions aligned")


def check_evidence_docs(root: Path, result: AuditResult) -> None:
    evidence_path = root / "docs" / "scientific-claims.md"
    readiness_path = root / "docs" / "scientific-readiness.json"
    native_scope_path = root / "docs" / "native-comparator-scope.md"
    for path in [evidence_path, readiness_path, native_scope_path, root / "docs" / "release-process.md"]:
        if not path.is_file():
            result.failures.append(f"{path.relative_to(root).as_posix()} missing")
        elif not path.read_text(encoding="utf-8").strip():
            result.failures.append(f"{path.relative_to(root).as_posix()} is empty")
    if not any("docs/" in failure for failure in result.failures):
        result.passed.append("evidence and release docs present")


def check_readme_distribution_status(root: Path, result: AuditResult) -> None:
    failures_before = len(result.failures)

    if not (root / "docs" / "packaging.md").is_file():
        result.failures.append("docs/packaging.md missing")
    try:
        json.loads((root / "docs" / "distribution-release.json").read_text(encoding="utf-8"))
    except Exception as exc:
        result.failures.append(f"docs/distribution-release.json could not be parsed: {exc}")

    if len(result.failures) == failures_before:
        result.passed.append("distribution metadata present")


def check_manifest(root: Path, result: AuditResult) -> None:
    manifest = (root / "MANIFEST.in").read_text(encoding="utf-8")
    for required in ["include/qdalign.h", "src/qdalign.c"]:
        if required not in manifest:
            result.failures.append(f"MANIFEST.in must include {required}")
    if "include/qdalign.h" in manifest and "src/qdalign.c" in manifest:
        result.passed.append("sdist native sources listed")


def check_pull_request_template(root: Path, result: AuditResult) -> None:
    template = (root / ".github" / "PULL_REQUEST_TEMPLATE.md").read_text(encoding="utf-8")
    if "make test" not in template:
        result.failures.append("pull request template must require make test evidence")
    if "make cli-test" not in template:
        result.failures.append("pull request template must require make cli-test evidence")
    if "make python-test" not in template:
        result.failures.append("pull request template must require make python-test evidence")
    if "make pretag-ready" not in template:
        result.failures.append("pull request template must mention make pretag-ready for release-surface changes")
    if "make asan" not in template:
        result.failures.append("pull request template must mention make asan for native safety changes")
    if "make scientific-readiness-ready" not in template:
        result.failures.append("pull request template must mention make scientific-readiness-ready for claim/evidence changes")
    if "make repository-ready" not in template:
        result.failures.append(
            "pull request template must mention make repository-ready for governance/security/licensing/trademark/public-doc changes"
        )
    if not any("pull request template" in failure for failure in result.failures):
        result.passed.append("pull request template evidence checklist present")


def check_pages_workflow(root: Path, result: AuditResult) -> None:
    failures_before = len(result.failures)
    path = root / ".github" / "workflows" / "pages.yml"
    try:
        workflow = path.read_text(encoding="utf-8")
    except OSError as exc:
        result.failures.append(f".github/workflows/pages.yml could not be read: {exc}")
        return

    for fragment in [
        "actions/configure-pages@v6",
        "actions/upload-pages-artifact@v5",
        "actions/deploy-pages@v4",
        "pages: write",
        "id-token: write",
    ]:
        if fragment not in workflow:
            result.failures.append(f"pages workflow must include: {fragment}")
    if "repos/${GITHUB_REPOSITORY}/pages/deployments" in workflow:
        result.failures.append("pages workflow must use actions/deploy-pages instead of a custom deployment API call")
    if len(result.failures) == failures_before:
        result.passed.append("GitHub Pages workflow uses the supported deployment action")


def check_release_workflow(root: Path, result: AuditResult) -> None:
    path = root / ".github" / "workflows" / "release.yml"
    try:
        workflow = path.read_text(encoding="utf-8")
    except OSError as exc:
        result.failures.append(f".github/workflows/release.yml could not be read: {exc}")
        return

    if "sha256sum *.whl *.tar.gz > SHA256SUMS.txt" not in workflow:
        result.failures.append(
            "release workflow must checksum only downloadable wheel and source-distribution assets"
        )
        return
    if "sha256sum * > SHA256SUMS.txt" in workflow:
        result.failures.append("release workflow must not checksum internal build records")
        return
    result.passed.append("release checksum manifest contains only downloadable artifacts")


def _require_text(path: Path, needles: list[str], result: AuditResult) -> None:
    rel_path = path.relative_to(path.parents[1]).as_posix() if path.parent.name == "docs" else path.name
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        result.failures.append(f"{rel_path} could not be read: {exc}")
        return
    for needle in needles:
        if needle not in text:
            result.failures.append(f"{rel_path} must include: {needle}")


def check_open_core_governance(root: Path, result: AuditResult) -> None:
    failures_before = len(result.failures)

    try:
        license_text = (root / "LICENSE").read_text(encoding="utf-8")
        if "Apache License" not in license_text or "Version 2.0" not in license_text:
            result.failures.append("LICENSE must remain Apache License 2.0")
    except OSError as exc:
        result.failures.append(f"LICENSE could not be read: {exc}")

    _require_text(
        root / "CONTRIBUTING.md",
        ["same Apache-2.0 terms", "unless otherwise agreed"],
        result,
    )
    _require_text(
        root / "SECURITY.md",
        [
            "Do not attach real FASTQ, BAM, BCL, customer assay data",
            "Use synthetic or minimized reproductions",
            "Report security and data-leak issues privately",
        ],
        result,
    )
    _require_text(
        root / "TRADEMARKS.md",
        [
            "DotMatch",
            "DotMatch Pro",
            "do not change the Apache-2.0 license",
        ],
        result,
    )
    _require_text(
        root / "NOTICE",
        ["Apache License, Version 2.0", "DotMatch Pro"],
        result,
    )
    _require_text(
        root / "README.md",
        [
            "dotmatch-community",
            "Apache-2.0",
        ],
        result,
    )
    _require_text(
        root / "docs" / "commercial-boundary.md",
        [
            "the deterministic assignment engine and CLI",
            "hosted or team workspaces",
            "run registries",
            "signed reports",
            "private assay registries",
            "enterprise connectors",
            "commercial support",
            "license present and still Apache-2.0",
            "security policy present",
            "no raw customer assay data",
            "docs build passes with `make docs-ready`",
            "tests pass",
        ],
        result,
    )
    _require_text(
        root / "docs" / "evidence-packet-v1.md",
        ["unique", "ambiguous", "none", "invalid", "private customer FASTQ/BAM/BCL"],
        result,
    )

    if len(result.failures) == failures_before:
        result.passed.append("open-core governance files present")


def check_repository_tree(root: Path, result: AuditResult) -> None:
    files = repository_files(root)
    total = 0
    for path in files:
        relative = rel(path, root)
        size = path.stat().st_size
        total += size
        if path.name == ".DS_Store" or "__MACOSX" in path.parts or path.name.startswith("._"):
            result.failures.append(f"generated macOS metadata is tracked: {relative}")
        if size > MAX_FILE_BYTES:
            result.failures.append(f"repository file exceeds {MAX_FILE_BYTES} bytes: {relative}")
        if any(part in GENERATED_PATH_PARTS for part in path.relative_to(root).parts):
            result.failures.append(f"generated path is tracked: {relative}")
    if total > MAX_TOTAL_BYTES:
        result.failures.append(f"repository tree exceeds {MAX_TOTAL_BYTES} bytes")
    if not any("repository file exceeds" in failure or "generated" in failure for failure in result.failures):
        result.passed.append(f"repository tree size ok ({total} bytes)")


def check_no_local_absolute_paths(root: Path, result: AuditResult) -> None:
    offenders: list[str] = []
    for path in repository_files(root):
        relative = rel(path, root)
        try:
            data = path.read_bytes()
        except OSError as exc:
            result.failures.append(f"could not read repository file {relative}: {exc}")
            continue
        for pattern in LOCAL_ABSOLUTE_PATH_PATTERNS:
            if pattern.search(data):
                offenders.append(relative)
                break
    if offenders:
        for relative in offenders[:20]:
            result.failures.append(f"repository file contains local absolute path: {relative}")
        if len(offenders) > 20:
            result.failures.append(f"repository file contains local absolute path: {len(offenders) - 20} more files")
    else:
        result.passed.append("no local absolute paths in repository files")


def _has_raw_data_suffix(relative: str) -> bool:
    return relative.endswith(RAW_DATA_SUFFIXES)


def check_no_private_raw_data(root: Path, result: AuditResult) -> None:
    offenders: list[str] = []
    marker_offenders: list[str] = []
    for path in repository_files(root):
        relative = rel(path, root)
        if not _has_raw_data_suffix(relative):
            continue
        if not relative.startswith(ALLOWED_RAW_DATA_PREFIXES):
            offenders.append(relative)
            continue
        try:
            sample = path.read_bytes()[:4096]
        except OSError as exc:
            result.failures.append(f"could not read raw-data fixture {relative}: {exc}")
            continue
        if any(pattern.search(sample) for pattern in PRIVATE_DATA_MARKERS):
            marker_offenders.append(relative)

    for relative in offenders[:20]:
        result.failures.append(f"raw biological data fixture is outside approved public/synthetic paths: {relative}")
    if len(offenders) > 20:
        result.failures.append(f"raw biological data fixture is outside approved paths: {len(offenders) - 20} more files")
    for relative in marker_offenders[:20]:
        result.failures.append(f"raw-data fixture contains private/customer marker: {relative}")
    if len(marker_offenders) > 20:
        result.failures.append(f"raw-data fixture contains private/customer marker: {len(marker_offenders) - 20} more files")
    if not offenders and not marker_offenders:
        result.passed.append("no private raw-data fixtures detected")


def audit(root: Path) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    check_required_files(root, result)
    check_metadata(root, result)
    check_release_versions(root, result)
    check_evidence_docs(root, result)
    check_readme_distribution_status(root, result)
    check_manifest(root, result)
    check_pull_request_template(root, result)
    check_pages_workflow(root, result)
    check_release_workflow(root, result)
    check_open_core_governance(root, result)
    check_repository_tree(root, result)
    check_no_local_absolute_paths(root, result)
    check_no_private_raw_data(root, result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit DotMatch GitHub repository readiness.")
    parser.add_argument("--root", default=".", help="repository root")
    args = parser.parse_args()

    result = audit(Path(args.root))
    for item in result.passed:
        print(f"PASS: {item}")
    for item in result.warnings:
        print(f"WARN: {item}")
    for item in result.failures:
        print(f"FAIL: {item}")
    if result.ok:
        print("REPOSITORY READINESS: PASS")
        return 0
    print("REPOSITORY READINESS: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
