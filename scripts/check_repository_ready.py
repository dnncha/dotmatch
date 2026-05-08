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
    "CITATION.cff",
    "codemeta.json",
    "CONTRIBUTING.md",
    "CODE_OF_CONDUCT.md",
    "SECURITY.md",
    "SUPPORT.md",
    "pyproject.toml",
    "package.json",
    "MANIFEST.in",
    "setup.py",
    ".gitignore",
    ".gitattributes",
    ".editorconfig",
    ".github/workflows/ci.yml",
    ".github/workflows/codeql.yml",
    ".github/workflows/release.yml",
    ".github/PULL_REQUEST_TEMPLATE.md",
    ".github/ISSUE_TEMPLATE/bug_report.yml",
    ".github/ISSUE_TEMPLATE/feature_request.yml",
    ".github/ISSUE_TEMPLATE/benchmark_evidence.yml",
    ".zenodo.json",
    "docs/scientific-claims.md",
    "docs/assay-evidence.json",
    "docs/distribution-release.json",
    "docs/distribution-submission.md",
    "docs/workflow-adoption.json",
    "docs/release-process.md",
    "docs/methods-and-citation.md",
    "docs/native-comparator-scope.md",
    "docs/packaging.md",
    "docs/schemas.md",
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
    "scripts/check_alphabet_policy.py",
    "scripts/check_citation_metadata.py",
    "scripts/check_distribution_channels.py",
    "scripts/check_distribution_record.py",
    "scripts/check_distribution_submission.py",
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
        return [path for item in proc.stdout.split(b"\0") if item and (path := root / item.decode()).exists()]

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
    readme = (root / "README.md").read_text(encoding="utf-8")
    evidence_path = root / "docs" / "scientific-claims.md"
    evidence = evidence_path.read_text(encoding="utf-8") if evidence_path.exists() else ""
    native_scope_path = root / "docs" / "native-comparator-scope.md"
    native_scope = native_scope_path.read_text(encoding="utf-8") if native_scope_path.exists() else ""

    if "docs/scientific-claims.md" not in readme:
        result.failures.append("README.md must link to docs/scientific-claims.md")
    if "make pretag-ready" not in readme:
        result.failures.append("README.md must document the consolidated local pre-tag release gate")
    if "make distribution-channels" not in readme or "make workflow-adoption-status" not in readme:
        result.failures.append("README.md must keep post-release and workflow-adoption gates separate from pretag-ready")
    if "barcode-comparison-gate" not in evidence or "requires real-data rows" not in evidence:
        result.failures.append("docs/scientific-claims.md must document the barcode-comparison-gate evidence boundary")
    if "bcl-comparison-gate" not in evidence or "requires real run folders" not in evidence:
        result.failures.append("docs/scientific-claims.md must document the bcl-comparison-gate evidence boundary")
    if "not a genome aligner" not in evidence and "General aligner replacement" not in evidence:
        result.failures.append("docs/scientific-claims.md must document the general-aligner evidence boundary")
    native_scope_fragments = [
        "SeqAn and Parasail comparisons are not claimed as completed release evidence yet",
        "equivalent global edit-distance or documented semi-global scoring semantics",
        "fixed threshold `k`",
        "zero assignment mismatches",
        "limited to Edlib exhaustive global edit-distance assignment scans",
    ]
    for fragment in native_scope_fragments:
        if fragment not in native_scope:
            result.failures.append(f"docs/native-comparator-scope.md must document native comparator scope: {fragment}")
    if not any("evidence boundary" in failure for failure in result.failures):
        result.passed.append("evidence boundaries documented")


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
    if "does not broaden README/docs claims beyond checked evidence" not in template:
        result.failures.append("pull request template must preserve claim-boundary checklist")
    if not any("pull request template" in failure for failure in result.failures):
        result.passed.append("pull request template evidence checklist present")


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


def audit(root: Path) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    check_required_files(root, result)
    check_metadata(root, result)
    check_release_versions(root, result)
    check_evidence_docs(root, result)
    check_manifest(root, result)
    check_pull_request_template(root, result)
    check_repository_tree(root, result)
    check_no_local_absolute_paths(root, result)
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
