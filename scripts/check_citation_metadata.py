#!/usr/bin/env python3
"""Audit citation and discovery metadata for release consistency."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


REQUIRED_KEYWORDS = [
    "bioinformatics",
    "computational biology",
    "CRISPR",
    "FASTQ",
    "known-target assignment",
]
REPOSITORY_URL = "https://github.com/dnncha/dotmatch"
AUTHOR_GIVEN = "Donncha"
AUTHOR_FAMILY = "O'Toole"


class AuditResult:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _read_json(path: Path) -> dict:
    return json.loads(_read(path))


def _project_version(root: Path) -> str:
    text = _read(root / "pyproject.toml")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, flags=re.M)
    return match.group(1) if match else ""


def _project_license(root: Path) -> str:
    text = _read(root / "pyproject.toml")
    match = re.search(r'^license\s*=\s*"([^"]+)"', text, flags=re.M)
    return match.group(1) if match else ""


def _toml_section(text: str, section: str) -> str:
    match = re.search(rf"^\[{re.escape(section)}\]\s*$", text, flags=re.M)
    if not match:
        return ""
    start = match.end()
    next_section = re.search(r"^\[[^\]]+\]\s*$", text[start:], flags=re.M)
    end = start + next_section.start() if next_section else len(text)
    return text[start:end]


def _toml_scalar(section: str, key: str) -> str:
    match = re.search(rf'^{re.escape(key)}\s*=\s*"([^"]*)"', section, flags=re.M)
    return match.group(1) if match else ""


def _toml_array(section: str, key: str) -> list[str]:
    match = re.search(rf"^{re.escape(key)}\s*=\s*\[(.*?)\]", section, flags=re.M | re.S)
    if not match:
        return []
    return re.findall(r'"([^"]+)"', match.group(1))


def _toml_table(section: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for key, value in re.findall(r'^([A-Za-z0-9_-]+)\s*=\s*"([^"]*)"', section, flags=re.M):
        values[key] = value
    return values


def _unquote(value: str) -> str:
    text = value.strip()
    if len(text) >= 2 and text[0] in {"'", '"'} and text[-1] == text[0]:
        return text[1:-1]
    return text


def _cff_scalar(text: str, key: str) -> str:
    match = re.search(rf"^{re.escape(key)}:\s*(.+)$", text, flags=re.M)
    return _unquote(match.group(1)) if match else ""


def _cff_keywords(text: str) -> list[str]:
    match = re.search(r"^keywords:\s*$((?:\n\s+-\s+.+)+)", text, flags=re.M)
    if not match:
        return []
    values = []
    for line in match.group(1).splitlines():
        item = line.strip()
        if item.startswith("- "):
            values.append(_unquote(item[2:]))
    return values


def _contains_unminted_doi_claim(path: Path, *, codemeta_context: bool = False) -> bool:
    if path.name == "codemeta.json" and codemeta_context:
        data = _read_json(path)
        data.pop("@context", None)
        text = json.dumps(data)
    else:
        text = _read(path)
    return "doi.org" in text.lower() or re.search(r"^\s*doi\s*:", text, flags=re.I | re.M) is not None


def _check_keywords(source: str, keywords: list[str], result: AuditResult) -> None:
    observed = {keyword.lower() for keyword in keywords}
    for keyword in REQUIRED_KEYWORDS:
        if keyword.lower() not in observed:
            result.failures.append(f"{source} missing discovery keyword: {keyword}")


def _check_pyproject_discovery(root: Path, result: AuditResult) -> None:
    text = _read(root / "pyproject.toml")
    project = _toml_section(text, "project")
    urls = _toml_table(_toml_section(text, "project.urls"))
    description = _toml_scalar(project, "description")
    keywords = _toml_array(project, "keywords")
    classifiers = _toml_array(project, "classifiers")

    _check_keywords("pyproject.toml", keywords, result)
    description_lc = description.lower()
    for fragment, label in [
        ("known-target short-dna assignment", "known-target short-DNA assignment"),
        ("crispr", "CRISPR"),
        ("barcode", "barcode"),
        ("fastq", "FASTQ"),
    ]:
        if fragment not in description_lc:
            result.failures.append(f"pyproject.toml description must mention {label}")
    if "Topic :: Scientific/Engineering :: Bio-Informatics" not in classifiers:
        result.failures.append("pyproject.toml classifiers must include Topic :: Scientific/Engineering :: Bio-Informatics")
    if "Intended Audience :: Science/Research" not in classifiers:
        result.failures.append("pyproject.toml classifiers must include Intended Audience :: Science/Research")
    for key in ["Homepage", "Repository"]:
        if urls.get(key) != REPOSITORY_URL:
            result.failures.append(f"pyproject.toml project URLs must include {key}")
    if urls.get("Issues") != f"{REPOSITORY_URL}/issues":
        result.failures.append("pyproject.toml project URLs must include Issues")
    if not urls.get("Documentation", "").startswith(REPOSITORY_URL):
        result.failures.append("pyproject.toml project URLs must include Documentation")


def _check_versions(root: Path, citation: dict, codemeta: dict, zenodo: dict, result: AuditResult) -> None:
    versions = {
        "pyproject.toml": _project_version(root),
        "CITATION.cff": str(citation.get("version") or ""),
        "codemeta.json version": str(codemeta.get("version") or ""),
        "codemeta.json softwareVersion": str(codemeta.get("softwareVersion") or ""),
        ".zenodo.json": str(zenodo.get("version") or ""),
    }
    missing = [name for name, version in versions.items() if not version]
    for name in missing:
        result.failures.append(f"{name} must declare version")
    unique = sorted({version for version in versions.values() if version})
    if len(unique) > 1:
        detail = ", ".join(f"{name}={version}" for name, version in sorted(versions.items()) if version)
        result.failures.append(f"citation metadata version mismatch: {detail}")


def audit(root: Path) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    try:
        citation_text = _read(root / "CITATION.cff")
        citation = {
            "title": _cff_scalar(citation_text, "title"),
            "version": _cff_scalar(citation_text, "version"),
            "repository-code": _cff_scalar(citation_text, "repository-code"),
            "license": _cff_scalar(citation_text, "license"),
            "abstract": _cff_scalar(citation_text, "abstract"),
            "keywords": _cff_keywords(citation_text),
        }
        codemeta = _read_json(root / "codemeta.json")
        zenodo = _read_json(root / ".zenodo.json")
    except Exception as exc:
        result.failures.append(f"citation metadata could not be read: {exc}")
        return result

    _check_versions(root, citation, codemeta, zenodo, result)

    cff_title = str(citation.get("title") or "")
    zenodo_title = str(zenodo.get("title") or "")
    if not cff_title or not zenodo_title or cff_title != zenodo_title:
        result.failures.append(f"citation metadata title mismatch: CITATION.cff={cff_title!r}, .zenodo.json={zenodo_title!r}")

    pyproject_license = _project_license(root)
    if pyproject_license != "Apache-2.0":
        result.failures.append("pyproject.toml license must be Apache-2.0")
    if citation.get("license") != "Apache-2.0":
        result.failures.append("CITATION.cff license must be Apache-2.0")
    if str(codemeta.get("license") or "").find("Apache-2.0") == -1:
        result.failures.append("codemeta.json license must reference Apache-2.0")
    if zenodo.get("license") != "Apache-2.0":
        result.failures.append(".zenodo.json license must be Apache-2.0")

    if citation.get("repository-code") != REPOSITORY_URL:
        result.failures.append("CITATION.cff repository-code must point to the canonical repository")
    if codemeta.get("codeRepository") != REPOSITORY_URL or codemeta.get("url") != REPOSITORY_URL:
        result.failures.append("codemeta.json repository URLs must point to the canonical repository")
    related_urls = [item.get("identifier") for item in zenodo.get("related_identifiers", []) if isinstance(item, dict)]
    if REPOSITORY_URL not in related_urls:
        result.failures.append(".zenodo.json related_identifiers must include the canonical repository")

    if AUTHOR_GIVEN not in citation_text or AUTHOR_FAMILY not in citation_text:
        result.failures.append("CITATION.cff must include the release author")
    authors = codemeta.get("author") or []
    if not any(item.get("givenName") == AUTHOR_GIVEN and item.get("familyName") == AUTHOR_FAMILY for item in authors if isinstance(item, dict)):
        result.failures.append("codemeta.json must include the release author")
    creators = zenodo.get("creators") or []
    if not any(AUTHOR_FAMILY in str(item.get("name") or "") and AUTHOR_GIVEN in str(item.get("name") or "") for item in creators if isinstance(item, dict)):
        result.failures.append(".zenodo.json must include the release creator")

    _check_keywords("CITATION.cff", list(citation.get("keywords") or []), result)
    _check_keywords("codemeta.json", [str(value) for value in codemeta.get("keywords") or []], result)
    _check_keywords(".zenodo.json", [str(value) for value in zenodo.get("keywords") or []], result)
    _check_pyproject_discovery(root, result)

    abstract = str(citation.get("abstract") or "")
    for fragment in ["known-target short-DNA assignment", "CRISPR", "barcode"]:
        if fragment.lower() not in abstract.lower():
            result.failures.append(f"CITATION.cff abstract must mention {fragment}")

    for path, allow_context in [
        (root / "CITATION.cff", False),
        (root / "codemeta.json", True),
        (root / ".zenodo.json", False),
    ]:
        if _contains_unminted_doi_claim(path, codemeta_context=allow_context):
            result.failures.append(f"{path.name} must not claim a DOI before Zenodo release")

    if result.ok:
        result.passed.append("citation metadata aligned")
        result.passed.append("citation discovery metadata complete")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="repository root")
    args = parser.parse_args()

    result = audit(Path(args.root))
    for item in result.passed:
        print(f"PASS: {item}")
    for item in result.failures:
        print(f"FAIL: {item}")
    if result.ok:
        print("CITATION METADATA: PASS")
        return 0
    print("CITATION METADATA: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
