import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_citation_metadata.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_citation_metadata", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_repo(root: Path, *, version: str = "0.1.0") -> None:
    files = {
        "pyproject.toml": (
            "[project]\n"
            'name = "dotmatch"\n'
            f'version = "{version}"\n'
            'description = "Deterministic known-target short-DNA assignment for CRISPR guide counting and barcode FASTQ workflows"\n'
            'license = "Apache-2.0"\n'
            'keywords = ["bioinformatics", "computational biology", "CRISPR", "FASTQ", "known-target assignment"]\n'
            "classifiers = [\n"
            '    "Intended Audience :: Science/Research",\n'
            '    "Topic :: Scientific/Engineering :: Bio-Informatics",\n'
            "]\n"
            "[project.urls]\n"
            'Homepage = "https://github.com/dnncha/dotmatch"\n'
            'Repository = "https://github.com/dnncha/dotmatch"\n'
            'Issues = "https://github.com/dnncha/dotmatch/issues"\n'
            'Documentation = "https://github.com/dnncha/dotmatch#readme"\n'
        ),
        "CITATION.cff": (
            "cff-version: 1.2.0\n"
            'title: "DotMatch: Streaming Exact One-Edit Barcode and Guide Assignment Without Exhaustive Scanning"\n'
            'message: "If you use DotMatch, please cite this software release."\n'
            "type: software\n"
            "authors:\n"
            '  - given-names: "Donncha"\n'
            '    family-names: "O\'Toole"\n'
            'repository-code: "https://github.com/dnncha/dotmatch"\n'
            "license: Apache-2.0\n"
            f'version: "{version}"\n'
            'abstract: "DotMatch is a deterministic known-target short-DNA assignment engine for CRISPR guide counting, barcode demultiplexing, and fixed-target short-read assays."\n'
            "keywords:\n"
            "  - bioinformatics\n"
            "  - computational biology\n"
            "  - CRISPR\n"
            "  - FASTQ\n"
            "  - known-target assignment\n"
        ),
        "codemeta.json": json.dumps(
            {
                "@context": "https://doi.org/10.5063/schema/codemeta-2.0",
                "@type": "SoftwareSourceCode",
                "name": "DotMatch",
                "description": "Fast exact short-DNA known-target assignment for CRISPR guides and barcodes.",
                "url": "https://github.com/dnncha/dotmatch",
                "codeRepository": "https://github.com/dnncha/dotmatch",
                "license": "https://spdx.org/licenses/Apache-2.0",
                "version": version,
                "softwareVersion": version,
                "keywords": [
                    "bioinformatics",
                    "computational biology",
                    "CRISPR",
                    "FASTQ",
                    "known-target assignment",
                ],
                "author": [{"@type": "Person", "givenName": "Donncha", "familyName": "O'Toole"}],
            },
            indent=2,
        )
        + "\n",
        ".zenodo.json": json.dumps(
            {
                "title": "DotMatch: Streaming Exact One-Edit Barcode and Guide Assignment Without Exhaustive Scanning",
                "upload_type": "software",
                "version": version,
                "creators": [{"name": "O'Toole, Donncha"}],
                "description": "DotMatch is a deterministic known-target short-DNA assignment engine for CRISPR guides and barcodes.",
                "license": "Apache-2.0",
                "access_right": "open",
                "keywords": [
                    "bioinformatics",
                    "computational biology",
                    "CRISPR",
                    "FASTQ",
                    "known-target assignment",
                ],
                "related_identifiers": [
                    {
                        "identifier": "https://github.com/dnncha/dotmatch",
                        "relation": "isSupplementTo",
                        "scheme": "url",
                        "resource_type": "software",
                    }
                ],
            },
            indent=2,
        )
        + "\n",
    }
    for path, text in files.items():
        full = root / path
        full.parent.mkdir(parents=True, exist_ok=True)
        full.write_text(text, encoding="utf-8")


def test_citation_metadata_accepts_aligned_metadata(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)

    result = checker.audit(tmp_path)

    assert result.failures == []
    assert any("citation metadata aligned" in item for item in result.passed)


def test_citation_metadata_rejects_version_mismatch(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    data = json.loads((tmp_path / "codemeta.json").read_text(encoding="utf-8"))
    data["softwareVersion"] = "0.2.0"
    (tmp_path / "codemeta.json").write_text(json.dumps(data), encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("version mismatch" in failure for failure in result.failures)


def test_citation_metadata_rejects_title_mismatch(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    data = json.loads((tmp_path / ".zenodo.json").read_text(encoding="utf-8"))
    data["title"] = "Wrong Title"
    (tmp_path / ".zenodo.json").write_text(json.dumps(data), encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("title mismatch" in failure for failure in result.failures)


def test_citation_metadata_rejects_missing_discovery_keyword(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    data = json.loads((tmp_path / "codemeta.json").read_text(encoding="utf-8"))
    data["keywords"] = ["bioinformatics"]
    (tmp_path / "codemeta.json").write_text(json.dumps(data), encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("codemeta.json missing discovery keyword" in failure for failure in result.failures)


def test_citation_metadata_rejects_weak_pyproject_discovery_metadata(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    (tmp_path / "pyproject.toml").write_text(
        "[project]\n"
        'name = "dotmatch"\n'
        'version = "0.1.0"\n'
        'description = "Fast matching"\n'
        'license = "Apache-2.0"\n'
        'keywords = ["bioinformatics"]\n'
        "classifiers = []\n"
        "[project.urls]\n"
        'Homepage = "https://github.com/dnncha/dotmatch"\n',
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("pyproject.toml missing discovery keyword: CRISPR" in failure for failure in result.failures)
    assert any("pyproject.toml description must mention known-target short-DNA assignment" in failure for failure in result.failures)
    assert any("pyproject.toml classifiers must include Topic :: Scientific/Engineering :: Bio-Informatics" in failure for failure in result.failures)
    assert any("pyproject.toml project URLs must include Repository" in failure for failure in result.failures)


def test_citation_metadata_rejects_unminted_doi_field(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    citation = (tmp_path / "CITATION.cff").read_text(encoding="utf-8") + "doi: 10.5281/zenodo.123\n"
    (tmp_path / "CITATION.cff").write_text(citation, encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("has a DOI field before Zenodo release" in failure for failure in result.failures)
