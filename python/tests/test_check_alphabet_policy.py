import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_alphabet_policy.py"
POLICY = "literal-byte; A/C/G/T/N/IUPAC symbols are ordinary byte symbols; no wildcard expansion"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_alphabet_policy", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["tool", "workflow", "alphabet_policy"]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_repo(root: Path) -> None:
    files = {
        "include/qdalign.h": f'#define QDALN_ALPHABET_POLICY "{POLICY}"\n',
        "src/qdalign.c": 'const char *qdaln_alphabet_policy(void) {\n    return QDALN_ALPHABET_POLICY;\n}\n',
        "README.md": f"# DotMatch\n\nN and IUPAC use `{POLICY}`.\n",
        "docs/schemas.md": "N and IUPAC ambiguity symbols are literal byte symbols, not wildcard expansions.\n",
        "docs/scientific-claims.md": "current `N`/IUPAC behavior is literal-byte matching, not wildcard expansion semantics.\n",
        "docs/methods-and-citation.md": "The current alphabet policy is literal-byte, with N and IUPAC not expanded as wildcards.\n",
    }
    for path, text in files.items():
        full = root / path
        full.parent.mkdir(parents=True, exist_ok=True)
        full.write_text(text, encoding="utf-8")
    _write_csv(
        root / "benchmarks" / "raw" / "assay.csv",
        [
            {"tool": "dotmatch_count", "workflow": "public_assay", "alphabet_policy": POLICY},
            {"tool": "exact_slice_hash", "workflow": "public_assay", "alphabet_policy": ""},
        ],
    )


def test_alphabet_policy_accepts_literal_byte_contract(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)

    result = checker.audit(tmp_path)

    assert result.failures == []
    assert any("alphabet policy contract documented" in item for item in result.passed)


def test_alphabet_policy_rejects_header_source_mismatch(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    (tmp_path / "src" / "qdalign.c").write_text(
        'const char *qdaln_alphabet_policy(void) {\n    return "wildcard";\n}\n',
        encoding="utf-8",
    )

    result = checker.audit(tmp_path)

    assert any("qdaln_alphabet_policy must return QDALN_ALPHABET_POLICY" in failure for failure in result.failures)


def test_alphabet_policy_rejects_stale_invalid_wording(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    stale = "ACGT only; N and IUPAC are " + "invalid for matching"
    (tmp_path / "docs" / "schemas.md").write_text(stale + "\n", encoding="utf-8")

    result = checker.audit(tmp_path)

    assert any("stale N/IUPAC policy wording" in failure for failure in result.failures)


def test_alphabet_policy_requires_dotmatch_rows_to_record_policy(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    _write_csv(
        tmp_path / "benchmarks" / "raw" / "assay.csv",
        [{"tool": "dotmatch_count", "workflow": "public_assay", "alphabet_policy": ""}],
    )

    result = checker.audit(tmp_path)

    assert any("dotmatch_count row in benchmarks/raw/assay.csv must record alphabet_policy" in failure for failure in result.failures)


def test_alphabet_policy_rejects_dotmatch_rows_with_wrong_policy(tmp_path):
    checker = _load_checker()
    _write_repo(tmp_path)
    _write_csv(
        tmp_path / "benchmarks" / "raw" / "assay.csv",
        [{"tool": "dotmatch_count", "workflow": "public_assay", "alphabet_policy": "wildcard"}],
    )

    result = checker.audit(tmp_path)

    assert any("must use literal-byte alphabet_policy" in failure for failure in result.failures)
