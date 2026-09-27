import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "check_mageck_yusa_reference_rows.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("check_mageck_yusa_reference_rows", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_reference_spotcheck_requires_both_sample_counts_and_gene(tmp_path):
    checker = _load_script()
    reference = _write(
        tmp_path / "reference.csv",
        "guide_id,gene,plasmid,ESC1\nsg1,GENE1,13,32\nsg2,GENE2,94,108\n",
    )
    counts = _write(
        tmp_path / "counts.tsv",
        "sgRNA\tGene\tplasmid\tESC1\nsg1\tGENE1\t13\t32\nsg2\tGENE2\t94\t108\n",
    )

    details = checker.compare_rows(reference, counts)

    assert [row["status"] for row in details] == ["match", "match"]


def test_reference_spotcheck_reports_missing_or_changed_rows(tmp_path):
    checker = _load_script()
    reference = _write(
        tmp_path / "reference.csv",
        "guide_id,gene,plasmid,ESC1\nsg1,GENE1,13,32\nsg2,GENE2,94,108\n",
    )
    counts = _write(
        tmp_path / "counts.tsv",
        "sgRNA\tGene\tplasmid\tESC1\nsg1\tGENE1\t13\t31\n",
    )

    details = checker.compare_rows(reference, counts)

    assert [row["status"] for row in details] == ["mismatch", "mismatch"]


def test_reference_spotcheck_rejects_duplicate_guides(tmp_path):
    checker = _load_script()
    counts = _write(
        tmp_path / "counts.tsv",
        "sgRNA\tGene\tplasmid\tESC1\nsg1\tGENE1\t13\t32\nsg1\tGENE1\t13\t32\n",
    )

    try:
        checker.read_unique_rows(counts, "sgRNA")
    except ValueError as exc:
        assert "duplicate" in str(exc)
    else:
        raise AssertionError("duplicate guide should fail closed")
