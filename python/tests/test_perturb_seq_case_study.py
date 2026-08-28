from __future__ import annotations

import copy
import importlib.util
import sys
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "run_perturb_seq_gse146194.py"
SPEC = importlib.util.spec_from_file_location("dotmatch_perturb_seq_case_study", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
case_study = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = case_study
SPEC.loader.exec_module(case_study)


def test_access_record_does_not_overstate_source_rights() -> None:
    record = case_study.access_and_reuse_record()
    assert record["sequence_data"]["policy_url"] == (
        "https://www.ncbi.nlm.nih.gov/home/about/policies/"
    )
    assert record["publisher_supplement"]["license_status"] == (
        "source terms apply; this workflow asserts no redistribution license"
    )
    assert record["repository_redistribution"]["raw_reads"] is False
    assert record["repository_redistribution"]["publisher_workbook"] is False


def test_fixture_oracle_exercises_all_outcomes() -> None:
    fixture = ROOT / "examples" / "perturb_seq_gse146194" / "fixture"
    targets = case_study.load_targets(fixture / "targets.tsv")
    rows = [
        case_study.oracle_row(record, targets, 0, 18, 1)
        for record in case_study.iter_fastq(fixture / "reads.fastq")
    ]
    summary = case_study.summarize_assignments(rows)
    assert summary["assigned_exact"] == 1
    assert summary["assigned_corrected"] == 1
    assert summary["ambiguous"] == 2
    assert summary["unmatched"] == 1
    assert summary["invalid"] == 1
    assert case_study.library_audit(targets)["minimum_pairwise_hamming_distance"] == 1


def test_window_discovery_uses_declared_tie_breaks() -> None:
    fixture = ROOT / "examples" / "perturb_seq_gse146194" / "fixture"
    targets = case_study.load_targets(fixture / "targets.tsv")
    records = [
        case_study.FastqRecord(
            target.target_id,
            f"@{target.target_id}",
            f"NNNNN{target.sequence}AAAA",
            "+",
            "I" * 27,
        )
        for target in targets
    ]
    protocol = copy.deepcopy(case_study.load_json(case_study.DEFAULT_PROTOCOL))
    discovery = protocol["analysis_plan"]["window_discovery"]
    discovery["minimum_distinct_exact_targets"] = 2
    discovery["minimum_exact_assignment_fraction"] = 0.5
    oriented, result = case_study.discover_window(records, targets, protocol)
    assert result["orientation"] == "forward"
    assert result["target_start"] == 5
    assert result["distinct_exact_targets"] == 4
    assert [target.sequence for target in oriented] == [target.sequence for target in targets]


def test_primary_workbook_table_extraction_contract(tmp_path: Path) -> None:
    workbook = tmp_path / "fixture.xlsx"
    strings = [
        "Supplementary Table 2",
        "sgRNA",
        "Target gene name",
        "Sequence of GBC (only used for GBC Perturb-seq)",
        "GuideA",
        "GENEA",
        "ACGTACGTACGTACGTAC",
        "GuideB",
        "GENEB",
        "TTTTCCCCAAAAGGGGTT",
    ]
    shared = "".join(f"<si><t>{value}</t></si>" for value in strings)
    rows = """
      <row r="1"><c r="A1" t="s"><v>0</v></c></row>
      <row r="2"><c r="A2" t="s"><v>1</v></c><c r="B2" t="s"><v>2</v></c><c r="C2" t="s"><v>3</v></c></row>
      <row r="3"><c r="A3" t="s"><v>4</v></c><c r="B3" t="s"><v>5</v></c><c r="C3" t="s"><v>6</v></c></row>
      <row r="4"><c r="A4" t="s"><v>7</v></c><c r="B4" t="s"><v>8</v></c><c r="C4" t="s"><v>9</v></c></row>
    """
    with zipfile.ZipFile(workbook, "w") as archive:
        archive.writestr(
            "xl/workbook.xml",
            """<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Supplementary Table 2" sheetId="1" r:id="rId1"/></sheets></workbook>""",
        )
        archive.writestr(
            "xl/_rels/workbook.xml.rels",
            """<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/></Relationships>""",
        )
        archive.writestr(
            "xl/sharedStrings.xml",
            f"<sst xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" count=\"{len(strings)}\" uniqueCount=\"{len(strings)}\">{shared}</sst>",
        )
        archive.writestr(
            "xl/worksheets/sheet1.xml",
            f"<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\"><sheetData>{rows}</sheetData></worksheet>",
        )
    protocol = {
        "inputs": {
            "guide_library": {
                "worksheet": "Supplementary Table 2",
                "id_column": "sgRNA",
                "sequence_column": "Sequence of GBC (only used for GBC Perturb-seq)",
                "target_length": 18,
            }
        },
        "dataset": {"guide_count_expected": 2},
    }
    output = tmp_path / "targets.tsv"
    targets = case_study.extract_targets(workbook, protocol, output)
    assert [(target.target_id, target.sequence, target.gene) for target in targets] == [
        ("GuideA", "ACGTACGTACGTACGTAC", "GENEA"),
        ("GuideB", "TTTTCCCCAAAAGGGGTT", "GENEB"),
    ]
    assert output.read_text(encoding="utf-8").startswith("target_id\ttarget_seq\tgene\n")
