#!/usr/bin/env python3
"""Validate the GSE146194 Perturb-seq case-study contracts and artifacts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "examples" / "perturb_seq_gse146194"
PROTOCOL = BASE / "protocol.json"
EXPECTED = BASE / "expected-results.json"
RESULTS = BASE / "results.json"
PROVENANCE = BASE / "provenance.json"
FIXTURE_EXPECTED = BASE / "fixture" / "expected.json"
FIXTURE_RESULTS = BASE / "work" / "fixture" / "results.json"
CSV_PATH = ROOT / "benchmarks" / "raw" / "perturb_seq_gse146194.csv"
REPORT = ROOT / "docs" / "benchmarks" / "perturb_seq_gse146194" / "README.md"
HTML_REPORT = ROOT / "docs" / "benchmarks" / "perturb_seq_gse146194" / "report.html"


def read_json(path: Path, failures: list[str]) -> dict:
    if not path.is_file():
        failures.append(f"missing JSON artifact: {path.relative_to(ROOT)}")
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"invalid JSON artifact {path.relative_to(ROOT)}: {exc}")
        return {}
    if not isinstance(value, dict):
        failures.append(f"JSON artifact must contain an object: {path.relative_to(ROOT)}")
        return {}
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_equal(actual, expected, label: str, failures: list[str]) -> None:
    if actual != expected:
        failures.append(f"{label}: expected {expected!r}, observed {actual!r}")


def public_gate(failures: list[str], require_work: bool) -> None:
    protocol = read_json(PROTOCOL, failures)
    expected = read_json(EXPECTED, failures)
    results = read_json(RESULTS, failures)
    provenance = read_json(PROVENANCE, failures)
    if not all([protocol, expected, results, provenance]):
        return
    require_equal(
        protocol.get("protocol_status"),
        "frozen_before_read_inspection",
        "protocol status",
        failures,
    )
    require_equal(sha256(PROTOCOL), expected.get("protocol_sha256"), "protocol SHA-256", failures)
    require_equal(
        provenance.get("protocol_sha256"),
        expected.get("protocol_sha256"),
        "provenance protocol SHA-256",
        failures,
    )
    require_equal(
        results.get("case_study_id"),
        protocol.get("case_study_id"),
        "case-study identifier",
        failures,
    )
    accessions = results.get("dataset_accessions", {})
    for key, value in protocol.get("dataset", {}).get("accessions", {}).items():
        require_equal(accessions.get(key), value, f"dataset accession {key}", failures)

    access = results.get("access_and_reuse", {})
    require_equal(
        access.get("sequence_data", {}).get("policy_url"),
        "https://www.ncbi.nlm.nih.gov/home/about/policies/",
        "sequence-data policy URL",
        failures,
    )
    require_equal(
        access.get("publisher_supplement", {}).get("license_status"),
        "source terms apply; this workflow asserts no redistribution license",
        "publisher supplement license boundary",
        failures,
    )
    require_equal(
        access.get("repository_redistribution", {}).get("raw_reads"),
        False,
        "raw-read redistribution boundary",
        failures,
    )
    require_equal(
        provenance.get("access_and_reuse"),
        access,
        "provenance access and reuse record",
        failures,
    )

    guide = results.get("guide_library", {})
    expected_guide = expected.get("guide_library", {})
    for key in ["target_count", "workbook_sha256", "targets_sha256"]:
        require_equal(guide.get(key), expected_guide.get(key), f"guide library {key}", failures)
    require_equal(
        guide.get("audit", {}).get("minimum_pairwise_hamming_distance"),
        expected_guide.get("minimum_pairwise_hamming_distance"),
        "guide library minimum pairwise Hamming distance",
        failures,
    )

    input_slice = results.get("input_slice", {})
    for key, value in expected.get("input_slice", {}).items():
        require_equal(input_slice.get(key), value, f"input slice {key}", failures)
    require_equal(
        input_slice.get("full_archive_md5_status"),
        "registry_value_recorded_not_reverified",
        "full archive MD5 evidence status",
        failures,
    )

    window = results.get("window_discovery", {})
    for key, value in expected.get("window_discovery", {}).items():
        require_equal(window.get(key), value, f"window discovery {key}", failures)

    runs = {
        f"k{run.get('k')}": run
        for run in results.get("runs", [])
        if isinstance(run, dict)
    }
    for run_id, wanted in expected.get("runs", {}).items():
        observed = runs.get(run_id, {})
        dotmatch = observed.get("dotmatch", {})
        agreement = observed.get("agreement", {})
        for key, value in wanted.items():
            source = (
                agreement
                if key in {"checked_records", "validation_mismatches", "agreement_rate"}
                else dotmatch
            )
            require_equal(source.get(key), value, f"{run_id} {key}", failures)
        require_equal(
            observed.get("ambiguity_policy"),
            "radius",
            f"{run_id} ambiguity policy",
            failures,
        )
        require_equal(observed.get("metric"), "hamming", f"{run_id} metric", failures)
        comparator = observed.get("comparator", {})
        require_equal(
            comparator.get("total_reads"),
            dotmatch.get("total_reads"),
            f"{run_id} comparator total",
            failures,
        )
        require_equal(
            comparator.get("assigned_unique"),
            dotmatch.get("assigned_unique"),
            f"{run_id} comparator unique",
            failures,
        )
    if set(runs) != {"k0", "k1"}:
        failures.append("public results must contain exactly k0 and k1 runs")

    resources = results.get("resources", {})
    if float(resources.get("wall_seconds", 0)) <= 0:
        failures.append("resource record must contain positive wall_seconds")
    if float(resources.get("peak_rss_mib", 0)) <= 0:
        failures.append("resource record must contain positive peak_rss_mib")
    if not provenance.get("dotmatch_version"):
        failures.append("provenance must record the DotMatch version")
    if len(provenance.get("commands", [])) != 2:
        failures.append("provenance must record both public DotMatch commands")

    if not CSV_PATH.is_file():
        failures.append("missing public raw CSV artifact")
    else:
        with CSV_PATH.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        if {row.get("k") for row in rows} != {"0", "1"}:
            failures.append("public raw CSV must contain k=0 and k=1 rows")
        for row in rows:
            if (
                row.get("validation_mismatches") != "0"
                or row.get("agreement_rate") != "1.00000000"
            ):
                failures.append(
                    "public raw CSV rows must record zero mismatches and full agreement"
                )

    if not REPORT.is_file():
        failures.append(f"missing rendered report artifact: {REPORT.relative_to(ROOT)}")
    else:
        text = REPORT.read_text(encoding="utf-8").lower()
        for required in [
            "48,000",
            "zero held-out per-read differences",
            "the unmatched and ambiguous columns are part of the result",
            "does not prove guide-per-cell accuracy",
        ]:
            if required.lower() not in text:
                failures.append(
                    f"{REPORT.relative_to(ROOT)} is missing evidence text: {required}"
                )
    if not HTML_REPORT.is_file():
        failures.append(
            f"missing rendered report artifact: {HTML_REPORT.relative_to(ROOT)}"
        )
    else:
        html_text = HTML_REPORT.read_text(encoding="utf-8").lower()
        for required in ["48,000 held-out", "evidence boundary", "srr11214031"]:
            if required not in html_text:
                failures.append(
                    f"{HTML_REPORT.relative_to(ROOT)} is missing rendered evidence text: {required}"
                )
    readme = BASE / "README.md"
    readme_text = readme.read_text(encoding="utf-8") if readme.is_file() else ""
    report_text = REPORT.read_text(encoding="utf-8") if REPORT.is_file() else ""
    combined = (report_text + readme_text).lower()
    for required in [
        "does not support cell-barcode correction",
        "full archive md5",
        "not locally reverified",
        "asserts no redistribution license",
        "next step",
    ]:
        if required not in combined:
            failures.append(f"case-study documentation is missing boundary text: {required}")

    if require_work:
        work = BASE / "work"
        checks = {
            "targets.tsv": expected_guide.get("targets_sha256"),
            "evaluation.fastq.gz": expected.get("input_slice", {}).get(
                "evaluation_gzip_sha256"
            ),
        }
        for name, digest in checks.items():
            path = work / name
            if not path.is_file():
                failures.append(
                    f"missing regenerated work artifact: {path.relative_to(ROOT)}"
                )
            elif sha256(path) != digest:
                failures.append(
                    f"regenerated work artifact checksum differs: {path.relative_to(ROOT)}"
                )
        for k in (0, 1):
            for stem in [
                f"dotmatch.k{k}.assignments.tsv",
                f"oracle.k{k}.assignments.tsv",
            ]:
                if not (work / stem).is_file():
                    failures.append(
                        f"missing regenerated assignment artifact: {(work / stem).relative_to(ROOT)}"
                    )


def fixture_gate(failures: list[str]) -> None:
    expected = read_json(FIXTURE_EXPECTED, failures)
    results = read_json(FIXTURE_RESULTS, failures)
    if not expected or not results:
        return
    if results.get("failures"):
        failures.append(f"fixture run reported failures: {results['failures']}")
    runs = {
        f"k{run.get('k')}": run
        for run in results.get("runs", [])
        if isinstance(run, dict)
    }
    for run_id, wanted in expected.get("runs", {}).items():
        observed = runs.get(run_id, {})
        for key, value in wanted.items():
            source = (
                observed.get("agreement", {})
                if key == "validation_mismatches"
                else observed.get("dotmatch", {})
            )
            require_equal(source.get(key), value, f"fixture {run_id} {key}", failures)
    k1 = expected.get("runs", {}).get("k1", {})
    for key in ["assigned_corrected", "ambiguous", "unmatched", "invalid"]:
        if int(k1.get(key, 0)) < 1:
            failures.append(f"fixture expected contract must exercise {key}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-public", action="store_true")
    parser.add_argument("--require-work", action="store_true")
    parser.add_argument("--require-fixture-work", action="store_true")
    args = parser.parse_args(argv)
    failures: list[str] = []
    if args.require_public:
        public_gate(failures, args.require_work)
    if args.require_fixture_work:
        fixture_gate(failures)
    if not args.require_public and not args.require_fixture_work:
        public_gate(failures, False)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        print("PERTURB-SEQ GSE146194: FAIL")
        return 1
    print("PERTURB-SEQ GSE146194: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
