#!/usr/bin/env python3

import argparse
import json
import re
from pathlib import Path


REQUIRED_CONTROLS = {
    "input_integrity",
    "memory_safety",
    "oracle_validation",
    "public_assay_evidence",
    "public_claim_wording",
    "distribution_reproducibility",
    "release_governance",
}

PUBLIC_CLAIM_FILES = [
    "README.md",
    "pyproject.toml",
    "codemeta.json",
    ".zenodo.json",
    "packaging/bioconda/meta.yaml",
    "examples/workflows/nf-core/modules/local/dotmatch/crispr_count/meta.yml",
    "examples/workflows/nf-core/modules/local/dotmatch/assay_run/meta.yml",
]

FORBIDDEN_PUBLIC_CLAIM_PATTERNS = [
    (
        re.compile(r"\bfast exact short-dna known-target assignment\b", re.IGNORECASE),
        "use deterministic, evidence-backed assignment language instead of broad fast-exact positioning",
    ),
    (
        re.compile(r"\b(?:sota|state[- ]of[- ]the[- ]art)\b", re.IGNORECASE),
        "SOTA language requires a scoped benchmark report and must not appear in broad public metadata",
    ),
]

SCIENTIFIC_CLAIMS_REQUIRED_FRAGMENTS = [
    "## Strongest Scoped Performance Evidence",
    "global short-read",
    "all demultiplexing tasks",
    "`make native-exact-gate`",
    "large-library fixed-length `k=2` substitution rows",
    "Levenshtein `k=2` insertion/deletion rows",
    "`make crispr-comparison-gate` requires two real public",
    "DotMatch-vs-Bowtie 1 Hamming `k=2`/`k=3`",
    "fair guide-counter compatibility lane is Hamming `k=1`, no indels",
    "calibrated posterior or likelihood-based assignment is not currently claimed",
    "not a target posterior probability",
    "`make public-crispr-evidence-gate`",
    "not universal CRISPR superiority",
    "`make barcode-comparison-gate`",
    "Levenshtein one-edit barcode lane remains synthetic fixture evidence",
]

README_REQUIRED_EVIDENCE_FRAGMENTS = [
    "Evidence boundary:",
    "DotMatch Evidence Notes",
    "strongest",
    "native fixed-window indexed assignment",
    "public CRISPR",
    "guide-counting comparisons",
    "checked public inline-barcode lanes",
    "broader",
    "BCL replacement claims need their",
]

DOCS_INDEX_REQUIRED_EVIDENCE_FRAGMENTS = [
    "Evidence boundary:",
    "DotMatch Evidence Notes",
    "strongest",
    "native fixed-window indexed assignment",
    "public CRISPR",
    "guide-counting comparisons",
    "checked public inline-barcode lanes",
    "broader",
    "BCL replacement claims",
]


class AuditResult:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def _make_targets(root: Path) -> set[str]:
    text = (root / "Makefile").read_text(encoding="utf-8")
    targets: set[str] = set()
    for line in text.splitlines():
        match = re.match(r"^([A-Za-z0-9_.-]+):(?:\s|$)", line)
        if match:
            targets.add(match.group(1))
    return targets


def _check_path(root: Path, field: str, value: str, result: AuditResult) -> None:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        result.failures.append(f"{field} must be repository-relative: {value}")
        return
    if not (root / path).exists():
        result.failures.append(f"missing {field}: {value}")


def _check_gate(gate: str, make_targets: set[str], result: AuditResult) -> None:
    parts = gate.split()
    if len(parts) != 2 or parts[0] != "make":
        result.failures.append(f"gate must be a simple make target command: {gate}")
        return
    if parts[1] not in make_targets:
        result.failures.append(f"missing make target for scientific readiness gate: {parts[1]}")


def _check_public_claim_wording(root: Path, result: AuditResult) -> None:
    for rel_path in PUBLIC_CLAIM_FILES:
        path = root / rel_path
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        for pattern, message in FORBIDDEN_PUBLIC_CLAIM_PATTERNS:
            for match in pattern.finditer(text):
                line_no = text.count("\n", 0, match.start()) + 1
                result.failures.append(f"{rel_path}:{line_no}: {message}")


def _check_scientific_claims_contract(root: Path, result: AuditResult) -> None:
    path = root / "docs" / "scientific-claims.md"
    if not path.exists():
        result.failures.append("missing scientific claims document: docs/scientific-claims.md")
        return
    text = path.read_text(encoding="utf-8")
    for fragment in SCIENTIFIC_CLAIMS_REQUIRED_FRAGMENTS:
        if fragment not in text:
            result.failures.append(f"docs/scientific-claims.md must retain scoped evidence statement: {fragment}")
    _check_native_claim_values(root, text, result)
    _check_crispr_claim_values(root, text, result)


def _parse_native_gated_rows(report_text: str) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    in_table = False
    fields: list[str] = []
    for line in report_text.splitlines():
        stripped = line.strip()
        if stripped == "## Gated Native Scaling Claims":
            in_table = True
            continue
        if in_table and stripped.startswith("## "):
            break
        if not in_table or not stripped.startswith("|"):
            continue
        cells = [cell.strip() for cell in stripped.strip("|").split("|")]
        if cells and all(set(cell) <= {"-"} for cell in cells):
            continue
        if cells and cells[0] == "claim":
            fields = cells
            continue
        if fields and len(cells) == len(fields):
            row = dict(zip(fields, cells))
            claim = row.get("claim", "")
            if claim:
                rows[claim] = row
    return rows


def _check_native_claim_values(root: Path, scientific_text: str, result: AuditResult) -> None:
    report = root / "docs" / "benchmarks" / "native" / "README.md"
    if not report.exists():
        result.failures.append("missing native benchmark report: docs/benchmarks/native/README.md")
        return
    rows = _parse_native_gated_rows(report.read_text(encoding="utf-8"))
    for claim in [
        "k=2 substitution indexed rows",
        "Levenshtein k=2 insertion/deletion rows",
    ]:
        row = rows.get(claim)
        if not row:
            result.failures.append(f"native benchmark report missing gated claim row: {claim}")
            continue
        speedup = row.get("min_speedup_vs_edlib", "")
        verified = row.get("max_verified_per_read", "")
        if speedup and f"minimum `{speedup}x` speedup over Edlib" not in scientific_text:
            result.failures.append(
                "docs/scientific-claims.md must cite current native benchmark minimum "
                f"for {claim}: {speedup}x"
            )
        if verified and f"at most `{verified}` verified" not in scientific_text:
            result.failures.append(
                "docs/scientific-claims.md must cite current native benchmark verified-candidate "
                f"bound for {claim}: {verified}"
            )


def _parse_crispr_hamming_rows(report_text: str) -> dict[str, str]:
    rows: dict[str, str] = {}
    in_table = False
    fields: list[str] = []
    for line in report_text.splitlines():
        stripped = line.strip()
        if stripped == "## Hamming k2/k3 External Comparator Rows":
            in_table = True
            continue
        if in_table and stripped.startswith("## "):
            break
        if not in_table or not stripped.startswith("|"):
            continue
        cells = [cell.strip() for cell in stripped.strip("|").split("|")]
        if cells and all(set(cell) <= {"-"} for cell in cells):
            continue
        if cells and cells[0] == "dataset":
            fields = cells
            continue
        if fields and len(cells) == len(fields):
            row = dict(zip(fields, cells))
            k = row.get("k", "")
            speedup = row.get("speedup", "")
            if k and speedup:
                rows[k] = speedup
    return rows


def _check_crispr_claim_values(root: Path, scientific_text: str, result: AuditResult) -> None:
    report = root / "docs" / "benchmarks" / "crispr_comparison" / "README.md"
    if not report.exists():
        result.failures.append("missing CRISPR comparison report: docs/benchmarks/crispr_comparison/README.md")
        return
    normalized_scientific_text = re.sub(r"\s+", " ", scientific_text)
    rows = _parse_crispr_hamming_rows(report.read_text(encoding="utf-8"))
    for k in ["2", "3"]:
        speedup = rows.get(k)
        if not speedup:
            result.failures.append(f"CRISPR comparison report missing Hamming k={k} observed speedup")
            continue
        expected = f"current Bowtie 1 artifact records `{speedup}x` for Hamming `k={k}`"
        if expected not in normalized_scientific_text:
            result.failures.append(
                "docs/scientific-claims.md must cite current CRISPR Bowtie 1 observed speedup "
                f"for Hamming k={k}: {speedup}x"
            )


def _check_readme_evidence_boundary(root: Path, result: AuditResult) -> None:
    path = root / "README.md"
    if not path.exists():
        result.failures.append("missing README.md")
        return
    text = path.read_text(encoding="utf-8")
    for fragment in README_REQUIRED_EVIDENCE_FRAGMENTS:
        if fragment not in text:
            result.failures.append(f"README.md must retain evidence-boundary statement: {fragment}")


def _check_docs_index_evidence_boundary(root: Path, result: AuditResult) -> None:
    path = root / "docs" / "index.md"
    if not path.exists():
        result.failures.append("missing docs/index.md")
        return
    text = path.read_text(encoding="utf-8")
    for fragment in DOCS_INDEX_REQUIRED_EVIDENCE_FRAGMENTS:
        if fragment not in text:
            result.failures.append(f"docs/index.md must retain evidence-boundary statement: {fragment}")


def audit(root: Path) -> AuditResult:
    root = root.resolve()
    result = AuditResult()
    path = root / "docs" / "scientific-readiness.json"
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        result.failures.append(f"docs/scientific-readiness.json could not be read: {exc}")
        return result

    if manifest.get("schema_version") != 1:
        result.failures.append("docs/scientific-readiness.json must declare schema_version 1")
    if manifest.get("status") != "evidence_bounded":
        result.failures.append("scientific readiness status must be evidence_bounded")
    if not manifest.get("scope"):
        result.failures.append("scientific readiness manifest must declare scope")
    not_validated = manifest.get("not_validated_for")
    if not isinstance(not_validated, list) or not not_validated:
        result.failures.append("scientific readiness manifest must declare not_validated_for boundaries")

    controls = manifest.get("controls")
    if not isinstance(controls, list):
        result.failures.append("scientific readiness manifest must contain a controls list")
        return result

    make_targets = _make_targets(root)
    seen: set[str] = set()
    for control in controls:
        if not isinstance(control, dict):
            result.failures.append("scientific readiness controls must be objects")
            continue
        control_id = str(control.get("id") or "")
        if not control_id:
            result.failures.append("scientific readiness control missing id")
            continue
        if control_id in seen:
            result.failures.append(f"duplicate scientific readiness control: {control_id}")
        seen.add(control_id)
        if control.get("status") != "required":
            result.failures.append(f"{control_id} status must be required")
        if not control.get("acceptance"):
            result.failures.append(f"{control_id} must declare acceptance")
        evidence = control.get("evidence")
        if not isinstance(evidence, list) or not evidence:
            result.failures.append(f"{control_id} must list evidence")
        else:
            for value in evidence:
                _check_path(root, f"{control_id} evidence", str(value), result)
        gates = control.get("gates")
        if not isinstance(gates, list) or not gates:
            result.failures.append(f"{control_id} must list gates")
        else:
            for gate in gates:
                _check_gate(str(gate), make_targets, result)

    missing = REQUIRED_CONTROLS - seen
    for control_id in sorted(missing):
        result.failures.append(f"missing required scientific readiness control: {control_id}")
    _check_public_claim_wording(root, result)
    _check_scientific_claims_contract(root, result)
    _check_readme_evidence_boundary(root, result)
    _check_docs_index_evidence_boundary(root, result)
    if not missing and not result.failures:
        result.passed.append("scientific readiness controls valid")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit DotMatch scientific readiness controls.")
    parser.add_argument("--root", default=".", help="repository root")
    args = parser.parse_args()

    result = audit(Path(args.root))
    for item in result.passed:
        print(f"PASS: {item}")
    for item in result.failures:
        print(f"FAIL: {item}")
    if result.ok:
        print("SCIENTIFIC READINESS: PASS")
        return 0
    print("SCIENTIFIC READINESS: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
