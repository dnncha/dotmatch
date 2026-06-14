#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


DEFAULT_NATIVE_CSV = Path("benchmarks/raw/native_edlib_assignment.csv")
DEFAULT_NATIVE_REPORT = Path("docs/benchmarks/native/README.md")
KEY_FIELDS = (
    "workload",
    "error_mode",
    "n_reads",
    "n_targets",
    "len",
    "k",
    "err",
    "indel_rate",
    "repeat",
)


class AuditResult:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def _fnum(row: dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, "") or "0")
    except ValueError:
        return 0.0


def _key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row.get(field, "") for field in KEY_FIELDS)


def _best_by_key(rows: list[dict[str, str]]) -> dict[tuple[str, ...], dict[str, str]]:
    best: dict[tuple[str, ...], dict[str, str]] = {}
    for row in rows:
        key = _key(row)
        current = best.get(key)
        if current is None or _fnum(row, "reads_per_sec") > _fnum(current, "reads_per_sec"):
            best[key] = row
    return best


def _int(row: dict[str, str], key: str) -> int:
    try:
        return int(row.get(key, "") or "0")
    except ValueError:
        return 0


def _metric_row(claim: str, ratios: list[float], verified: list[float], min_required: float, max_verified: float) -> dict[str, str] | None:
    if not ratios or not verified:
        return None
    return {
        "claim": claim,
        "large_library_rows": str(len(ratios)),
        "min_speedup_vs_edlib": f"{min(ratios):.2f}",
        "median_speedup_vs_edlib": f"{statistics.median(ratios):.2f}",
        "max_verified_per_read": f"{max(verified):.2f}",
        "min_speedup_required": f"{min_required:.2f}",
        "max_verified_required": f"{max_verified:.2f}",
    }


def report_gate(report: Path, gated_rows: list[dict[str, str]], result: AuditResult) -> None:
    if not report.is_file():
        result.failures.append(f"missing native benchmark report: {report.as_posix()}")
        return
    text = report.read_text(encoding="utf-8")
    for required in [
        "## Gated Native Scaling Claims",
        "not end-to-end workflow speed claims",
        "This remains scoped to packed A/C/G/T fixed-window assignment up to 32 bases",
    ]:
        if required not in text:
            result.failures.append(f"native benchmark report must mention: {required}")
    for row in gated_rows:
        expected = (
            f"| {row['claim']} | {row['large_library_rows']} | {row['min_speedup_vs_edlib']} | "
            f"{row['median_speedup_vs_edlib']} | {row['max_verified_per_read']} | "
            f"{row['min_speedup_required']} | {row['max_verified_required']} |"
        )
        if expected not in text:
            result.failures.append(f"native benchmark report missing gated scaling row: {expected}")


def audit(path: Path, report: Path | None = None) -> AuditResult:
    result = AuditResult()
    if not path.is_file():
        result.failures.append(f"missing native exact benchmark CSV: {path.as_posix()}")
        return result

    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        result.failures.append(f"empty native exact benchmark CSV: {path.as_posix()}")
        return result

    direct = [row for row in rows if row.get("tool") in {"dotmatch_exact_direct", "dotmatch_exact_batch"} and row.get("k") == "0"]
    exact_hash = [row for row in rows if row.get("tool") == "exact_hash_lookup" and row.get("k") == "0"]
    edlib = [row for row in rows if row.get("tool") == "edlib_native_scan" and row.get("k") == "0"]
    if not direct:
        result.failures.append("native exact benchmark must include dotmatch_exact_direct or dotmatch_exact_batch k=0 rows")
    if not exact_hash:
        result.failures.append("native exact benchmark must include exact_hash_lookup k=0 rows")
    if not edlib:
        result.failures.append("native exact benchmark must include edlib_native_scan k=0 rows")

    for row in direct + exact_hash + edlib:
        if row.get("mismatches") not in {"0", "0.0", ""}:
            result.failures.append(f"{row.get('tool')} row has assignment mismatches for {_key(row)}")
        if _fnum(row, "reads_per_sec") <= 0.0:
            result.failures.append(f"{row.get('tool')} row has non-positive throughput for {_key(row)}")

    hash_by_key = {_key(row): row for row in exact_hash}
    direct_by_key = _best_by_key(direct)
    ratios: list[float] = []
    large_library_ratios: list[float] = []
    missing_pairs = 0
    for row in direct_by_key.values():
        key = _key(row)
        hash_row = hash_by_key.get(key)
        if hash_row is None:
            missing_pairs += 1
            continue
        hash_rps = _fnum(hash_row, "reads_per_sec")
        direct_rps = _fnum(row, "reads_per_sec")
        if hash_rps > 0.0 and direct_rps > 0.0:
            ratio = direct_rps / hash_rps
            ratios.append(ratio)
            n_targets = _int(row, "n_targets")
            if n_targets >= 4096:
                large_library_ratios.append(ratio)
    if direct and missing_pairs:
        result.failures.append(f"native exact benchmark is missing {missing_pairs} exact_hash_lookup comparison rows")
    if ratios:
        result.passed.append(f"dotmatch exact/exact_hash_lookup median ratio: {sorted(ratios)[len(ratios) // 2]:.3f}")
    if not large_library_ratios:
        result.failures.append("native exact benchmark must include n_targets>=4096 exact-hash comparison rows")
    elif min(large_library_ratios) <= 1.0:
        result.failures.append(
            f"large-library native exact rows must beat exact_hash_lookup; minimum ratio {min(large_library_ratios):.3f}"
        )
    else:
        result.passed.append(
            f"large-library dotmatch exact/exact_hash_lookup minimum ratio: {min(large_library_ratios):.3f}"
        )

    indexed_k1 = [
        row for row in rows
        if row.get("tool") == "dotmatch_indexed"
        and row.get("k") == "1"
        and row.get("error_mode") in {"exact", "one_substitution"}
    ]
    edlib_k1 = [row for row in rows if row.get("tool") == "edlib_native_scan" and row.get("k") == "1"]
    approximate = [
        row for row in rows
        if row.get("tool") in {"bk_tree", "neighbor_lookup"} and row.get("k") == "1"
    ]
    if not indexed_k1:
        result.failures.append("native assignment benchmark must include dotmatch_indexed k=1 rows")
    if not edlib_k1:
        result.failures.append("native assignment benchmark must include edlib_native_scan k=1 rows")
    if not approximate:
        result.failures.append("native assignment benchmark must include k=1 BK-tree or neighbor lookup baseline rows")

    edlib_by_key = {_key(row): row for row in edlib_k1}
    approximate_by_key: dict[tuple[str, ...], dict[str, str]] = {}
    for row in approximate:
        key = _key(row)
        current = approximate_by_key.get(key)
        if current is None or _fnum(row, "reads_per_sec") > _fnum(current, "reads_per_sec"):
            approximate_by_key[key] = row

    edlib_ratios: list[float] = []
    approximate_ratios: list[float] = []
    large_library_approximate_ratios: list[float] = []
    large_library_verified: list[float] = []
    large_library_k1_substitution_ratios: list[float] = []
    large_library_k1_substitution_verified: list[float] = []
    missing_edlib = 0
    missing_approximate = 0
    for row in indexed_k1:
        if row.get("mismatches") not in {"0", "0.0", ""}:
            result.failures.append(f"dotmatch_indexed k=1 row has assignment mismatches for {_key(row)}")
        if _fnum(row, "reads_per_sec") <= 0.0:
            result.failures.append(f"dotmatch_indexed k=1 row has non-positive throughput for {_key(row)}")
        edlib_row = edlib_by_key.get(_key(row))
        if edlib_row is None:
            missing_edlib += 1
        elif _fnum(edlib_row, "reads_per_sec") > 0.0:
            edlib_ratio = _fnum(row, "reads_per_sec") / _fnum(edlib_row, "reads_per_sec")
            edlib_ratios.append(edlib_ratio)
            if _int(row, "n_targets") >= 4096 and row.get("error_mode") == "one_substitution":
                large_library_k1_substitution_ratios.append(edlib_ratio)
        baseline_row = approximate_by_key.get(_key(row))
        if baseline_row is None:
            missing_approximate += 1
        elif _fnum(baseline_row, "reads_per_sec") > 0.0:
            ratio = _fnum(row, "reads_per_sec") / _fnum(baseline_row, "reads_per_sec")
            approximate_ratios.append(ratio)
            if _int(row, "n_targets") >= 4096:
                large_library_approximate_ratios.append(ratio)
        if _int(row, "n_targets") >= 4096:
            large_library_verified.append(_fnum(row, "verified_per_read"))
            if row.get("error_mode") == "one_substitution":
                large_library_k1_substitution_verified.append(_fnum(row, "verified_per_read"))

    if indexed_k1 and missing_edlib:
        result.failures.append(f"native k=1 indexed benchmark is missing {missing_edlib} edlib comparison rows")
    if indexed_k1 and missing_approximate:
        result.failures.append(
            f"native k=1 indexed benchmark is missing {missing_approximate} BK-tree/neighbor comparison rows"
        )
    if edlib_ratios:
        result.passed.append(f"dotmatch k=1/edlib median ratio: {sorted(edlib_ratios)[len(edlib_ratios) // 2]:.3f}")
        if min(edlib_ratios) <= 10.0:
            result.failures.append(
                f"native k=1 indexed rows must beat exhaustive Edlib scan by >10x; minimum ratio {min(edlib_ratios):.3f}"
            )
    if approximate_ratios:
        result.passed.append(
            f"dotmatch k=1/best approximate-baseline median ratio: {sorted(approximate_ratios)[len(approximate_ratios) // 2]:.3f}"
        )
    if large_library_approximate_ratios:
        result.passed.append(
            "large-library dotmatch k=1/best approximate-baseline minimum ratio: "
            f"{min(large_library_approximate_ratios):.3f}"
        )
        if min(large_library_approximate_ratios) <= 1.0:
            result.failures.append(
                "large-library native k=1 indexed rows must beat the best BK-tree/neighbor baseline; "
                f"minimum ratio {min(large_library_approximate_ratios):.3f}"
            )
    elif indexed_k1:
        result.failures.append("native k=1 benchmark must include n_targets>=4096 approximate-baseline comparisons")
    if not large_library_verified:
        result.failures.append("native k=1 benchmark must include n_targets>=4096 indexed rows")
    elif max(large_library_verified) > 1.05:
        result.failures.append(
            "large-library native k=1 indexed rows must verify no more than 1.05 candidates/read; "
            f"maximum {max(large_library_verified):.3f}"
        )
    else:
        result.passed.append(
            f"large-library k=1 verified candidates/read maximum: {max(large_library_verified):.3f}"
        )
    indexed_k2 = [
        row for row in rows
        if row.get("tool") == "dotmatch_indexed"
        and row.get("k") == "2"
        and row.get("error_mode") in {"one_insertion", "one_deletion"}
    ]
    edlib_k2 = [row for row in rows if row.get("tool") == "edlib_native_scan" and row.get("k") == "2"]
    if not indexed_k2:
        result.failures.append("native assignment benchmark must include dotmatch_indexed Levenshtein k=2 insertion/deletion rows")
    if not edlib_k2:
        result.failures.append("native assignment benchmark must include edlib_native_scan Levenshtein k=2 rows")

    edlib_k2_by_key = {_key(row): row for row in edlib_k2}
    k2_ratios: list[float] = []
    large_library_k2_ratios: list[float] = []
    large_library_k2_verified: list[float] = []
    k2_modes = {row.get("error_mode", "") for row in indexed_k2}
    missing_k2_edlib = 0
    for row in indexed_k2:
        if row.get("mismatches") not in {"0", "0.0", ""}:
            result.failures.append(f"dotmatch_indexed k=2 row has assignment mismatches for {_key(row)}")
        if _fnum(row, "reads_per_sec") <= 0.0:
            result.failures.append(f"dotmatch_indexed k=2 row has non-positive throughput for {_key(row)}")
        edlib_row = edlib_k2_by_key.get(_key(row))
        if edlib_row is None:
            missing_k2_edlib += 1
            continue
        if _fnum(edlib_row, "reads_per_sec") > 0.0:
            ratio = _fnum(row, "reads_per_sec") / _fnum(edlib_row, "reads_per_sec")
            k2_ratios.append(ratio)
            if _int(row, "n_targets") >= 4096:
                large_library_k2_ratios.append(ratio)
                large_library_k2_verified.append(_fnum(row, "verified_per_read"))

    if indexed_k2 and missing_k2_edlib:
        result.failures.append(f"native k=2 indexed benchmark is missing {missing_k2_edlib} edlib comparison rows")
    if indexed_k2 and not {"one_insertion", "one_deletion"}.issubset(k2_modes):
        result.failures.append("native k=2 benchmark must include both one_insertion and one_deletion rows")
    if k2_ratios:
        result.passed.append(f"dotmatch k=2/edlib median ratio: {sorted(k2_ratios)[len(k2_ratios) // 2]:.3f}")
    if large_library_k2_ratios:
        result.passed.append(
            f"large-library dotmatch k=2/edlib minimum ratio: {min(large_library_k2_ratios):.3f}"
        )
        if min(large_library_k2_ratios) <= 8.0:
            result.failures.append(
                f"large-library native k=2 indexed rows must beat exhaustive Edlib scan by >8x; minimum ratio {min(large_library_k2_ratios):.3f}"
            )
    elif indexed_k2:
        result.failures.append("native k=2 benchmark must include n_targets>=4096 insertion/deletion rows")
    if large_library_k2_verified:
        if max(large_library_k2_verified) > 25.0:
            result.failures.append(
                "large-library native k=2 indexed rows must verify no more than 25 candidates/read; "
                f"maximum {max(large_library_k2_verified):.3f}"
            )
        else:
            result.passed.append(
                f"large-library k=2 verified candidates/read maximum: {max(large_library_k2_verified):.3f}"
            )
    indexed_k2_substitution = [
        row for row in rows
        if row.get("tool") == "dotmatch_indexed"
        and row.get("k") == "2"
        and row.get("error_mode") == "one_substitution"
    ]
    if not indexed_k2_substitution:
        result.failures.append("native assignment benchmark must include dotmatch_indexed Hamming-style k=2 substitution rows")
    k2_substitution_ratios: list[float] = []
    large_library_k2_substitution_ratios: list[float] = []
    large_library_k2_substitution_verified: list[float] = []
    missing_k2_substitution_edlib = 0
    for row in indexed_k2_substitution:
        if row.get("mismatches") not in {"0", "0.0", ""}:
            result.failures.append(f"dotmatch_indexed k=2 substitution row has assignment mismatches for {_key(row)}")
        if _fnum(row, "reads_per_sec") <= 0.0:
            result.failures.append(f"dotmatch_indexed k=2 substitution row has non-positive throughput for {_key(row)}")
        edlib_row = edlib_k2_by_key.get(_key(row))
        if edlib_row is None:
            missing_k2_substitution_edlib += 1
            continue
        if _fnum(edlib_row, "reads_per_sec") > 0.0:
            ratio = _fnum(row, "reads_per_sec") / _fnum(edlib_row, "reads_per_sec")
            k2_substitution_ratios.append(ratio)
            if _int(row, "n_targets") >= 4096:
                large_library_k2_substitution_ratios.append(ratio)
                large_library_k2_substitution_verified.append(_fnum(row, "verified_per_read"))
    if indexed_k2_substitution and missing_k2_substitution_edlib:
        result.failures.append(
            f"native k=2 substitution benchmark is missing {missing_k2_substitution_edlib} edlib comparison rows"
        )
    if k2_substitution_ratios:
        result.passed.append(
            "dotmatch k=2 substitution/edlib median ratio: "
            f"{sorted(k2_substitution_ratios)[len(k2_substitution_ratios) // 2]:.3f}"
        )
    if large_library_k2_substitution_ratios:
        result.passed.append(
            "large-library dotmatch k=2 substitution/edlib minimum ratio: "
            f"{min(large_library_k2_substitution_ratios):.3f}"
        )
        if min(large_library_k2_substitution_ratios) <= 8.0:
            result.failures.append(
                "large-library native k=2 substitution rows must beat exhaustive Edlib scan by >8x; "
                f"minimum ratio {min(large_library_k2_substitution_ratios):.3f}"
            )
    elif indexed_k2_substitution:
        result.failures.append("native k=2 substitution benchmark must include n_targets>=4096 rows")
    if large_library_k2_substitution_verified:
        if max(large_library_k2_substitution_verified) > 1.05:
            result.failures.append(
                "large-library native k=2 substitution rows must verify no more than 1.05 candidates/read; "
                f"maximum {max(large_library_k2_substitution_verified):.3f}"
            )
        else:
            result.passed.append(
                "large-library k=2 substitution verified candidates/read maximum: "
                f"{max(large_library_k2_substitution_verified):.3f}"
            )
    gated_report_rows = [
        row for row in [
            _metric_row(
                "k=1 substitution indexed rows",
                large_library_k1_substitution_ratios,
                large_library_k1_substitution_verified,
                10.0,
                1.05,
            ),
            _metric_row(
                "k=2 substitution indexed rows",
                large_library_k2_substitution_ratios,
                large_library_k2_substitution_verified,
                8.0,
                1.05,
            ),
            _metric_row(
                "Levenshtein k=2 insertion/deletion rows",
                large_library_k2_ratios,
                large_library_k2_verified,
                8.0,
                25.0,
            ),
        ]
        if row is not None
    ]
    if report is not None:
        report_gate(report, gated_report_rows, result)
    if not result.failures:
        result.passed.append("native exact, k=1 indexed, k=2 substitution, and k=2 indel benchmark rows are comparator-backed")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Check native exact-assignment benchmark evidence.")
    parser.add_argument("csv", nargs="?", default=str(DEFAULT_NATIVE_CSV), help="native benchmark CSV")
    parser.add_argument("--report", default=str(DEFAULT_NATIVE_REPORT), help="native benchmark report")
    parser.add_argument("--skip-report", action="store_true")
    args = parser.parse_args()

    result = audit(Path(args.csv), None if args.skip_report else Path(args.report))
    for item in result.passed:
        print(f"PASS: {item}")
    for item in result.failures:
        print(f"FAIL: {item}")
    if result.ok:
        print("NATIVE EXACT GATE: PASS")
        return 0
    print("NATIVE EXACT GATE: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
