#!/usr/bin/env python3
"""Fail unless CRISPR SOTA evidence has real datasets, repeats, competitors, and oracle validation."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks" / "raw"
DATASETS = ["mageck_yusa", "sanson_brunello"]
FULL_FASTQ_MIN_READS = {
    "mageck_yusa": 20394663,
    "sanson_brunello": 246950411,
}
FULL_FASTQ_SAMPLE_READS = {
    "mageck_yusa": {
        "plasmid": 10093905,
        "ESC1": 10300758,
    },
    "sanson_brunello": {
        "plasmid": 9821128,
        "RepA": 76471324,
        "RepB": 85301059,
        "RepC": 75356900,
    },
}


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"missing required artifact: {path}")
    with path.open() as fh:
        return list(csv.DictReader(fh))


def as_int(value: str | None, default: int = 0) -> int:
    if value is None or value == "":
        return default
    if value == "full":
        return 10**18
    return int(float(value))


def as_float(value: str | None, default: float = 0.0) -> float:
    if value is None or value == "":
        return default
    return float(value)


def as_read_count(value: str | None) -> int:
    if value is None or value == "":
        return 0
    try:
        return int(float(value))
    except ValueError:
        return 0


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def full_rows_for(ok: list[dict[str, str]], dataset: str, tool: str) -> list[dict[str, str]]:
    return [
        row for row in ok
        if row.get("dataset_id") == dataset
        and row.get("tool") == tool
        and row.get("requested_records_per_sample") == "full"
    ]


def has_full_depth_evidence(rows: list[dict[str, str]], dataset: str) -> bool:
    if any(as_read_count(row.get("n_reads")) >= FULL_FASTQ_MIN_READS[dataset] for row in rows):
        return True

    expected = FULL_FASTQ_SAMPLE_READS.get(dataset, {})
    if not expected:
        return False
    best_by_sample: dict[str, int] = {}
    for row in rows:
        sample_id = row.get("sample_id", "")
        if sample_id not in expected:
            continue
        best_by_sample[sample_id] = max(best_by_sample.get(sample_id, 0), as_read_count(row.get("n_reads")))
    return all(best_by_sample.get(sample_id, 0) >= reads for sample_id, reads in expected.items())


def full_depth_rate(rows: list[dict[str, str]], dataset: str) -> float | None:
    monolithic = [
        as_float(row.get("reads_per_sec"))
        for row in rows
        if as_read_count(row.get("n_reads")) >= FULL_FASTQ_MIN_READS[dataset]
        and as_float(row.get("reads_per_sec")) > 0
    ]
    if monolithic:
        return sum(monolithic) / len(monolithic)

    expected = FULL_FASTQ_SAMPLE_READS.get(dataset, {})
    if not expected:
        return None
    best_by_sample: dict[str, dict[str, str]] = {}
    for row in rows:
        sample_id = row.get("sample_id", "")
        if sample_id not in expected:
            continue
        if as_read_count(row.get("n_reads")) < expected[sample_id]:
            continue
        current = best_by_sample.get(sample_id)
        if current is None or as_float(row.get("reads_per_sec")) > as_float(current.get("reads_per_sec")):
            best_by_sample[sample_id] = row
    if set(best_by_sample) != set(expected):
        return None
    total_reads = sum(as_read_count(row.get("n_reads")) for row in best_by_sample.values())
    total_seconds = sum(as_float(row.get("seconds")) for row in best_by_sample.values())
    return total_reads / total_seconds if total_seconds > 0 else None


def repeated_gate(rows: list[dict[str, str]], min_records: int, min_repeats: int,
                  require_full: bool, require_mageck: bool, require_guide_counter: bool, failures: list[str],
                  min_guide_counter_speedup: float = 1.0) -> None:
    ok = [r for r in rows if r.get("exit_code") == "0"]
    require(bool(ok), "crispr_sota_repeated.csv has no successful rows", failures)
    required_tools = ["dotmatch_exact_k0", "dotmatch_hamming_k1", "dotmatch_levenshtein_k1"]
    if require_mageck:
        required_tools.append("mageck_count_exact")
    if require_guide_counter:
        required_tools.append("guide_counter_one_mismatch")

    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in ok:
        grouped[(row.get("dataset_id", ""), row.get("tool", ""), row.get("requested_records_per_sample", ""))].append(row)

    for dataset in DATASETS:
        require(any(r.get("dataset_id") == dataset for r in ok), f"missing successful rows for {dataset}", failures)
        for tool in required_tools:
            best = max(
                (
                    len(group)
                    for (ds, t, requested), group in grouped.items()
                    if ds == dataset and t == tool and as_int(requested) >= min_records
                ),
                default=0,
            )
            require(best >= min_repeats,
                    f"{dataset}:{tool} needs >= {min_repeats} repeats at >= {min_records} records/sample; found {best}",
                    failures)
        if require_full:
            for tool in ["dotmatch_exact_k0", "dotmatch_hamming_k1", "dotmatch_levenshtein_k1"]:
                full_group = full_rows_for(ok, dataset, tool)
                require(has_full_depth_evidence(full_group, dataset),
                        f"{dataset}:{tool} needs at least one full FASTQ timing row", failures)
                if full_group and not has_full_depth_evidence(full_group, dataset):
                    max_reads = max(as_read_count(row.get("n_reads")) for row in full_group)
                    min_full_reads = FULL_FASTQ_MIN_READS[dataset]
                    require(max_reads >= min_full_reads,
                            f"{dataset}:{tool} full FASTQ row has too few reads: "
                            f"max_n_reads={max_reads}; need >= {min_full_reads}",
                            failures)

    for row in ok:
        if row.get("tool") == "dotmatch_levenshtein_k1" and row.get("verified_per_read"):
            verified_per_read = as_float(row.get("verified_per_read"))
            n_targets = as_int(row.get("n_targets"))
            collapse_limit = max(5.0, 0.001 * float(n_targets))
            require(verified_per_read <= collapse_limit,
                    f"{row.get('dataset_id')} Levenshtein candidate collapse is too weak: "
                    f"verified_per_read={verified_per_read:.4f}, n_targets={n_targets}, "
                    f"limit={collapse_limit:.4f}",
                    failures)

    if require_guide_counter:
        for dataset in DATASETS:
            dotmatch = [
                as_float(r.get("reads_per_sec"))
                for r in ok
                if r.get("dataset_id") == dataset
                and r.get("tool") == "dotmatch_hamming_k1"
                and r.get("requested_records_per_sample") != "full"
                and as_int(r.get("requested_records_per_sample")) >= min_records
            ]
            guide_counter = [
                as_float(r.get("reads_per_sec"))
                for r in ok
                if r.get("dataset_id") == dataset
                and r.get("tool") == "guide_counter_one_mismatch"
                and r.get("requested_records_per_sample") != "full"
                and as_int(r.get("requested_records_per_sample")) >= min_records
            ]
            if not dotmatch or not guide_counter:
                require(False, f"{dataset} missing rows for DotMatch-vs-guide-counter speedup gate", failures)
                continue
            dm_mean = sum(dotmatch) / len(dotmatch)
            gc_mean = sum(guide_counter) / len(guide_counter)
            speedup = dm_mean / gc_mean if gc_mean > 0 else 0.0
            require(speedup >= min_guide_counter_speedup,
                    f"{dataset} DotMatch Hamming mean speedup vs guide-counter is {speedup:.2f}x; "
                    f"need >= {min_guide_counter_speedup:.2f}x",
                    failures)
            if require_full:
                full_dm_mean = full_depth_rate(full_rows_for(ok, dataset, "dotmatch_hamming_k1"), dataset)
                full_gc_mean = full_depth_rate(full_rows_for(ok, dataset, "guide_counter_one_mismatch"), dataset)
                if full_dm_mean is None or full_gc_mean is None:
                    require(False, f"{dataset} missing full rows for DotMatch-vs-guide-counter speedup gate", failures)
                    continue
                full_speedup = full_dm_mean / full_gc_mean if full_gc_mean > 0 else 0.0
                require(full_speedup >= min_guide_counter_speedup,
                        f"{dataset} full DotMatch Hamming speedup vs guide-counter is {full_speedup:.2f}x; "
                        f"need >= {min_guide_counter_speedup:.2f}x",
                        failures)


def validation_gate(rows: list[dict[str, str]], min_checked: int, failures: list[str]) -> None:
    require(bool(rows), "crispr_sota_edlib_validation.csv is empty", failures)
    seen = set()
    for row in rows:
        dataset = row.get("dataset", "")
        seen.add(dataset)
        require(as_int(row.get("mismatches")) == 0,
                f"Edlib mismatch for {dataset}:{row.get('sample', '')}", failures)
        require(as_int(row.get("checked_reads")) >= min_checked,
                f"Edlib checked_reads below {min_checked} for {dataset}:{row.get('sample', '')}: {row.get('checked_reads')}",
                failures)
    for dataset in DATASETS:
        require(dataset in seen, f"missing Edlib validation rows for {dataset}", failures)


def agreement_gate(rows: list[dict[str, str]], require_guide_counter: bool, failures: list[str],
                   min_count_total: int = 100000) -> None:
    by_dataset: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_dataset[row.get("dataset", "")].append(row)
    for dataset in DATASETS:
        rows_for_dataset = by_dataset.get(dataset, [])
        require(bool(rows_for_dataset), f"missing count agreement rows for {dataset}", failures)
        exact = [r for r in rows_for_dataset if r.get("comparison", "").endswith("dotmatch_exact_vs_mageck_exact")]
        require(bool(exact), f"missing exact DotMatch-vs-MAGeCK agreement row for {dataset}", failures)
        if dataset == "mageck_yusa":
            require(bool(exact) and exact[0].get("status") == "ok",
                    f"missing exact DotMatch-vs-MAGeCK agreement for {dataset}", failures)
        if dataset == "mageck_yusa" and exact and exact[0].get("status") == "ok":
            require(as_int(exact[0].get("total_delta")) == 0,
                    f"{dataset} exact total differs from MAGeCK", failures)
            require(as_int(exact[0].get("differing_guides")) == 0,
                    f"{dataset} exact guide-level differences vs MAGeCK", failures)
        if require_guide_counter:
            ham = [r for r in rows_for_dataset if r.get("comparison", "").endswith("dotmatch_hamming_vs_guide_counter")]
            require(bool(ham) and ham[0].get("status") == "ok",
                    f"missing Hamming DotMatch-vs-guide-counter agreement for {dataset}", failures)
            if ham and ham[0].get("status") == "ok":
                left_total = as_int(ham[0].get("total_left"))
                right_total = as_int(ham[0].get("total_right"))
                require(
                    min(left_total, right_total) >= min_count_total,
                    f"{dataset} Hamming count agreement is below publication depth: "
                    f"left={left_total}, right={right_total}; need both >= {min_count_total}",
                    failures,
                )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeated", default=str(RAW / "crispr_sota_repeated.csv"))
    parser.add_argument("--validation", default=str(RAW / "crispr_sota_edlib_validation.csv"))
    parser.add_argument("--count-agreement", default=str(RAW / "crispr_sota_count_agreement_summary.csv"))
    parser.add_argument("--min-records", type=int, default=100000)
    parser.add_argument("--min-repeats", type=int, default=5)
    parser.add_argument("--min-edlib-checked", type=int, default=10000)
    parser.add_argument("--require-full", action="store_true", default=True)
    parser.add_argument("--no-full", action="store_false", dest="require_full")
    parser.add_argument("--require-mageck", action="store_true", default=True)
    parser.add_argument("--no-mageck", action="store_false", dest="require_mageck")
    parser.add_argument("--require-guide-counter", action="store_true", default=True)
    parser.add_argument("--no-guide-counter", action="store_false", dest="require_guide_counter")
    parser.add_argument("--min-guide-counter-speedup", type=float, default=1.0,
                        help="minimum mean DotMatch Hamming speedup over guide-counter for strict speed claim")
    parser.add_argument("--skip-count-agreement", action="store_true")
    parser.add_argument("--smoke", action="store_true", help="lower thresholds for local graph plumbing only")
    args = parser.parse_args()
    if args.smoke:
        args.min_records = 1
        args.min_repeats = 1
        args.min_edlib_checked = 1
        args.require_full = False
        args.require_mageck = False
        args.require_guide_counter = False
        args.skip_count_agreement = True

    failures: list[str] = []
    repeated_gate(read_rows(Path(args.repeated)), args.min_records, args.min_repeats,
                  args.require_full, args.require_mageck, args.require_guide_counter, failures,
                  args.min_guide_counter_speedup)
    validation_gate(read_rows(Path(args.validation)), args.min_edlib_checked, failures)
    if not args.skip_count_agreement:
        agreement_gate(read_rows(Path(args.count_agreement)), args.require_guide_counter, failures, args.min_records)

    if failures:
        print("CRISPR SOTA GATE: FAIL")
        for failure in failures:
            print(f"- {failure}")
        raise SystemExit(1)
    print("CRISPR SOTA GATE: PASS")


if __name__ == "__main__":
    main()
