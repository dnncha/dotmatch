#!/usr/bin/env python3
"""Regenerate and audit the paired Hamming seed-layout evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "benchmarks/raw"
REPORT = ROOT / "docs/benchmarks/hamming_seed_layout/README.md"


def load(name: str, *, cli: bool = False):
    path = RAW / f"{name}.csv"
    provenance = json.loads(path.with_suffix(".json").read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == provenance["csv_sha256"], "CSV/provenance hash mismatch"
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows, "empty evidence"
    groups = defaultdict(dict)
    for row in rows:
        key = (row["workload"], int(row["k"])) + ((row["policy"],) if cli else ())
        pair = (row["tool"], int(row["repeat"]))
        assert pair not in groups[key], "duplicate repeat"
        assert pair[0] in {"baseline", "candidate"}, "unexpected engine"
        assert math.isfinite(float(row["seconds"])) and float(row["seconds"]) > 0
        assert math.isclose(float(row["reads_per_second"]) * float(row["seconds"]), int(row["n_reads"]), rel_tol=1e-9)
        if cli:
            assert int(row["exit_code"]) == 0 and int(row["count_mismatches"]) == 0
        else:
            assert int(row["oracle_reads"]) > 0
            assert int(row["candidates_considered"]) == int(row["candidates_verified"])
            if row["tool"] == "candidate":
                assert int(row["oracle_mismatches"]) == 0, "candidate oracle failure"
            if row["workload"] != "dense_literal_reads_1024_8" or int(row["k"]) != 3:
                assert int(row["baseline_mismatches"]) == int(row["oracle_mismatches"]) == 0
        groups[key][pair] = row
    for group in groups.values():
        repeats = {repeat for _, repeat in group}
        assert len(repeats) >= 3, "need at least three paired repeats"
        assert len(group) == 2 * len(repeats), "missing engine/repeat"
        for repeat in repeats:
            before, after = group["baseline", repeat], group["candidate", repeat]
            fields = ["n_reads", "n_targets", "length"]
            fields += ["counts_sha256", "outcomes_sha256"] if cli else [
                "oracle_reads", "reads_sha256", "targets_sha256", "candidates_verified", "baseline_mismatches"
            ]
            if not cli and int(after["baseline_mismatches"]) == 0:
                fields.append("results_sha256")
            assert all(before[field] == after[field] for field in fields), "unpaired input/result"
    if not cli:
        assert {key[0] for key in groups} == {
            "random_2048_20", "random_65536_20", "random_65536_32", "dense_shared_seed_1024_8",
            "dense_literal_reads_1024_8", "yusa_library_simulated_reads"
        }
        assert all((name, k) in groups for name, _ in groups for k in range(4)), "missing radius"
        regression = provenance["unknown_read_regression"]
        assert regression["candidate"] == regression["oracle"] != regression["baseline"]
        literal = groups["dense_literal_reads_1024_8", 3]["candidate", 0]
        assert int(literal["baseline_mismatches"]) > 0 and int(literal["oracle_reads"]) == int(literal["n_reads"])
    else:
        assert {key[1:] for key in groups} == {(2, "best"), (2, "radius"), (3, "best"), (3, "radius")}
    return groups, provenance


def summary(group):
    repeats = sorted({repeat for _, repeat in group})
    before = [float(group["baseline", repeat]["seconds"]) for repeat in repeats]
    after = [float(group["candidate", repeat]["seconds"]) for repeat in repeats]
    ratios = [a / b for a, b in zip(before, after)]
    return (statistics.median(before), statistics.median(after), statistics.median(ratios), min(ratios), max(ratios))


def render() -> str:
    native, provenance = load("hamming_seed_layout")
    cli, cli_provenance = load("hamming_seed_layout_cli", cli=True)
    oracle_checks = sum(int(g["candidate", 0]["oracle_reads"]) for g in native.values())
    checked_results = sum(int(g["candidate", 0]["n_reads"]) for g in native.values())
    literal = native["dense_literal_reads_1024_8", 3]["candidate", 0]
    public = native["yusa_library_simulated_reads", 3]["candidate", 0]
    lines = [
        "# Hamming seed-layout improvement", "",
        "This change accelerates indexed fixed-window Hamming assignment and fixes literal unknown-read verification.",
        f"The baseline is DotMatch `{provenance['baseline_ref']}`. Both engines use the same targets, reads, radii and complete result contract.",
        "These are before/after measurements on one Linux runner, not a new competitor comparison or evidence of overall state of the art.", "",
        "## What changed", "",
        "Seed buckets now occupy contiguous ranges. Each record carries the packed target code and length, avoiding linked-list traversal and separate target-array loads.",
        "Each target is verified only in its first matching partition, so packed A/C/G/T read queries need no growing duplicate set or heap allocation.",
        "The 64-bit runner uses 32 bytes per seed record instead of 40; this reduces seed-record storage by 20%, not total process memory by 20%.", "",
        "Completeness follows from the pigeonhole principle: with at most k substitutions across k+1 disjoint partitions, at least one partition must match exactly.",
        "Choosing the first matching partition visits every candidate once. Duplicate target rows remain distinct; exact matches, best ties, second-best distances and radius match counts remain explicit.",
        "Long sequences and literal symbols retain their supported fallback paths. The index remains immutable during queries, and the public C/Python interfaces are unchanged.", "",
        "## Correctness correction", "",
        "The old unknown-read route generated Hamming candidates but verified them with Levenshtein distance. This could understate distances and create a false best-match tie.",
        "For read `ANATGGAT` and targets `AATCGGAT`, `AAACGGAT`, the correct Hamming distances are 3 and 2. The old engine reported an ambiguous distance-2 tie.",
        "The new engine reports target 1 uniquely under `best`, with best distance 2, second-best distance 3 and two targets inside radius 3. The default `radius` policy correctly remains ambiguous.",
        "The same regression covers literal `R` and lowercase `a`; none is a wildcard. Historical Hamming k=3 results involving non-ACGT reads should be regenerated before reuse.", "",
        f"Across {checked_results:,} read/radius results, the new engine passed {oracle_checks:,} independent exhaustive byte-oracle checks with zero mismatches.",
        f"All A/C/G/T read results matched the baseline in every field. The fully oracle-checked dense literal corpus corrected {literal['baseline_mismatches']} results, including {literal['baseline_status_changes']} status/unique-assignment changes; baseline oracle errors are retained in the CSV.",
        "Native tests also exhaust all A/C/G/T words through length 5 at k=1..3, including duplicate target rows; test every length 0..33, mixed target lengths and dense shared seeds; and assert zero query heap allocations on a dense packed workload.",
        "Python and full-count CLI regressions exercise both ambiguity policies. This validates assignment semantics, not biological source identity or gene-hit accuracy.", "",
        "## Complete counting commands", "",
        f"The public Yusa library contains {int(public['n_targets']):,} 19-base guide rows. FASTQs are simulated A/C/G/T reads: equal exact, one-, two-, three-substitution and random-decoy groups.",
        "Each pair runs a fresh single-thread subprocess, including library/index construction, gzip FASTQ parsing, counting and MAGeCK TSV/JSON output. Five paired repeats alternate engine order.",
        "All guide counts and read-outcome/candidate totals agree between engines and with the validated native batch results. File hashes are recorded; timing excludes input generation and output validation.", "",
        "| Radius k | Policy | Baseline seconds | New seconds | Paired median speedup (min–max) |",
        "| --- | --- | --- | --- | --- |",
    ]
    for (_, k, policy), group in sorted(cli.items()):
        before, after, ratio, low, high = summary(group)
        lines.append(f"| {k} | {policy} | {before:.3f} | {after:.3f} | {ratio:.2f}× ({low:.2f}–{high:.2f}) |")
    lines += ["", "## Native batch assignment", "",
              "Times include statistics collection and exclude index construction, Python marshalling and output checks. Each engine is warmed once before the five alternating paired repeats.",
              "Large libraries use 256 uniformly sampled reads per radius for the independent oracle. Both dense 1,024-target workloads exhaust every read. All target rows participate in each oracle scan.", "",
              "| Workload | k | Reads | Baseline reads/s | New reads/s | Paired median speedup (min–max) |",
              "| --- | --- | --- | --- | --- | --- |"]
    for (name, k), group in sorted(native.items()):
        if k == 0:
            continue
        before, after, ratio, low, high = summary(group)
        n = int(group["candidate", 0]["n_reads"])
        lines.append(f"| {name} | {k} | {n:,} | {n / before:,.0f} | {n / after:,.0f} | {ratio:.2f}× ({low:.2f}–{high:.2f}) |")
    controls = [summary(group)[2] for (name, k), group in native.items() if k == 0]
    lines += ["", f"Exact k=0 control medians ranged from {min(controls):.2f}× to {max(controls):.2f}×; no exact-kernel improvement is claimed.",
              "Index construction is measured once per engine/workload and recorded separately in the native CSV. It is included in every complete-count timing above.", "",
              "## Reproduce and audit", "",
              "Recorded build: Linux x86-64, `cc -O3 -std=c11 -Wall -Wextra -Wpedantic -Iinclude -mavx2` (plus `-fPIC` for the shared object), with identical flags for baseline and new engine.",
              f"Compiler: `{provenance['compiler']}`. Platform: `{provenance['platform']}`. NumPy oracle: `{provenance['numpy']}`.",
              "Download the public library with the existing example fetcher if it is absent. Its small FASTQ download is not used by these benchmarks.", "", "```bash",
              "python3 scripts/fetch_mageck_demo.py --out benchmarks/real/data --subsample 25",
              "mkdir -p benchmarks/work/hamming-seed-baseline",
              f"git archive {provenance['baseline_ref']} Makefile pyproject.toml include src | tar -x -C benchmarks/work/hamming-seed-baseline",
              "make -C benchmarks/work/hamming-seed-baseline shared dotmatch CC=cc",
              "make shared dotmatch CC=cc",
              "python3 scripts/bench_hamming_seed_layout.py \\",
              "  --baseline-lib benchmarks/work/hamming-seed-baseline/libdotmatch.so \\",
              "  --baseline-source benchmarks/work/hamming-seed-baseline/src/qdalign.c \\",
              f"  --baseline-ref {provenance['baseline_ref']} \\",
              "  --public-library benchmarks/real/data/yusa_library.csv",
              "python3 scripts/bench_hamming_count_layout.py \\",
              "  --baseline-bin benchmarks/work/hamming-seed-baseline/dotmatch \\",
              "  --library benchmarks/real/data/yusa_library.csv",
              "python3 scripts/report_hamming_seed_layout.py",
              "python3 scripts/report_hamming_seed_layout.py --check --check-source",
              "make test cli-test python-test python-package-test repository-ready",
              "make coverage", "```", "",
              "Raw native [CSV](https://github.com/dnncha/dotmatch/blob/main/benchmarks/raw/hamming_seed_layout.csv) and [provenance](https://github.com/dnncha/dotmatch/blob/main/benchmarks/raw/hamming_seed_layout.json); complete-count [CSV](https://github.com/dnncha/dotmatch/blob/main/benchmarks/raw/hamming_seed_layout_cli.csv) and [provenance](https://github.com/dnncha/dotmatch/blob/main/benchmarks/raw/hamming_seed_layout_cli.json).",
              "The audit checks pairing, coverage of workloads/radii/policies, matching input/result hashes, positive timings, zero candidate oracle errors, complete count agreement and the recorded unknown-read correction. It also requires this report to match the CSV exactly. `--check-source` checks the current kernel and benchmark-script snapshots.", "",
              "## Limits", "",
              "These performance inputs are simulated reads, including the workload against the real public guide library. This pass does not retime public FASTQs, competitors, multi-thread scaling, indel matching or downstream inference.",
              "The repeat ranges describe runner variability, not confidence intervals for biological screens. The large-library oracle is sampled; zero observed mismatches is not a guarantee for every possible input.",
              "Existing Bowtie/Edlib comparison reports retain their historical commands and measurements. These before/after ratios must not be multiplied by historical competitor ratios or generalized to other error models.", ""]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--check-source", action="store_true")
    args = parser.parse_args()
    text = render()
    if args.check_source:
        metadata = json.loads((RAW / "hamming_seed_layout.json").read_text())
        recorded_script = next(Path(p) for p in metadata["files"] if p.endswith("/scripts/bench_hamming_seed_layout.py"))
        recorded_root = recorded_script.parents[1]
        for current in (ROOT / "src/qdalign.c", ROOT / "scripts/bench_hamming_seed_layout.py"):
            recorded = recorded_root / current.relative_to(ROOT)
            assert metadata["files"][str(recorded)] == hashlib.sha256(current.read_bytes()).hexdigest(), f"source changed: {current}"
        metadata = json.loads((RAW / "hamming_seed_layout_cli.json").read_text())
        current = ROOT / "scripts/bench_hamming_count_layout.py"
        recorded = next(p for p in metadata["files"] if p.endswith("/scripts/bench_hamming_count_layout.py"))
        assert metadata["files"][recorded] == hashlib.sha256(current.read_bytes()).hexdigest(), f"source changed: {current}"
    if args.check:
        assert REPORT.read_text() == text, "report does not match raw evidence; regenerate it"
        print("Hamming seed-layout evidence: PASS")
    else:
        REPORT.parent.mkdir(parents=True, exist_ok=True)
        REPORT.write_text(text)


if __name__ == "__main__":
    main()
