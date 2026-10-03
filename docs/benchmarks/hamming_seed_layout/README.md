# Hamming seed-layout improvement

This change accelerates indexed fixed-window Hamming assignment and fixes literal unknown-read verification.
The baseline is DotMatch `076b1644f130f3810fc3e5b64936658ee44e2096`. Both engines use the same targets, reads, radii and complete result contract.
These are before/after measurements on one Linux runner, not a new competitor comparison or evidence of overall state of the art.

## What changed

Seed buckets now occupy contiguous ranges. Each record carries the packed target code and length, avoiding linked-list traversal and separate target-array loads.
Each target is verified only in its first matching partition, so packed A/C/G/T read queries need no growing duplicate set or heap allocation.
The 64-bit runner uses 32 bytes per seed record instead of 40; this reduces seed-record storage by 20%, not total process memory by 20%.

Completeness follows from the pigeonhole principle: with at most k substitutions across k+1 disjoint partitions, at least one partition must match exactly.
Choosing the first matching partition visits every candidate once. Duplicate target rows remain distinct; exact matches, best ties, second-best distances and radius match counts remain explicit.
Long sequences and literal symbols retain their supported fallback paths. The index remains immutable during queries, and the public C/Python interfaces are unchanged.

## Correctness correction

The old unknown-read route generated Hamming candidates but verified them with Levenshtein distance. This could understate distances and create a false best-match tie.
For read `ANATGGAT` and targets `AATCGGAT`, `AAACGGAT`, the correct Hamming distances are 3 and 2. The old engine reported an ambiguous distance-2 tie.
The new engine reports target 1 uniquely under `best`, with best distance 2, second-best distance 3 and two targets inside radius 3. The default `radius` policy correctly remains ambiguous.
The same regression covers literal `R` and lowercase `a`; none is a wildcard. Historical Hamming k=3 results involving non-ACGT reads should be regenerated before reuse.

Across 612,000 read/radius results, the new engine passed 136,096 independent exhaustive byte-oracle checks with zero mismatches.
All A/C/G/T read results matched the baseline in every field. The fully oracle-checked dense literal corpus corrected 10 results, including 4 status/unique-assignment changes; baseline oracle errors are retained in the CSV.
Native tests also exhaust all A/C/G/T words through length 5 at k=1..3, including duplicate target rows; test every length 0..33, mixed target lengths and dense shared seeds; and assert zero query heap allocations on a dense packed workload.
Python and full-count CLI regressions exercise both ambiguity policies. This validates assignment semantics, not biological source identity or gene-hit accuracy.

## Complete counting commands

The public Yusa library contains 87,437 19-base guide rows. FASTQs are simulated A/C/G/T reads: equal exact, one-, two-, three-substitution and random-decoy groups.
Each pair runs a fresh single-thread subprocess, including library/index construction, gzip FASTQ parsing, counting and MAGeCK TSV/JSON output. Five paired repeats alternate engine order.
All guide counts and read-outcome/candidate totals agree between engines and with the validated native batch results. File hashes are recorded; timing excludes input generation and output validation.

| Radius k | Policy | Baseline seconds | New seconds | Paired median speedup (min–max) |
| --- | --- | --- | --- | --- |
| 2 | best | 0.491 | 0.308 | 1.57× (1.22–1.68) |
| 2 | radius | 0.455 | 0.304 | 1.50× (0.85–1.51) |
| 3 | best | 2.399 | 0.445 | 5.24× (4.35–6.74) |
| 3 | radius | 2.299 | 0.387 | 5.69× (4.82–6.01) |

## Native batch assignment

Times include statistics collection and exclude index construction, Python marshalling and output checks. Each engine is warmed once before the five alternating paired repeats.
Large libraries use 256 uniformly sampled reads per radius for the independent oracle. Both dense 1,024-target workloads exhaust every read. All target rows participate in each oracle scan.

| Workload | k | Reads | Baseline reads/s | New reads/s | Paired median speedup (min–max) |
| --- | --- | --- | --- | --- | --- |
| dense_literal_reads_1024_8 | 1 | 3,000 | 662,658 | 1,503,030 | 2.53× (1.07–3.24) |
| dense_literal_reads_1024_8 | 2 | 3,000 | 118,994 | 489,150 | 4.11× (3.22–4.91) |
| dense_literal_reads_1024_8 | 3 | 3,000 | 41,314 | 104,871 | 2.54× (2.15–3.79) |
| dense_shared_seed_1024_8 | 1 | 30,000 | 724,241 | 2,272,987 | 3.04× (2.98–3.50) |
| dense_shared_seed_1024_8 | 2 | 30,000 | 119,442 | 602,261 | 5.02× (4.21–5.11) |
| dense_shared_seed_1024_8 | 3 | 30,000 | 46,289 | 193,610 | 4.29× (4.07–4.45) |
| random_2048_20 | 1 | 30,000 | 1,793,235 | 2,027,716 | 1.17× (1.02–1.39) |
| random_2048_20 | 2 | 30,000 | 1,684,664 | 1,863,026 | 1.11× (0.92–1.26) |
| random_2048_20 | 3 | 30,000 | 1,343,242 | 1,492,253 | 1.10× (0.75–1.28) |
| random_65536_20 | 1 | 30,000 | 1,086,510 | 1,386,844 | 1.42× (0.74–1.52) |
| random_65536_20 | 2 | 30,000 | 336,933 | 962,210 | 3.10× (2.51–3.11) |
| random_65536_20 | 3 | 30,000 | 35,192 | 434,735 | 12.33× (9.69–13.75) |
| random_65536_32 | 1 | 30,000 | 951,471 | 1,055,857 | 1.18× (0.97–1.22) |
| random_65536_32 | 2 | 30,000 | 745,582 | 1,035,847 | 1.45× (1.12–1.50) |
| random_65536_32 | 3 | 30,000 | 514,992 | 795,031 | 1.54× (1.28–1.68) |
| yusa_library_simulated_reads | 1 | 30,000 | 1,360,611 | 2,120,819 | 1.44× (0.71–3.12) |
| yusa_library_simulated_reads | 2 | 30,000 | 141,712 | 787,859 | 4.98× (3.81–8.51) |
| yusa_library_simulated_reads | 3 | 30,000 | 12,710 | 228,107 | 18.18× (16.79–20.10) |

Exact k=0 control medians ranged from 0.89× to 1.09×; no exact-kernel improvement is claimed.
Index construction is measured once per engine/workload and recorded separately in the native CSV. It is included in every complete-count timing above.

## Reproduce and audit

Recorded build: Linux x86-64, `cc -O3 -std=c11 -Wall -Wextra -Wpedantic -Iinclude -mavx2` (plus `-fPIC` for the shared object), with identical flags for baseline and new engine.
Compiler: `cc (Debian 14.2.0-19) 14.2.0`. Platform: `Linux-6.18.44-x86_64-with-glibc2.41`. NumPy oracle: `2.3.5`.
Download the public library with the existing example fetcher if it is absent. Its small FASTQ download is not used by these benchmarks.

```bash
python3 scripts/fetch_mageck_demo.py --out benchmarks/real/data --subsample 25
mkdir -p benchmarks/work/hamming-seed-baseline
git archive 076b1644f130f3810fc3e5b64936658ee44e2096 Makefile pyproject.toml include src | tar -x -C benchmarks/work/hamming-seed-baseline
make -C benchmarks/work/hamming-seed-baseline shared dotmatch CC=cc
make shared dotmatch CC=cc
python3 scripts/bench_hamming_seed_layout.py \
  --baseline-lib benchmarks/work/hamming-seed-baseline/libdotmatch.so \
  --baseline-source benchmarks/work/hamming-seed-baseline/src/qdalign.c \
  --baseline-ref 076b1644f130f3810fc3e5b64936658ee44e2096 \
  --public-library benchmarks/real/data/yusa_library.csv
python3 scripts/bench_hamming_count_layout.py \
  --baseline-bin benchmarks/work/hamming-seed-baseline/dotmatch \
  --library benchmarks/real/data/yusa_library.csv
python3 scripts/report_hamming_seed_layout.py
python3 scripts/report_hamming_seed_layout.py --check --check-source
make test cli-test python-test python-package-test repository-ready
make coverage
```

Raw native [CSV](https://github.com/dnncha/dotmatch/blob/main/benchmarks/raw/hamming_seed_layout.csv) and [provenance](https://github.com/dnncha/dotmatch/blob/main/benchmarks/raw/hamming_seed_layout.json); complete-count [CSV](https://github.com/dnncha/dotmatch/blob/main/benchmarks/raw/hamming_seed_layout_cli.csv) and [provenance](https://github.com/dnncha/dotmatch/blob/main/benchmarks/raw/hamming_seed_layout_cli.json).
The audit checks pairing, coverage of workloads/radii/policies, matching input/result hashes, positive timings, zero candidate oracle errors, complete count agreement and the recorded unknown-read correction. It also requires this report to match the CSV exactly. `--check-source` checks the current kernel and benchmark-script snapshots.

## Limits

These performance inputs are simulated reads, including the workload against the real public guide library. This pass does not retime public FASTQs, competitors, multi-thread scaling, indel matching or downstream inference.
The repeat ranges describe runner variability, not confidence intervals for biological screens. The large-library oracle is sampled; zero observed mismatches is not a guarantee for every possible input.
Existing Bowtie/Edlib comparison reports retain their historical commands and measurements. These before/after ratios must not be multiplied by historical competitor ratios or generalized to other error models.
