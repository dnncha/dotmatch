# Native Edlib Benchmark Report

- Platform: `macOS-26.2-arm64-arm-64bit`
- Python: `3.9.6`
- Reads per benchmark case: `10`
- Repetitions per benchmark case: `2`
- Comparator: native Edlib C/C++ API, `EDLIB_MODE_NW`, `EDLIB_TASK_DISTANCE`, fixed threshold `k`.
- Additional baselines: exact hash lookup for `k=0`; BK-tree and neighbor lookup approximate baselines for `k=1`.
- Assignment mismatches recorded across all rows: `0`.
- Every benchmark run aborts on assignment disagreement between DotMatch and native Edlib scan.

![Native speedup vs Edlib](native_speedup_vs_edlib.svg)

![Native candidates per read](native_candidates_per_read.svg)

![Native assignment throughput](native_assignment_throughput.svg)

## Highest Observed Microbenchmark Speedups

| n_targets | len | k | error_mode | err | reads_per_sec_dotmatch | reads_per_sec_edlib | verified_per_read | peak_rss_kb | speedup_vs_edlib_native |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 4096 | 16 | 0 | one_substitution | 0.005 | 4166669.30 | 248.10 | 0.00 | 2960 | 16953.44 |
| 4096 | 16 | 0 | one_substitution | 0.000 | 3333340.30 | 256.85 | 0.00 | 2880 | 12991.50 |
| 4096 | 16 | 0 | one_substitution | 0.030 | 3749998.75 | 299.85 | 0.00 | 3072 | 12214.04 |
| 4096 | 16 | 0 | exact | 0.000 | 2666668.35 | 251.00 | 1.00 | 2576 | 10711.79 |
| 4096 | 24 | 0 | one_substitution | 0.030 | 3499997.35 | 321.85 | 0.00 | 4176 | 10513.38 |
| 4096 | 24 | 0 | one_substitution | 0.005 | 2666671.25 | 304.65 | 0.00 | 4176 | 8718.31 |
| 4096 | 24 | 0 | one_substitution | 0.000 | 2666668.35 | 308.60 | 0.00 | 4176 | 8563.76 |
| 4096 | 16 | 0 | one_substitution | 0.010 | 2250000.70 | 264.95 | 0.00 | 3040 | 8489.23 |
| 4096 | 32 | 0 | one_substitution | 0.000 | 2499999.20 | 297.70 | 0.00 | 4912 | 8398.06 |
| 4096 | 24 | 0 | exact | 0.000 | 2666663.15 | 330.45 | 1.00 | 4080 | 8058.60 |
| 4096 | 32 | 0 | one_substitution | 0.005 | 1833331.25 | 230.90 | 0.00 | 4912 | 7948.68 |
| 4096 | 24 | 0 | one_substitution | 0.010 | 2666663.15 | 336.90 | 0.00 | 4176 | 7901.77 |

## Median Speedup Summary

| len | k | n_targets | error_mode | speedup_vs_edlib_native |
| --- | --- | --- | --- | --- |
| 16 | 0 | 4096 | one_substitution | 12747.04 |
| 16 | 0 | 4096 | exact | 10711.79 |
| 24 | 0 | 4096 | one_substitution | 8718.31 |
| 24 | 0 | 4096 | exact | 8058.60 |
| 32 | 0 | 4096 | one_substitution | 7099.98 |
| 32 | 0 | 4096 | exact | 6074.37 |
| 16 | 0 | 737 | exact | 2418.54 |
| 16 | 0 | 737 | one_substitution | 2316.79 |
| 24 | 0 | 737 | one_substitution | 1994.31 |
| 24 | 0 | 737 | exact | 1579.75 |
| 32 | 0 | 737 | one_substitution | 1464.31 |
| 32 | 0 | 737 | exact | 1372.04 |

## Repeated-Run Statistics

| tool | error_mode | n_targets | len | k | reads_per_sec_mean | reads_per_sec_p50 | reads_per_sec_p95 | reads_per_sec_cv | peak_rss_kb_max | mismatches_sum |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dotmatch_indexed | exact | 96 | 16 | 0 | 4166669.30 | 4166669.30 | 4916665.40 | 0.28 | 1472 | 0 |
| dotmatch_indexed | exact | 96 | 16 | 1 | 100399.60 | 100399.60 | 108941.05 | 0.13 | 1616 | 0 |
| dotmatch_indexed | exact | 96 | 24 | 0 | 3333324.10 | 3333324.10 | 3333324.10 | 0.00 | 3824 | 0 |
| dotmatch_indexed | exact | 96 | 24 | 1 | 63122.25 | 63122.25 | 64376.71 | 0.03 | 3888 | 0 |
| dotmatch_indexed | exact | 96 | 32 | 0 | 2666663.15 | 2666663.15 | 3266658.00 | 0.35 | 4400 | 0 |
| dotmatch_indexed | exact | 96 | 32 | 1 | 135257.80 | 135257.80 | 145878.70 | 0.12 | 4448 | 0 |
| dotmatch_indexed | exact | 737 | 16 | 0 | 4166669.30 | 4166669.30 | 4916665.40 | 0.28 | 1840 | 0 |
| dotmatch_indexed | exact | 737 | 16 | 1 | 87746.30 | 87746.30 | 89131.76 | 0.02 | 1936 | 0 |
| dotmatch_indexed | exact | 737 | 24 | 0 | 2666671.25 | 2666671.25 | 3266673.39 | 0.35 | 4016 | 0 |
| dotmatch_indexed | exact | 737 | 24 | 1 | 60086.60 | 60086.60 | 61223.39 | 0.03 | 4064 | 0 |
| dotmatch_indexed | exact | 737 | 32 | 0 | 2916669.75 | 2916669.75 | 3291673.24 | 0.20 | 4528 | 0 |
| dotmatch_indexed | exact | 737 | 32 | 1 | 120490.60 | 120490.60 | 128932.15 | 0.11 | 4528 | 0 |
| dotmatch_indexed | exact | 4096 | 16 | 0 | 2666668.35 | 2666668.35 | 3266673.10 | 0.35 | 2576 | 0 |
| dotmatch_indexed | exact | 4096 | 16 | 1 | 36357.60 | 36357.60 | 55359.93 | 0.82 | 3104 | 0 |
| dotmatch_indexed | exact | 4096 | 24 | 0 | 2666663.15 | 2666663.15 | 3266658.00 | 0.35 | 4080 | 0 |
| dotmatch_indexed | exact | 4096 | 24 | 1 | 31952.15 | 31952.15 | 33600.64 | 0.08 | 4176 | 0 |

## Evidence Boundary

These are native Edlib scan microbenchmarks for exact short-DNA assignment workloads, plus simple exact-hash and BK-tree baselines. The largest rows are useful for understanding algorithmic scaling against exhaustive scan, but they are not end-to-end workflow speed claims. Exact `k=0` lookup should be judged against hash-table baselines. For `k=1`, the indexed path is reported only when it has zero correctness disagreements against the exhaustive comparator. Levenshtein `k=2` uses the exhaustive assignment path. The direct Levenshtein `k=1` threshold helper is separately covered by `make test` through `tests/test_qdalign_threshold_alloc.c`, which asserts exact/substitution/insertion/deletion threshold cases do not allocate heap memory.
