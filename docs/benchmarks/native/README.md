# Native Edlib Benchmark Report

- Platform: `macOS-26.2-arm64-arm-64bit`
- Python: `3.9.6`
- Reads per benchmark case: `1000`
- Repetitions per benchmark case: `3`
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
| 4096 | 16 | 0 | one_substitution | 0.005 | 7874020.80 | 650.10 | 0.00 | 13680 | 18096.75 |
| 4096 | 16 | 0 | one_substitution | 0.000 | 10526313.00 | 797.20 | 0.00 | 13680 | 14933.23 |
| 4096 | 16 | 0 | exact | 0.000 | 11764696.90 | 828.70 | 1.00 | 10128 | 14196.57 |
| 4096 | 16 | 0 | one_substitution | 0.010 | 11764696.90 | 852.00 | 0.00 | 13680 | 13972.73 |
| 4096 | 16 | 0 | one_substitution | 0.030 | 11363640.00 | 855.90 | 0.00 | 13680 | 13276.83 |
| 4096 | 24 | 0 | one_substitution | 0.010 | 7633590.70 | 640.00 | 0.00 | 22112 | 10698.80 |
| 4096 | 32 | 0 | one_substitution | 0.000 | 6024097.00 | 613.70 | 0.00 | 23952 | 10445.26 |
| 4096 | 24 | 0 | exact | 0.000 | 7194244.10 | 735.00 | 1.00 | 15280 | 10310.81 |
| 4096 | 24 | 0 | one_substitution | 0.000 | 7407408.50 | 752.90 | 0.00 | 19040 | 10094.78 |
| 4096 | 24 | 0 | one_substitution | 0.005 | 7812497.40 | 699.50 | 0.00 | 22112 | 10065.06 |
| 4096 | 24 | 0 | one_substitution | 0.030 | 7462681.60 | 664.60 | 0.00 | 23152 | 9699.35 |
| 4096 | 32 | 0 | exact | 0.000 | 6369430.20 | 695.70 | 1.00 | 23744 | 9394.78 |

## Median Speedup Summary

| len | k | n_targets | error_mode | speedup_vs_edlib_native |
| --- | --- | --- | --- | --- |
| 16 | 0 | 4096 | one_substitution | 14452.98 |
| 16 | 0 | 4096 | exact | 14196.57 |
| 24 | 0 | 4096 | exact | 10310.81 |
| 24 | 0 | 4096 | one_substitution | 10079.92 |
| 32 | 0 | 4096 | exact | 9394.78 |
| 32 | 0 | 4096 | one_substitution | 8284.86 |
| 16 | 0 | 737 | exact | 3615.65 |
| 16 | 0 | 737 | one_substitution | 2561.36 |
| 24 | 0 | 737 | exact | 1975.52 |
| 24 | 0 | 737 | one_substitution | 1848.80 |
| 32 | 0 | 737 | exact | 1675.21 |
| 32 | 0 | 737 | one_substitution | 1591.65 |

## Repeated-Run Statistics

| tool | error_mode | n_targets | len | k | reads_per_sec_mean | reads_per_sec_p50 | reads_per_sec_p95 | reads_per_sec_cv | peak_rss_kb_max | mismatches_sum |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dotmatch_indexed | exact | 96 | 16 | 0 | 15394023.83 | 12195133.00 | 21219516.25 | 0.38 | 1824 | 0 |
| dotmatch_indexed | exact | 96 | 16 | 1 | 1119531.03 | 1414427.10 | 1463028.45 | 0.50 | 2016 | 0 |
| dotmatch_indexed | exact | 96 | 24 | 0 | 14041685.13 | 14925363.20 | 15128892.35 | 0.12 | 14064 | 0 |
| dotmatch_indexed | exact | 96 | 24 | 1 | 655481.13 | 645161.30 | 679690.43 | 0.04 | 14096 | 0 |
| dotmatch_indexed | exact | 96 | 32 | 0 | 12736209.47 | 12499995.80 | 13412164.79 | 0.05 | 23344 | 0 |
| dotmatch_indexed | exact | 96 | 32 | 1 | 394514.23 | 428632.70 | 451768.28 | 0.21 | 23344 | 0 |
| dotmatch_indexed | exact | 737 | 16 | 0 | 9291860.43 | 9615386.70 | 11188814.67 | 0.24 | 3200 | 0 |
| dotmatch_indexed | exact | 737 | 16 | 1 | 630874.47 | 664451.80 | 986690.59 | 0.65 | 3984 | 0 |
| dotmatch_indexed | exact | 737 | 24 | 0 | 8314953.80 | 8333334.60 | 8525643.81 | 0.03 | 14208 | 0 |
| dotmatch_indexed | exact | 737 | 24 | 1 | 528465.50 | 554939.00 | 573628.67 | 0.12 | 14208 | 0 |
| dotmatch_indexed | exact | 737 | 32 | 0 | 6034443.80 | 6329117.00 | 7299573.59 | 0.26 | 23344 | 0 |
| dotmatch_indexed | exact | 737 | 32 | 1 | 310177.47 | 340136.10 | 388344.33 | 0.33 | 23344 | 0 |
| dotmatch_indexed | exact | 4096 | 16 | 0 | 10767109.50 | 11764696.90 | 11764696.90 | 0.16 | 10128 | 0 |
| dotmatch_indexed | exact | 4096 | 16 | 1 | 454141.77 | 531349.70 | 554527.67 | 0.34 | 13680 | 0 |
| dotmatch_indexed | exact | 4096 | 24 | 0 | 7383201.60 | 7194244.10 | 7750678.46 | 0.05 | 15280 | 0 |
| dotmatch_indexed | exact | 4096 | 24 | 1 | 270206.53 | 331895.10 | 392470.95 | 0.62 | 23152 | 0 |

## Evidence Boundary

These are native Edlib scan microbenchmarks for exact short-DNA assignment workloads, plus simple exact-hash and BK-tree baselines. The largest rows are useful for understanding algorithmic scaling against exhaustive scan, but they are not end-to-end workflow speed claims. Exact `k=0` lookup should be judged against hash-table baselines. For `k=1`, the indexed path is reported only when it has zero correctness disagreements against the exhaustive comparator. Levenshtein `k=2` has native hash-neighborhood pruning coverage for packed A/C/G/T windows up to 32 bases in this regenerated report; use only the recorded rows for `k=2` throughput statements, scoped to fixed-window short-DNA assignment.
