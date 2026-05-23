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

## Evidence Boundary

These are native Edlib scan microbenchmarks for exact short-DNA assignment workloads, plus simple exact-hash and BK-tree baselines. The largest rows are useful for understanding algorithmic scaling against exhaustive scan, but they are not end-to-end workflow speed claims. Exact `k=0` lookup should be judged against hash-table baselines. For `k=1`, the indexed path is reported only when it has zero correctness disagreements against the exhaustive comparator. Levenshtein `k=2` has native hash-neighborhood pruning coverage for packed A/C/G/T windows up to 32 bases in this regenerated report; use only the recorded rows for `k=2` throughput statements, scoped to fixed-window short-DNA assignment.

## Highest Observed Exhaustive-Scan Microbenchmark Speedups

| n_targets | len | k | error_mode | err | reads_per_sec_dotmatch | reads_per_sec_edlib | verified_per_read | peak_rss_kb | speedup_vs_edlib_native |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 4096 | 16 | 0 | one_substitution | 0.030 | 6289306.30 | 387.70 | 0.00 | 10976 | 16222.10 |
| 4096 | 16 | 0 | one_substitution | 0.010 | 6369425.50 | 395.40 | 0.00 | 7728 | 16212.08 |
| 4096 | 16 | 0 | exact | 0.000 | 6211182.20 | 384.30 | 1.00 | 7488 | 16162.33 |
| 4096 | 16 | 0 | one_substitution | 0.000 | 6289310.90 | 420.00 | 0.00 | 7728 | 15118.52 |
| 4096 | 16 | 0 | one_substitution | 0.005 | 6369425.50 | 423.70 | 0.00 | 7728 | 15032.87 |
| 4096 | 24 | 0 | exact | 0.000 | 4347826.00 | 365.20 | 1.00 | 13072 | 11378.76 |
| 4096 | 24 | 0 | one_substitution | 0.000 | 4273505.80 | 386.00 | 0.00 | 13296 | 10778.07 |
| 4096 | 24 | 0 | one_substitution | 0.005 | 4347826.00 | 423.40 | 0.00 | 14352 | 10180.98 |
| 4096 | 24 | 0 | one_substitution | 0.010 | 4385965.10 | 437.20 | 0.00 | 14432 | 10031.94 |
| 4096 | 24 | 0 | one_substitution | 0.030 | 4347826.00 | 437.80 | 0.00 | 14432 | 10018.19 |
| 4096 | 32 | 0 | exact | 0.000 | 3344480.90 | 340.10 | 1.00 | 14560 | 9801.04 |
| 4096 | 32 | 0 | one_substitution | 0.030 | 3344482.20 | 357.80 | 0.00 | 16336 | 9442.08 |

## Median Speedup Summary

| len | k | n_targets | error_mode | speedup_vs_edlib_native |
| --- | --- | --- | --- | --- |
| 16 | 0 | 4096 | exact | 16162.33 |
| 16 | 0 | 4096 | one_substitution | 15365.43 |
| 24 | 0 | 4096 | exact | 11378.76 |
| 24 | 0 | 4096 | one_substitution | 10048.70 |
| 32 | 0 | 4096 | exact | 9801.04 |
| 32 | 0 | 4096 | one_substitution | 8133.79 |
| 16 | 0 | 737 | exact | 3188.70 |
| 16 | 0 | 737 | one_substitution | 2538.68 |
| 24 | 0 | 737 | exact | 1872.99 |
| 24 | 0 | 737 | one_substitution | 1852.02 |
| 32 | 0 | 737 | exact | 1761.31 |
| 32 | 0 | 737 | one_substitution | 1572.20 |

## Repeated-Run Statistics

| tool | error_mode | n_targets | len | k | reads_per_sec_mean | reads_per_sec_p50 | reads_per_sec_p95 | reads_per_sec_cv | peak_rss_kb_max | mismatches_sum |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dotmatch_indexed | exact | 96 | 16 | 0 | 10275655.60 | 10309286.50 | 10405936.87 | 0.02 | 1712 | 0 |
| dotmatch_indexed | exact | 96 | 16 | 1 | 758549.83 | 761035.00 | 761035.00 | 0.01 | 1888 | 0 |
| dotmatch_indexed | exact | 96 | 24 | 0 | 8074783.20 | 8695656.40 | 8764306.87 | 0.14 | 11936 | 0 |
| dotmatch_indexed | exact | 96 | 24 | 1 | 353918.87 | 348432.10 | 371543.56 | 0.05 | 11968 | 0 |
| dotmatch_indexed | exact | 96 | 32 | 0 | 5968578.10 | 6134970.80 | 6419948.06 | 0.10 | 14560 | 0 |
| dotmatch_indexed | exact | 96 | 32 | 1 | 179151.83 | 189000.20 | 222565.97 | 0.29 | 14560 | 0 |
| dotmatch_indexed | exact | 737 | 16 | 0 | 6714006.10 | 6756753.10 | 6840061.06 | 0.02 | 2704 | 0 |
| dotmatch_indexed | exact | 737 | 16 | 1 | 483200.10 | 473484.80 | 502123.61 | 0.04 | 2800 | 0 |
| dotmatch_indexed | exact | 737 | 24 | 0 | 4469156.17 | 4608294.50 | 4808656.64 | 0.10 | 12032 | 0 |
| dotmatch_indexed | exact | 737 | 24 | 1 | 290842.07 | 303398.10 | 305400.96 | 0.08 | 12064 | 0 |
| dotmatch_indexed | exact | 737 | 32 | 0 | 3645371.97 | 3663003.90 | 3663003.90 | 0.01 | 14560 | 0 |
| dotmatch_indexed | exact | 737 | 32 | 1 | 178016.97 | 182681.80 | 184749.91 | 0.06 | 14560 | 0 |
| dotmatch_indexed | exact | 4096 | 16 | 0 | 6224120.77 | 6211182.20 | 6246120.38 | 0.00 | 7488 | 0 |
| dotmatch_indexed | exact | 4096 | 16 | 1 | 194791.67 | 175963.40 | 258173.72 | 0.33 | 11024 | 0 |
| dotmatch_indexed | exact | 4096 | 24 | 0 | 4068141.47 | 4347826.00 | 4347826.00 | 0.12 | 13072 | 0 |
| dotmatch_indexed | exact | 4096 | 24 | 1 | 201508.43 | 199481.30 | 206902.88 | 0.03 | 14432 | 0 |
