# DotMatch Performance Improvements - 2026-07-05

This benchmark note covers practical performance improvements in the current
codebase: Python callers were not reaching native exact/Hamming kernels, several
higher-level Python workflows ignored their Hamming metric, and the threaded
native exact-count path mishandled offset windows.

## Implemented Improvements

1. Bound `qdaln_index_assign_hamming_stats` in the Python `ctypes` layer.
2. Bound `qdaln_index_lookup_exact_many_stats` in the Python `ctypes` layer.
3. Bound `qdaln_index_lookup_exact_ascii_many_stats` in the Python `ctypes` layer.
4. Bound `qdaln_index_assign_status_stats` in the Python `ctypes` layer.
5. Added top-level `dotmatch.assign_hamming(...)`.
6. Added top-level `dotmatch.assign_exact(...)`.
7. Added `Matcher.assign_hamming(...)`.
8. Added `Matcher.assign_hamming_with_stats(...)`.
9. Added `Matcher.assign_exact(...)`.
10. Added `Matcher.assign_exact_with_stats(...)`.
11. Added `Matcher.assign_status_with_stats(...)` for status-only early-stop use cases.
12. Routed `assign_dataframe(..., metric="hamming")` through the Hamming kernel.
13. Routed `assign_dataframe(..., metric="exact")` through the exact lookup table.
14. Routed `stream_assign(..., metric="hamming")` through the Hamming kernel.
15. Routed `stream_assign(..., metric="exact")` through the exact lookup table.
16. Routed `dotmatch.tl.assign_features(..., metric="hamming"|"exact")` through the native fast paths.
17. Routed `dotmatch.tl.feature_counts(..., metric="hamming"|"exact")` through the native fast paths.
18. Switched AssaySpec fixed-window start scoring to the Hamming kernel.
19. Replaced streamed assignment TSV `DictWriter` rows with direct fixed-column writes.
20. Fixed threaded native Hamming `k=0` offset-window counting so batch mode evaluates extracted windows, records totals, and falls back to the threaded direct worker for large buffers.
21. Replaced the generic exact-assignment fallback in native offset detection with the exact ASCII lookup path.
22. Routed `inspect-unmatched --k 0` primary fixed-window assignment through exact ASCII lookup.
23. Routed `inspect-unmatched --offset-window ... --k 0` shifted-window hints through exact ASCII lookup.
24. Changed `assignment_summary(...)` to accumulate local integer counters instead of mutating a dictionary per row.
25. Changed `write_assignments_tsv(...)` to keep local summary counters while writing rows.

## Benchmarks

Hardware/environment: local macOS build, `make dotmatch`, `PYTHONPATH=python`,
local `libdotmatch.dylib`, synthetic fixed-window DNA workloads.

### Python matcher, 100,000 reads, 4,096 targets, length 20

| Path | Seconds | Reads/sec | Notes |
| --- | ---: | ---: | --- |
| `Matcher.assign_with_stats(..., k=1)` | 0.367143 | 272,374 | General Levenshtein path |
| `Matcher.assign_hamming_with_stats(..., k=1)` | 0.200717 | 498,214 | 1.83x faster |
| `Matcher.assign_exact_with_stats(...)` | 0.194195 | 514,946 | 1.89x faster for exact windows |
| `Matcher.assign_status_with_stats(..., k=1)` | 0.207492 | 481,946 | 1.77x faster when lower-bound ambiguity details are acceptable |

### Python streaming FASTQ, 200,000 reads, 1,024 targets, length 20

| Path | Seconds | Reads/sec | Notes |
| --- | ---: | ---: | --- |
| `stream_assign(..., metric="levenshtein", k=1)` | 1.666834 | 119,988 | Previous default route |
| `stream_assign(..., metric="hamming", k=1)` | 1.323489 | 151,116 | 1.26x faster with identical assignments for this substitution-only workload |

### Native threaded exact count offset regression

Workload: 200,000 reads, target window starts at offset 5, target length 20,
Hamming `k=0`, 1,024 targets. Before this pass, `--threads 4` fed whole reads
to exact lookup and reported `total_reads=0`, `assigned_unique=0`,
`unmatched=200000`.

| Path | Count phase seconds | Total reads | Assigned exact | Unmatched |
| --- | ---: | ---: | ---: | ---: |
| `--threads 1` | 0.214056 | 200000 | 190000 | 10000 |
| `--threads 4` | 0.206341 | 200000 | 190000 | 10000 |

The threaded path now matches single-thread counts and is slightly faster on the
counting phase for this workload.

### Additional exact-window and streaming-summary checks

These checks compare the patched tree against a temporary baseline binary built
from the checked-in pre-change `src/qda.c`. Outputs were compared byte-for-byte
for the native commands.

| Workload | Baseline seconds | Patched seconds | Impact |
| --- | ---: | ---: | --- |
| `count --metric levenshtein --k 0 --auto-offset 8`, 300,000 reads, 4,096 targets | 0.763317 | 0.684115 | 1.12x faster |
| `inspect-unmatched --k 0`, 300,000 reads, 4,096 targets | 0.118202 | 0.116536 | 1.01x faster |
| `assignment_summary(...)`, 1,000,000 streamed rows | 0.401176 | 0.152804 | 2.63x faster |
