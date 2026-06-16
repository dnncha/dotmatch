# Production CRISPR CPU vs Metal

This report benchmarks the shipping `dotmatch crispr-count` CLI path, not the standalone `bench_gpu_crispr_metal` microbench. CPU remains the assignment authority. Metal is opt-in via `--backend gpu-metal-experimental` and should be paired with `--metal-validate` before production use.

## Decision Rule

- Use CPU when the workload is not Metal-eligible, when `auto` is selected, or when validation fails.
- Treat any Metal speedup as advisory until `metal_validation=passed` and guide-by-guide counts match the CPU shadow run.
- Do not use Metal for Sanson/Brunello-style multi-offset counting; the production path rejects or blocks that lane today.

## Yusa Speedup Figure

![Production CRISPR CPU vs Metal](../../../benchmarks/figures/crispr_cpu_metal_speedup.svg)

## CPU vs Metal Speedup

| dataset | records_per_sample | metal_validate | cpu_reads_per_sec | metal_reads_per_sec | metal_vs_cpu_wall | metal_vs_cpu_throughput | count_match_cpu | metal_validation | count_engine |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mageck_yusa | 10000 | 0 | 317460.7 | 183165.3 | 0.58x | 0.58x | 0 | None | hamming_metal_seed_index |
| mageck_yusa | 100000 | 0 | 2285296.4 | 1580563.3 | 0.69x | 0.69x | 0 | None | hamming_metal_seed_index |


## Raw Rows

| dataset | backend | metal_validate | status | records_per_sample | total_reads | n_targets | offset_mode | auto_offset | wall_seconds | reads_per_sec | backend_effective | count_engine | metal_validation | count_match_cpu | exit_code | device | notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mageck_yusa | cpu | 0 | ok | 10000 | 20000 | 87437 | fixed | 0 | 0.063000 | 317460.7 | cpu | hamming_lookup_direct_single_offset | None |  | 0 | macOS-26.2-arm64-arm-64bit |  |
| mageck_yusa | gpu-metal-experimental | 0 | ok | 10000 | 20000 | 87437 | fixed | 0 | 0.109191 | 183165.3 | gpu-metal-experimental | hamming_metal_seed_index | None | 0 | 0 | macOS-26.2-arm64-arm-64bit |  |
| mageck_yusa | gpu-metal-experimental | 1 | validation_failed | 10000 | 0 |  | fixed | 0 | 0.090965 |  |  |  |  |  | 1 | macOS-26.2-arm64-arm-64bit | Metal validation failed against CPU authority checksum |
| mageck_yusa | cpu | 0 | ok | 100000 | 200000 | 87437 | fixed | 0 | 0.087516 | 2285296.4 | cpu | hamming_lookup_direct_single_offset | None |  | 0 | macOS-26.2-arm64-arm-64bit |  |
| mageck_yusa | gpu-metal-experimental | 0 | ok | 100000 | 200000 | 87437 | fixed | 0 | 0.126537 | 1580563.3 | gpu-metal-experimental | hamming_metal_seed_index | None | 0 | 0 | macOS-26.2-arm64-arm-64bit |  |
| mageck_yusa | gpu-metal-experimental | 1 | validation_failed | 100000 | 0 |  | fixed | 0 | 0.151587 |  |  |  |  |  | 1 | macOS-26.2-arm64-arm-64bit | Metal validation failed against CPU authority checksum |
| sanson_brunello | cpu | 0 | ok | 10000 | 40000 | 77441 | multi | 20 | 0.237331 | 168540.9 | cpu | hamming_lookup_direct | None |  | 0 | macOS-26.2-arm64-arm-64bit |  |
| sanson_brunello | gpu-metal-experimental | 0 | ineligible | 10000 |  |  | multi | 20 |  |  |  |  |  |  |  | macOS-26.2-arm64-arm-64bit | Metal requires fixed single-offset Hamming k=1; offset-mode multi blocks the production Metal path |
| sanson_brunello | gpu-metal-experimental | 1 | ineligible | 10000 |  |  | multi | 20 |  |  |  |  |  |  |  | macOS-26.2-arm64-arm-64bit | Metal requires fixed single-offset Hamming k=1; offset-mode multi blocks the production Metal path |
| sanson_brunello | cpu | 0 | ok | 100000 | 400000 | 77441 | multi | 20 | 0.800146 | 499908.9 | cpu | hamming_lookup_direct | None |  | 0 | macOS-26.2-arm64-arm-64bit |  |
| sanson_brunello | gpu-metal-experimental | 0 | ineligible | 100000 |  |  | multi | 20 |  |  |  |  |  |  |  | macOS-26.2-arm64-arm-64bit | Metal requires fixed single-offset Hamming k=1; offset-mode multi blocks the production Metal path |
| sanson_brunello | gpu-metal-experimental | 1 | ineligible | 100000 |  |  | multi | 20 |  |  |  |  |  |  |  | macOS-26.2-arm64-arm-64bit | Metal requires fixed single-offset Hamming k=1; offset-mode multi blocks the production Metal path |


## Ineligible Or Failed Rows

| dataset | backend | metal_validate | status | records_per_sample | offset_mode | auto_offset | notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| mageck_yusa | gpu-metal-experimental | 1 | validation_failed | 10000 | fixed | 0 | Metal validation failed against CPU authority checksum |
| mageck_yusa | gpu-metal-experimental | 1 | validation_failed | 100000 | fixed | 0 | Metal validation failed against CPU authority checksum |
| sanson_brunello | gpu-metal-experimental | 0 | ineligible | 10000 | multi | 20 | Metal requires fixed single-offset Hamming k=1; offset-mode multi blocks the production Metal path |
| sanson_brunello | gpu-metal-experimental | 1 | ineligible | 10000 | multi | 20 | Metal requires fixed single-offset Hamming k=1; offset-mode multi blocks the production Metal path |
| sanson_brunello | gpu-metal-experimental | 0 | ineligible | 100000 | multi | 20 | Metal requires fixed single-offset Hamming k=1; offset-mode multi blocks the production Metal path |
| sanson_brunello | gpu-metal-experimental | 1 | ineligible | 100000 | multi | 20 | Metal requires fixed single-offset Hamming k=1; offset-mode multi blocks the production Metal path |


## Interpretation

- **Yusa Metal without `--metal-validate` can be faster but is not count-identical to CPU** on the recorded rows (`count_match_cpu=0`). Do not use those counts for downstream analysis.
- **`--metal-validate` fails on Yusa today** because the experimental Metal path disagrees with the CPU authority checksum. The CLI exits non-zero and deletes trust in Metal counts until the mismatch is fixed.
- **At 100k reads/sample on this host, CPU is faster than Metal even before validation overhead.** Production `auto` staying on CPU is consistent with these measurements.
- **Sanson/Brunello remains CPU-only** because `--offset-mode multi` is outside the production Metal eligibility contract.

## Workload Notes

### MAGeCK/Yusa (`mageck_yusa`)

- Fixed guide window: `--guide-start 23 --guide-length 19`.
- Hamming `k=1` with `--ambiguity-policy best`.
- No auto-offset, so the workload stays in the Metal-eligible single-offset lane.
- `auto` intentionally stays on CPU in the shipping CLI; use explicit `--backend gpu-metal-experimental` to test Metal.

### Sanson/Brunello (`sanson_brunello`)

- Uses `--offset-mode multi` and auto-offset detection, matching the public guide-counter lane.
- That multi-offset path is CPU-only today; Metal rows are recorded as `ineligible`.

## Commands

```bash
make bench-crispr-cpu-metal
make crispr-cpu-metal-report
```

Optional full FASTQs:

```bash
DOTMATCH_CPU_METAL_FULL=1 make bench-crispr-cpu-metal
```
