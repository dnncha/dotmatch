# Experimental GPU Acceleration Benchmark

This report is a skunk-works GPU evidence lane. It is intentionally not a production speed claim: the Metal path brute-forces packed Hamming `k=1` distances and is compared against DotMatch's existing CPU indexed Hamming assignment with identical output checks.

The decision rule is simple: GPU rows must have zero mismatches before any speed result is considered, and CPU-indexed throughput remains the production baseline unless the GPU path wins end-to-end on real workloads.

## Figure

![GPU speedup](../../../benchmarks/figures/gpu_metal_speedup.svg)

## CPU vs Metal Rows

| n_reads | n_targets | len | k | gpu_reads_per_sec | cpu_reads_per_sec | gpu_vs_cpu_kernel | gpu_vs_cpu_total | mismatches |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20000 | 737 | 20 | 1 | 4575960.0 | 456764.0 | 10.02x | 9.90x | 0 |
| 50000 | 4096 | 20 | 1 | 3757110.0 | 443491.0 | 8.47x | 8.50x | 0 |


## Raw Rows

| tool | backend | status | workload | n_reads | n_targets | len | k | prep_seconds | seconds | total_seconds | reads_per_sec | total_reads_per_sec | pairs_per_sec | checksum | mismatches | device | notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dotmatch_cpu_index | cpu | ok | synthetic_hamming | 20000 | 737 | 20 | 1 | 0.000362625 | 0.0437863 | 0.044149 | 456764 | 453012 | 3.36635e+08 | 124483641 | 0 | Apple M5 | CPU indexed Hamming k=1 baseline; prep is index build |
| dotmatch_gpu_metal | metal | ok | synthetic_hamming | 20000 | 737 | 20 | 1 | 8.9291e-05 | 0.00437067 | 0.00445996 | 4.57596e+06 | 4.48435e+06 | 3.37248e+09 | 124483641 | 0 | Apple M5 | Metal brute-force packed Hamming k=1; prep is shared-buffer allocation and copy |
| dotmatch_cpu_index | cpu | ok | synthetic_hamming | 50000 | 4096 | 20 | 1 | 0.00176629 | 0.112742 | 0.114508 | 443491 | 436650 | 1.81654e+09 | 1714574030 | 0 | Apple M5 | CPU indexed Hamming k=1 baseline; prep is index build |
| dotmatch_gpu_metal | metal | ok | synthetic_hamming | 50000 | 4096 | 20 | 1 | 0.000163958 | 0.0133081 | 0.013472 | 3.75711e+06 | 3.71139e+06 | 1.53891e+10 | 1714574030 | 0 | Apple M5 | Metal brute-force packed Hamming k=1; prep is shared-buffer allocation and copy |


## Scope

This lane tests whether GPU compute is worth productizing. It currently covers fixed-length packed A/C/G/T Hamming `k=1` assignment only. It does not cover Levenshtein indels, BCL conversion, GPU-resident FASTQ parsing, CUDA deployment, or production scheduling.

A future production GPU path needs a real workload gate that includes sequence extraction, packing, transfer or shared-memory preparation, dispatch, readback, and downstream count/QC generation.
