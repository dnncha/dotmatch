# DotMatch Publication Target

DotMatch should be published when it has a credible, reproducible proof that it is faster, simpler, and more predictable for real FASTQ barcode workloads than common demultiplexing paths.

This document is a target, not a current claim.

## Product Scope

DotMatch is a short-DNA matching engine for known target panels. The first publishable product surface is barcode and guide-style assignment:

- exact matching;
- Hamming distance <= 1 where substitutions-only semantics are intended;
- edit distance <= 1 and <= 2 where insertion/deletion errors matter;
- deterministic `unique`, `ambiguous`, and `no-match` results;
- FASTQ input with explicit barcode extraction rules;
- assignment output that can be consumed by larger sequencing pipelines.

DotMatch is not currently a genome-scale aligner. A future genome aligner can share the DotMatch brand, kernels, and benchmark discipline, but it will need a separate reference index and mapping layer.

The first high-impact wedge is `k=1` assignment. Exact `k=0` lookup is important for completeness, but speedups over Edlib scan are not scientifically interesting by themselves because a serious exact-only system should use a hash table.

## Competitors

Benchmark against tools people actually use or might reasonably reach for:

- Cutadapt for adapter/barcode-oriented workflows;
- Bowtie and Bowtie2 for general-aligner-based demultiplexing workflows;
- native Edlib assignment scans as an exact edit-distance baseline;
- exact hash lookup for `k=0`;
- BK-tree, trie, or related approximate lookup baselines for `k=1`;
- simple Python or shell exact-match scripts as a scripting baseline;
- additional demultiplexing tools when their semantics match the workload.

Bowtie and Bowtie2 are useful headline baselines, but Cutadapt is likely the more relevant external competitor for barcode demultiplexing. Do not frame DotMatch as a Bowtie replacement unless the benchmark is explicitly limited to barcode assignment.

## Benchmark Matrix

The minimum publishable benchmark suite should cover:

| Dimension | Values |
| --- | --- |
| Reads | 100k, 1M, 10M |
| Barcode count | 96, 384, 737, 4096 |
| Barcode length | 8, 12, 16, 24, 32 bp |
| Threshold | exact, Hamming <= 1, edit <= 1, edit <= 2 |
| Error rate | 0%, 0.1%, 0.5%, 1%, 3% |
| Error type | substitution-only, indel-only, mixed |
| Ambiguity | no collisions, 1% collisions, 5% collisions |
| Input | uncompressed FASTQ, gzipped FASTQ |

Synthetic data controls scaling. Real FASTQ data proves the system is useful outside toy inputs.

## Correctness Proof

Every benchmark must normalize outputs into a common assignment schema:

```text
read_id
observed_barcode
target_id
best_distance
second_best_distance
match_count
status
```

The correctness oracle is a slow dynamic-programming matcher over the extracted barcode region. The oracle defines `unique`, `ambiguous`, and `no-match` status. Tool outputs that do not expose these states directly must be postprocessed before comparison.

Pass criteria:

- exact mode is bit-for-bit equivalent to the oracle;
- edit-distance modes match the oracle's assignment status and target choice;
- ambiguous ties are reported as ambiguous, not silently broken;
- behavior for `N`, invalid reads, short reads, and low-complexity barcodes is documented;
- every benchmark row aborts on correctness disagreement.

## Speed Proof

Primary metric:

```text
reads_per_second = total FASTQ reads processed / wall_time
```

Secondary metrics:

- CPU seconds;
- peak resident memory;
- input MB/s;
- output MB/s;
- index build time;
- assignment time excluding FASTQ parsing;
- end-to-end time including FASTQ parsing and output;
- runtime variance across repeated runs.

Run at least ten repetitions for headline rows and report mean, p50, p95, and coefficient of variation.

## Simplicity Proof

Score each competitor on operational complexity:

| Metric | Meaning |
| --- | --- |
| Commands | Number of commands required for a complete demux run |
| Intermediate files | Number of generated files before final assignment |
| Native ambiguity status | Whether `ambiguous` is reported without custom postprocessing |
| Native edit-distance reporting | Whether the best and second-best distance are available |
| FASTQ-in assignment-out | Whether the tool can run the target workflow directly |
| Config surface | Number of flags/options needed for the benchmark command |

DotMatch should aim for one command, zero intermediate files, and native assignment semantics.

## Predictability Proof

Predictability means deterministic, measurable behavior:

- zero assignment nondeterminism across repeated runs;
- stable results under input reordering;
- explicit tie handling;
- documented no-match behavior;
- bounded memory growth in streaming mode;
- low runtime variance on warm-cache runs.

Headline rows should report runtime coefficient of variation. A practical target is below 5% on a quiet machine.

## Publication Claim Shape

Acceptable claim:

> DotMatch is a fast, deterministic short-DNA barcode assignment engine for known target panels, with exact `unique`, `ambiguous`, and `no-match` semantics validated against a dynamic-programming oracle and benchmarked on real FASTQ workloads.

Potential headline claim, only if supported:

> On real `k=1` barcode and guide-assignment workloads, DotMatch returns oracle-equivalent assignments while verifying only a small fraction of candidate targets, making it faster and simpler than applicable Cutadapt and Bowtie/Bowtie2-based workflows while preserving explicit ambiguity semantics.

Do not claim:

- broad state-of-the-art alignment performance;
- replacement for Bowtie, Bowtie2, BWA, minimap2, or STAR;
- genome-scale mapping performance;
- support for CIGAR, SAM/BAM, paired-end mapping, or reference-index search until those systems exist.

## Release Gate

DotMatch is ready to publish when:

- the FASTQ demultiplexing CLI exists;
- the benchmark harness can regenerate all headline plots from raw data;
- raw CSV results and exact commands are checked in or archived;
- `k=0` is compared against exact hash lookup, not only Edlib scan;
- `k=1` is compared against at least one approximate lookup baseline such as BK-tree or trie search;
- at least one realistic Cutadapt comparison is included;
- Bowtie/Bowtie2 comparisons are included where their semantics fit;
- all published rows pass oracle correctness checks;
- the README claim is narrower than the strongest benchmark result.
