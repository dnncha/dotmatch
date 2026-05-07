# DotMatch Publication Roadmap

This is the execution roadmap for the strongest `k=1` DotMatch paper. It is a product and benchmark checklist, not a current performance claim.

## Central Claim

DotMatch should make exact `k=1` Levenshtein assignment cheap enough that users no longer need to choose between exact matching, approximate heuristics, or full aligner workflows for known short-DNA target sets.

The strongest defensible paper claim is:

> DotMatch turns exact one-edit assignment from an exhaustive-alignment operation into an indexed lookup-like operation while preserving oracle-equivalent `unique`, `ambiguous`, and `no-match` semantics.

## Intended Title

Preferred:

> DotMatch: Streaming Exact One-Edit Barcode and Guide Assignment Without Exhaustive Scanning

Alternative:

> DotMatch: Exact One-Edit Assignment of Known Short-DNA Targets at FASTQ Scale

## Scope

DotMatch is a known-target short-DNA assignment engine. It is not a genome aligner and should not be positioned as a Bowtie2, BWA, minimap2, STAR, or general Edlib replacement.

Primary workflows:

- CRISPR guide counting;
- barcode demultiplexing;
- adapter/primer/amplicon-panel target assignment;
- whitelist-like barcode correction;
- synthetic library read assignment.

Current scientific emphasis:

- exact global Levenshtein distance;
- `k=1` assignment;
- substitutions, insertions, and deletions;
- deterministic ambiguity reporting;
- FASTQ/FASTQ.gz input;
- count-table output.

## Product Requirements

### Assignment Contract

The C API now exposes the stable assignment contract that every CLI, Python, report, and benchmark should use:

- [x] assignment status: unique, ambiguous, no-match, invalid;
- [x] edit provenance: exact, `k1_sub`, `k1_ins`, `k1_del`, `k2`, other, invalid;
- [x] best target index;
- [x] best and second-best distance;
- [x] number of best-distance targets;
- [x] number of targets within the requested radius;
- [x] ambiguity policy: `best` or `radius`.

Required invariant:

```text
indexed_result(query, targets, k, policy)
==
exhaustive_result(query, targets, k, policy)
```

Current implementation note: `qdaln_assign_many(...)` is the native exhaustive assignment oracle for this contract. The indexed paths must continue to match it exactly for `k <= 1` before claims are made.

### Required Commands

- [x] `dotmatch fastq-assign`
- [x] `dotmatch count`
- [x] `dotmatch audit-targets`
- [x] `dotmatch validate` against native exhaustive scan
- [x] optional native Edlib validation helper
- [x] `dotmatch count --format mageck`
- [x] native C `dotmatch count` for maximum FASTQ throughput
- [x] separate `ambiguous.tsv` and `unmatched.tsv` outputs
- [x] `dotmatch count --metric hamming|levenshtein`
- [x] guide-counter-fair one-mismatch/no-indel counting mode via `--metric hamming`
- [x] one-base Levenshtein indel counting windows via `--indel-window 1`
- [x] automatic guide-offset detection via `--auto-offset`
- [x] `dotmatch crispr-count` CRISPR-facing wrapper with `--library`, `--samples`, `--guide-start`, `--guide-length`, MAGeCK-compatible output, and QC JSON
- [x] `dotmatch demux` fixed-position inline barcode FASTQ/FASTQ.gz splitter
- [x] barcode demux smoke benchmark/report target
- [x] `dotmatch bcl-demux` first classic per-cycle BCL run-folder milestone
- [x] `dotmatch bcl-validate` FASTQ.gz output comparison helper
- [x] `dotmatch inspect-unmatched` top failed-assignment diagnosis
- [x] native `dotmatch count --report report.html` summary report
- [x] public schema documentation for count, QC, audit, unmatched, and summary outputs
- [x] open-core boundary documentation

### Ambiguity Semantics

DotMatch supports two assignment policies:

- `--ambiguity-policy best`: assign if exactly one target has the best distance within `k`.
- `--ambiguity-policy radius`: assign only if exactly one target is within the allowed radius.

Ambiguous reads must never be silently counted. Output handling remains:

- `--ambiguous discard`: exclude ambiguous reads from counts.
- `--ambiguous report`: exclude ambiguous reads from counts and emit diagnostics when requested.

Do not add random assignment for the first paper.

### Count Output

The count table must include:

```text
target_id
target_seq
gene
count_exact
count_corrected_substitution
count_corrected_insertion
count_corrected_deletion
count_corrected_other
count_total
ambiguous_nearby
```

Current native count outputs:

- [x] MAGeCK-compatible count matrix via `--format mageck`;
- [x] DotMatch count table with exact/corrected/other/total columns;
- [x] target provenance table via `--target-counts-long`;
- [x] sample QC table via `--sample-qc`;
- [x] summary JSON via `--summary`;
- [x] self-contained HTML report via `--report`;
- [x] optional ambiguous and unmatched read diagnostics.

The target provenance TSV records exact, substitution, insertion, deletion, other, total, and nearby-ambiguity counts per target. The sample QC TSV records assignment rates, exact/rescue breakdown, ambiguous/no-match/invalid counts, observed targets, zero-count targets, Gini index, top-1% read fraction, and verified candidate count.

### Target Audit

`dotmatch audit-targets --targets guides.tsv --k 1` should answer:

> Can this target library safely use one-edit correction?

Required outputs:

- number of targets;
- duplicates;
- pairs within `k`;
- unsafe target IDs;
- minimum pairwise distance where tractable;
- expected ambiguity risk summary.

Current native audit outputs:

- [x] `audit_summary.tsv`;
- [x] `audit_summary.json`;
- [x] `collision_pairs.tsv`;
- [x] `collision_clusters.tsv`;
- [x] `target_safety.tsv`;
- [x] duplicate detection;
- [x] minimum pairwise edit distance where tractable;
- [x] distance-0/1/2 pair counts;
- [x] k1/k2 safety flags;
- [x] nearest target and nearest distance per target;
- [x] exact ambiguous-query-variant enumeration for `k=1`;
- [x] variant-indexed `--audit-mode fast` for large-library `k=1` safety;
- [ ] exact/optimized `k=2` audit for very large 70k-120k guide libraries.

Current limitation: `--audit-mode fast` scales the `k=1` safety question using one-edit variant indexing, but reports `not_computed` for `k=2`-only risk fields. Use `--audit-mode exact` as the small-library oracle.

### Unmatched Diagnosis

`dotmatch inspect-unmatched` should turn failed assignment into actionable debugging.

Current native output:

- [x] top unassigned extracted sequences;
- [x] count per sequence;
- [x] nearest known target;
- [x] nearest edit distance;
- [x] nearest edit class;
- [x] reverse-complement nearest-target hint;
- [x] offset-shift hint via `--offset-window`;
- [x] adapter/primer hint via `--adapter`;
- [x] low-quality hint via `--low-quality-threshold`;
- [x] coarse reason labels: `near_known_target_above_k`, `reverse_complement_candidate`, `offset_shift_candidate`, `adapter_or_primer_candidate`, `low_quality_candidate`, `contains_N`, `wrong_length`, `unknown`.
- [x] optional integration into the HTML report through `--report-unmatched`.

Still needed:

- [ ] automatic adapter/primer discovery from unmatched reads.

### Validation

Users should be able to run:

```bash
dotmatch validate \
  --targets guides.tsv \
  --reads sample.fastq.gz \
  --target-start 0 \
  --target-length 20 \
  --k 1 \
  --sample 100000
```

Current validation compares indexed assignment against DotMatch's native exhaustive scan oracle. `--oracle edlib` is available when `make edlib-tools` has built the pinned native Edlib helper; `--indel-window 1` validates the same variable-length one-edit windows used by Levenshtein counting.

Native test coverage now includes deterministic 1k, 16k, and 65k target panels comparing indexed `k=0`/`k=1` assignment against exhaustive `qdaln_match_many(...)`, with exact, substitution, insertion, deletion, ambiguous, N-containing/non-ACGT, invalid, unequal-length, and 8/32 bp edge-length reads. This is correctness coverage, not a performance claim.

Current implementation note: the `k=1` index now uses exact neighbor-key lookup for A/C/G/T targets up to 32 bp, including substitutions and one-base insertions/deletions, followed by exact verification. This is the first candidate-collapse implementation; it should be benchmarked honestly against BK-tree, trie/automaton, and guide-counter-style mismatch methods before any algorithmic novelty claim.

## Benchmark Stories

### 1. Core Algorithm Benchmark

Panel sizes:

```text
96, 737, 4096, 16384, 65536 targets
```

Lengths:

```text
12, 16, 20, 24, 32 bp
```

Error modes:

```text
exact
one substitution
one insertion
one deletion
random no-match
ambiguous by design
```

Methods:

```text
DotMatch indexed
DotMatch exhaustive scan
native Edlib exhaustive scan
hash lookup for k=0
BK-tree k=1
trie / neighbor-generation baseline
Cutadapt where semantics fit
Bowtie2 / BWA / minimap2 as workflow comparators only
```

Metrics:

```text
reads/sec
peak RSS
index build time
candidate targets considered/read
candidate targets verified/read
gzip throughput
correctness agreement
ambiguous rate
no-match rate
```

The key figure must show that Edlib scan verification grows with target count while DotMatch verified candidates/read stays nearly flat.

### 2. CRISPR Guide-Counting Benchmark

Use a real or realistic sgRNA library:

```text
20 nt guides
70k-120k guides if possible
FASTQ.gz input
multiple samples
```

Required outputs:

```text
counts.tsv
summary.json
assignments.tsv
ambiguous/unmatched diagnostics
```

Comparators:

```text
MAGeCK count exact FASTQ mode
guide-counter one-mismatch CRISPR workflow
MAGeCK via aligner/BAM mismatch workflow
Bowtie2-based guide alignment workflow
native Edlib exhaustive scan on subset/full where feasible
```

Current gauntlet status:

- [x] pinned local competitor installer for MAGeCK, Cutadapt, Bowtie2, and guide-counter
- [x] public MAGeCK/Yusa smoke CSV with DotMatch Hamming, DotMatch Levenshtein, MAGeCK exact, guide-counter one-mismatch, Cutadapt, and Bowtie2 workflow rows
- [x] command/version/semantic lane recording for public CRISPR comparator rows
- [x] repeated-run harness with mean, p50, p95, CV, peak RSS, hardware metadata, and raw CSV output
- [x] guide-counter-equivalent count agreement analysis by guide ID
- [x] DotMatch exact-vs-MAGeCK exact count agreement analysis by guide ID
- [x] stratified public-data Edlib validation artifact for exact, corrected, ambiguous, unmatched, and N-containing reads
- [x] claim-grade repeated 10k and 100k-record/sample public FASTQ runs with installed MAGeCK and guide-counter comparators
- [x] Sanson/Brunello dataset fetcher with PRJNA508200/SRP172473 manifest, Brunello library fetch, source checksums, and one-FASTQ-per-sample concatenation
- [x] two-dataset CRISPR SOTA repeated-run harness and strict gate scaffold
- [x] resumable/checkpointed two-dataset repeated-run harness for long atlas jobs
- [x] resumable/checkpointed native Edlib validation harness with sample-level and helper-level parallelism
- [x] claim-grade repeated full-record public FASTQ rows recorded in the checked CRISPR SOTA artifacts
- [x] strict two-dataset CRISPR SOTA gate passing with full rows, count agreement, and Edlib validation artifacts

The adoption result should be:

> DotMatch recovers one-error guide reads directly from FASTQ.gz and emits a count matrix usable by downstream MAGeCK-style analysis.

### 3. Barcode / Whitelist Correction Benchmark

Use:

```text
sample barcode panels
cell-barcode whitelist-like panels
synthetic barcode sets with known edit-distance spacing
```

Current implementation:

- `dotmatch demux` supports fixed-position single-end inline barcodes;
- output is one FASTQ per uniquely assigned barcode plus optional ambiguous/unmatched FASTQs;
- `scripts/bench_barcode_demux.py` can benchmark DotMatch and Cutadapt on a supplied real barcode FASTQ/barcode file;
- `scripts/fetch_srp009896_barcode_demo.py --require-barcodes` records ENA metadata, checksums, barcode sheet metadata, barcode count, barcode length, and extraction position;
- `scripts/check_barcode_sota_gate.py` fails unless rows use a real FASTQ plus real barcode sheet and comparator rows;
- the built-in benchmark fixture is only a smoke test and must not support a state-of-the-art claim.

Report:

```text
valid exact
corrected at one edit
ambiguous at one edit
rejected
```

## Internal Publish Gates

For the main `k=1` claim:

| Requirement | Target |
| --- | ---: |
| Agreement with exhaustive Edlib scan | 100% |
| Median verified candidates/read | <= 1.2 on well-designed panels |
| P95 verified candidates/read | <= 3-5 |
| Speedup vs Edlib scan at 737 targets | >= 10x |
| Speedup vs Edlib scan at 4096 targets | >= 30x |
| Speedup vs Edlib scan at 16k+ targets | >= 50x |
| FASTQ.gz processing | near I/O-bound where possible |
| Ambiguity handling | exact and fully reported |

For `k=0`, compare against hash lookup and do not use Edlib-scan speedups as a headline.

## Figures

Required paper figures:

1. Problem definition: FASTQ read -> extracted guide/barcode -> DotMatch index -> unique/ambiguous/no-match.
2. Candidate collapse: verified candidates/read vs target count.
3. Throughput: reads/sec vs target count.
4. CRISPR guide-counting workflow: runtime, assigned reads, corrected reads, ambiguous reads, count agreement.
5. Ambiguity/rescue tradeoff: k=0 exact, k=1 rescued, k=1 ambiguous, k=1 rejected.

## Reproducibility Requirements

Targets:

```bash
make test
make coverage
make bench-small
make bench-paper
make figures
make public-crispr-smoke-gate
make public-crispr-claim-gate
```

Artifacts:

```text
benchmarks/raw/*.csv
benchmarks/figures/*.svg
benchmarks/figures/*.pdf
benchmarks/README.md
```

Every headline row must record:

- DotMatch version and commit;
- compiler and flags;
- OS and hardware;
- Edlib version;
- command;
- repetitions;
- checksum;
- correctness mismatches.

Claim discipline:

- [x] `public-crispr-smoke-gate` checks local artifact wiring with reduced thresholds.
- [x] `public-crispr-claim-gate` fails unless repeated public CRISPR rows, count agreement, and Edlib validation satisfy publication thresholds.
- [x] strict two-dataset Hamming-vs-guide-counter agreement gate rejects shallow smoke-count artifacts.
- [x] strict Levenshtein candidate-collapse gate checks verified candidates/read relative to target-panel size.
- [x] strict guide-counter speed gate uses a configurable threshold and defaults to "DotMatch must beat guide-counter."
- [x] strict claim gate passing on repeated 10k and 100k-record/sample MAGeCK/Yusa rows with installed MAGECK and guide-counter competitors.
- [x] native C/CLI coverage target with a 75% line-coverage threshold for `src/qdalign.c` and `src/qda.c`.
- [x] strict two-dataset CRISPR SOTA gate passing on the checked 10k/100k/full artifact matrix with installed comparator rows where required.

## Packaging Requirements

Before a serious public launch:

- [ ] source release;
- [ ] Python package;
- [ ] wheels or documented local shared-library install;
- [x] Dockerfile;
- [ ] Bioconda recipe plan;
- [x] `CITATION.cff`;
- [ ] Zenodo DOI;
- [ ] examples with expected outputs.

Target future methods sentence:

> Reads were assigned to the guide library using DotMatch v1.0 with Levenshtein distance <=1. Assignments were retained only when the best target was unique; ambiguous and unmatched reads were discarded.

## Current Honest Status

DotMatch is on track for this paper, but it is not publication-complete.

Already present:

- native exact edit-distance core;
- stable native assignment contract with status, edit provenance, ambiguity policy, best distance, second-best distance, best-target count, and within-radius count;
- indexed `k=0` and `k=1` assignment;
- deterministic unique/ambiguous/no-match semantics;
- best-match and radius-unique ambiguity policies;
- native Edlib assignment benchmark;
- FASTQ/FASTQ.gz assignment CLI;
- Python count-table CLI;
- native C count-table CLI;
- MAGeCK-compatible count output;
- native target audit command with summary, collision pairs, collision clusters, and per-target safety outputs;
- native target provenance TSV and sample QC TSV outputs;
- native HTML count report;
- native unmatched-sequence inspection command;
- public schemas and open-core boundary docs;
- indexed-vs-scan validation command;
- optional native Edlib validation helper;
- generated native Edlib graphs.
- real public MAGeCK/Yusa CRISPR benchmark report against native Edlib scan.
- strict public CRISPR claim gate passing on five repeated 10k and 100k-record/sample Yusa runs, with MAGeCK exact, guide-counter one-mismatch, count agreement, and 1,000-read/sample Edlib validation at zero mismatches.
- strict two-dataset CRISPR SOTA gate passing on checked MAGeCK/Yusa and Sanson/Brunello artifacts, including full-row evidence, count agreement, candidate-collapse checks, and Edlib validation artifacts.
- native fixed-position barcode demux command and barcode-demux report scaffold.

Still needed for the strongest paper:

- exact large-library `k=2` audit for 70k-120k guide panels;
- automatic adapter/primer discovery from unmatched reads;
- CBCL/NovaSeq-style raw run-folder support;
- real classic-BCL and CBCL run-folder benchmark rows against BCL Convert/bcl2fastq/CUDA-Demux with zero validation mismatches;
- real public inline barcode dataset run, starting with SRP009896/SRR391079-SRR391082 or another recognized demux dataset;
- strict barcode SOTA gate passing with a real barcode sheet and repeated Cutadapt plus second-comparator rows;
- fair barcode-demux competitor rows for Ultraplex, Je, deML, sabre/fastx-style splitters, and Illumina demux tools where their input model matches;
- paired-end/dual-index/index-read barcode modes if we want claims beyond fixed-position inline demux;
- native Edlib validation on larger stratified subsets;
- trie baseline beyond the current neighbor-generation lookup;
- benchmark rows for larger target panels up to 16k and 65k;
- packaging beyond the Dockerfile and citation metadata.
