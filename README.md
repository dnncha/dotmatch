# DotMatch

[![CI](https://github.com/Dnncha/dotmatch/actions/workflows/ci.yml/badge.svg)](https://github.com/Dnncha/dotmatch/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![Citation](https://img.shields.io/badge/cite-CITATION.cff-green.svg)](CITATION.cff)

Fast exact short-DNA edit distance, threshold matching, and barcode-style assignment.

The target is narrow on purpose:

> Exact global edit-distance matching for short DNA barcodes, primers, adapters, UMIs, and guide-style sequences.

DotMatch is starting as a short-DNA matching and barcode assignment engine. The long-term direction is a broader high-performance sequence matching toolkit, with genome-scale alignment treated as a separate future layer rather than a current claim.

This is not a general-purpose aligner today. It is a small C core with a stable low-level ABI, CLI, Python ctypes bindings, correctness tests, and reproducible benchmarks.

DotMatch Core is Apache-2.0 open source. The core assignment engine, CLI, audit, basic QC report, validation harnesses, and public schemas are intended to stay inspectable and citable. See [Open-Core Boundary](docs/open-core-strategy.md) for the planned split between the open scientific trust layer and future proprietary production/workbench infrastructure.

## Why Trust It

- Claims are tied to checked artifacts in the [Scientific Claim Ledger](docs/scientific-claims.md).
- Strict CRISPR gates currently pass from the committed raw benchmark evidence: `make public-crispr-claim-gate` and `make crispr-sota-gate`.
- Barcode and raw BCL state-of-the-art gates intentionally fail until their real-data and comparator requirements are met.
- Assignment semantics are deterministic: `unique`, `ambiguous`, `none`, and `invalid` outcomes are explicit rather than silently tie-broken.
- Public schemas, benchmark CSVs, figures, reports, validation harnesses, and contribution rules are in the repository.

## Current Status

`v0.1.0` includes:

- exact global edit distance;
- Myers 64-bit bit-vector kernel when one sequence is `<=64 bp`;
- thresholded `distance <= k` queries with early rejection;
- batch many-read vs many-target assignment;
- unique/ambiguous/no-match result semantics;
- tiny CLI for pairwise and file-based batch workflows;
- Python `dotmatch` console workflow for FASTQ/FASTQ.gz count tables;
- target-set auditing for one-edit ambiguity risk;
- indexed-vs-exhaustive native validation mode;
- native FASTQ/FASTQ.gz `count` with `--metric hamming|levenshtein`;
- optional one-base Levenshtein indel windows and guide-offset detection for CRISPR-style reads;
- Python `dotmatch` module via ctypes;
- Python source/wheel builds that bundle the native core on Linux and macOS;
- deterministic fuzz tests against the DP oracle;
- microbenchmarks and synthetic barcode batch benchmarks;
- optional Python Edlib comparison.

It does **not** yet include:

- semi-global/infix alignment;
- traceback/CIGAR;
- wildcard `N` semantics;
- native SeqAn/Parasail benchmark harnesses;
- SIMD/NEON-specific implementation;
- published PyPI/Bioconda release artifacts.

## Quickstart

```bash
git clone https://github.com/Dnncha/dotmatch.git
cd dotmatch
make

./dotmatch dist ACGT AGGT
./dotmatch leq 1 ACGT AGGT
make cli-test
```

Python source install:

```bash
python3 -m pip install .
python3 -c "import dotmatch; print(dotmatch.distance('ACGT', 'AGGT'))"
```

Release package verifier:

```bash
python3 -m pip install build
make python-package-test
make publication-ready
```

## Build

```bash
make
make shared
make test
make coverage
make asan
```

`make coverage` runs the native C tests and CLI fixture suite against an instrumented build, writes text/JSON/HTML reports under `build/coverage/`, and currently enforces at least 75% line coverage across `src/qdalign.c` and `src/qda.c`.

## CLI

Pairwise:

```bash
./dotmatch dist ACGT AGGT
# 1

./dotmatch leq 1 ACGT AGGT
# true
```

Batch assignment:

```bash
cat > barcodes.tsv <<'EOF'
bc0	ACGT
bc1	AGGT
bc2	ACGA
EOF

cat > reads.tsv <<'EOF'
r0	ACGT
r1	ACGC
r2	TTTT
EOF

./dotmatch assign 1 barcodes.tsv reads.tsv
```

Input files accept either one sequence per line or `id<TAB>sequence`.

FASTQ barcode assignment:

```bash
./dotmatch fastq-assign \
  --barcodes barcodes.tsv \
  --reads reads.fastq.gz \
  --barcode-start 0 \
  --barcode-length 16 \
  --k 1 \
  --out assignments.tsv
```

`fastq-assign` streams the FASTQ input, supports `.gz` reads via zlib, uses the reusable target index for `k=0` and `k=1`, and writes deterministic `unique`/`ambiguous`/`none`/`invalid` assignment rows.

FASTQ barcode demultiplexing:

```bash
./dotmatch demux \
  --barcodes barcodes.tsv \
  --reads pooled.fastq.gz \
  --barcode-start 0 \
  --barcode-length 8 \
  --k 1 \
  --metric hamming \
  --out-dir demuxed \
  --summary demux.qc.json \
  --assignments demux.assignments.tsv \
  --ambiguous-out ambiguous.fastq \
  --unmatched-out unmatched.fastq
```

`demux` writes one FASTQ per uniquely assigned barcode under `--out-dir`. Ambiguous reads are never assigned silently; use `--ambiguous-out` and `--unmatched-out` to inspect rejected reads. This demux surface targets fixed-position inline barcodes in single-end FASTQ/FASTQ.gz. Use `--barcode-length auto` when the barcode sheet contains multiple barcode lengths; prefix-overlapping exact matches are reported as ambiguous rather than silently assigned to the longest prefix.

Classic Illumina BCL demultiplexing:

```bash
./dotmatch bcl-demux \
  --run-folder 240101_TEST_RUN \
  --sample-sheet SampleSheet.csv \
  --out-dir bcl_demuxed \
  --barcode-mismatches 1 \
  --summary bcl.summary.json

./dotmatch bcl-validate \
  --dotmatch-out bcl_demuxed \
  --truth-out bclconvert_output
```

`bcl-demux` currently supports the first classic per-cycle BCL milestone: `RunInfo.xml`, sample-sheet v1/v2 data sections, `Data/Intensities/BaseCalls/L00*/C*.1/s_*_*.bcl(.gz)`, and `.filter` files. It writes sample FASTQ.gz, `Undetermined` FASTQ.gz, `Demultiplex_Stats.csv`, `SampleSheet.normalized.csv`, and summary JSON. CBCL/NovaSeq-style input is a planned next milestone, not a current claim.

Native count-table workflow:

```bash
./dotmatch count \
  --targets guides.csv \
  --reads sample_R1.fastq.gz \
  --target-start 0 \
  --target-length 20 \
  --k 1 \
  --metric levenshtein \
  --ambiguity-policy best \
  --indel-window 1 \
  --out counts.tsv \
  --target-counts-long target_counts.long.tsv \
  --sample-qc sample_qc.tsv \
  --assignments assignments.tsv \
  --summary summary.json \
  --report report.html \
  --report-audit-dir audit/ \
  --report-unmatched top_unmatched.tsv
```

Targets may be tab-separated `target_id<TAB>target_seq<TAB>gene` or MAGeCK-style CSV with headers such as `id,gRNA.sequence,Gene`. The count table reports exact reads, one-substitution corrections, one-insertion corrections, one-deletion corrections, other corrections, total assigned reads, and whether the target has a nearby target that can create `k`-edit ambiguity. `--target-counts-long` writes one row per sample/target with the same provenance fields, and `--sample-qc` writes assignment rate, exact/rescued/ambiguous/no-match rates, target coverage, zero-count targets, Gini index, top-1% dominance, and candidate-verification totals.

`--report report.html` writes a native HTML run summary with assignment rates, exact/rescued/ambiguous/no-match breakdowns, library coverage, candidate-verification totals, and warnings for high ambiguous or no-match rates. If `--report-audit-dir` or `--report-unmatched` are supplied, the report also embeds a library-audit summary and top-unmatched preview. It is intentionally deterministic and self-contained so it can be archived with the count matrix.

Use `--metric hamming` for a guide-counter-style one-mismatch/no-indel comparison. Use `--metric levenshtein --indel-window 1` when the workflow should recover one-base insertions/deletions around the extracted target window. `--auto-offset N` samples each FASTQ and chooses the best target start within `N` bases of `--target-start` using exact matches.

Assignment ambiguity is explicit. `--ambiguity-policy best` assigns a read when exactly one target has the best distance within `k`; `--ambiguity-policy radius` assigns only when exactly one target is within the whole radius. Ambiguous reads are discarded from counts by default. Use `--ambiguous report` to include ambiguous rows in `assignments.tsv` for diagnostics; they are still not silently counted for a target.

MAGeCK-compatible count table:

```bash
./dotmatch count \
  --targets yusa_library.csv \
  --reads ERR376998.fastq.gz \
  --reads ERR376999.fastq.gz \
  --sample-label plasmid,ESC1 \
  --target-start 23 \
  --target-length 19 \
  --k 1 \
  --metric levenshtein \
  --indel-window 1 \
  --format mageck \
  --out counts.mageck.tsv
```

CRISPR-focused multi-sample count matrix:

```bash
cat > samples.tsv <<'EOF'
sample_id	fastq
control_1	control_1.fastq.gz
control_2	control_2.fastq.gz
drug_1	drug_1.fastq.gz
drug_2	drug_2.fastq.gz
EOF

./dotmatch crispr-count \
  --library guides.csv \
  --samples samples.tsv \
  --guide-start 0 \
  --guide-length 20 \
  --k 1 \
  --metric levenshtein \
  --indel-window 1 \
  --out counts.mageck.tsv \
  --summary qc.json \
  --ambiguous discard
```

`crispr-count` is the CRISPR-facing wrapper around the native count engine. It defaults to MAGeCK-compatible count-matrix output and accepts the same deterministic ambiguity policy as `count`. The QC JSON reports exact assigned reads, one-edit rescued reads, ambiguous reads, unmatched reads, and candidate-verification totals per sample.

Target library audit:

```bash
./dotmatch audit \
  --targets guides.tsv \
  --k 1 \
  --audit-mode auto \
  --out-dir audit/
```

Native audit writes `audit_summary.tsv`, `audit_summary.json`, `collision_pairs.tsv`, `collision_clusters.tsv`, and `target_safety.tsv`. It reports duplicates, minimum edit distance, pairs at distances 0/1/2, whether the library is safe for `k=0`, `k=1`, and `k=2` correction, and per-target nearest-neighbor risk under the same edit-distance semantics used by assignment. `--audit-mode exact` uses exhaustive pairwise distances. `--audit-mode fast` uses one-edit variant indexing for large-library `k=1` safety and reports `not_computed`/`null` for `k=2`-only fields.

For `k=1`, audit also writes `ambiguous_variants.tsv` and records `ambiguous_query_variants_k1` in the summary. This enumerates exact one-edit query variants that would fall within distance 1 of multiple targets, which is the practical safety question behind one-edit rescue.

Top unmatched diagnosis:

```bash
./dotmatch inspect-unmatched \
  --targets guides.tsv \
  --reads sample_R1.fastq.gz \
  --target-start 0 \
  --target-length 20 \
  --k 1 \
  --offset-window 2 \
  --adapter ACGTTT \
  --low-quality-threshold 20 \
  --top 100 \
  --out top_unmatched.tsv
```

This reports the most frequent unassigned extracted sequences, nearest known target, nearest edit distance, edit class, reverse-complement nearest-target hint, optional offset-shift hint, adapter/primer hint when `--adapter` is supplied, low-quality hint when `--low-quality-threshold` is supplied, and a coarse reason such as `near_known_target_above_k`, `reverse_complement_candidate`, `offset_shift_candidate`, `adapter_or_primer_candidate`, `low_quality_candidate`, `contains_N`, or `wrong_length`.

Validation against the native exhaustive scan path:

```bash
DOTMATCH_LIB=$PWD/libdotmatch.dylib PYTHONPATH=$PWD/python python3 -m dotmatch.cli validate \
  --targets guides.tsv \
  --reads sample_R1.fastq.gz \
  --target-start 0 \
  --target-length 20 \
  --k 1 \
  --sample 100000
```

This validates DotMatch's indexed assignment against DotMatch's exact native scan oracle.

Optional native Edlib validator:

```bash
make edlib-tools
./dotmatch validate \
  --targets guides.tsv \
  --reads sample_R1.fastq.gz \
  --target-start 0 \
  --target-length 20 \
  --k 1 \
  --indel-window 1 \
  --oracle edlib \
  --sample 100000
```

## Python

```bash
make shared
DOTMATCH_LIB=$PWD/libdotmatch.dylib PYTHONPATH=$PWD/python python3
```

```python
import dotmatch

dotmatch.distance("ACGT", "AGGT")
# 1

dotmatch.distance_leq("ACGT", "AGGT", 1)
# True

dotmatch.assign(["ACGT", "ACGC"], ["ACGT", "AGGT", "ACGA"], k=1)
```

Reusable indexed assignment:

```python
matcher = dotmatch.Matcher(["ACGT", "AGGT", "ACGA"])
results, stats = matcher.assign_with_stats(["ACGT", "ACGC"], k=1)
stats.candidates_verified
```

The old `quickdna` Python package and `qda` CLI target are kept as transition aliases while the project moves to the DotMatch name.

## Benchmarks

Benchmark data policy: commit reproducible benchmark evidence, not large
working datasets. Keep curated CSVs under `benchmarks/raw/`, figures under
`benchmarks/figures/`, and reports under `docs/benchmarks/`. Generated run
folders, demultiplexed FASTQs, downloaded public datasets, and comparator
scratch output belong in ignored local paths such as `benchmarks/work/` and
`examples/*/data/`; recreate them with the fetch and benchmark scripts below.

Pairwise and threshold microbenchmark:

```bash
make bench
```

Synthetic barcode batch benchmark:

```bash
make bench-batch
# or a smaller smoke run:
./build/bench_batch 1000
```

Barcode demultiplexing benchmark/report:

```bash
make bench-barcode-demux

# For a real public or in-house inline-barcode dataset:
python3 scripts/bench_barcode_demux.py \
  --reads SRR391079.fastq.gz \
  --barcodes barcodes.tsv \
  --barcode-start 1 \
  --barcode-length auto \
  --k 0 \
  --run-cutadapt \
  --run-hash-splitter
python3 scripts/generate_barcode_demux_report.py
```

The default `make bench-barcode-demux` fixture is for smoke testing the graph/report pipeline and includes a simple exact prefix hash-splitter baseline. Barcode state-of-the-art claims require real public barcode FASTQ inputs and fair competitor rows, at minimum Cutadapt plus a second comparator such as Ultraplex, Je, or an exact hash splitter for the exact-prefix lane.

Claim-grade inline barcode evidence is gated separately:

```bash
python3 scripts/fetch_srp009896_barcode_demo.py \
  --metadata-only \
  --use-public-example-barcodes \
  --require-barcodes

# Or provide an already curated fixed-length barcode sheet:
export DOTMATCH_BARCODE_SOTA_BARCODES=/path/to/real_fixed_length_barcodes.tsv
# or: export DOTMATCH_BARCODE_SOTA_BARCODES_URL=https://...
# or for the public SRP009896 example barcode sheet:
export DOTMATCH_BARCODE_SOTA_USE_PUBLIC_EXAMPLE=1
export DOTMATCH_BARCODE_START=1
export DOTMATCH_BARCODE_K=0
make bench-barcode-sota
make barcode-sota-gate
```

`make barcode-sota-gate` intentionally fails if the SRP009896 workflow has no real barcode sheet, if rows are fixture-only, or if comparator rows are missing. This is distinct from raw BCL/CBCL demultiplexing.
The public SRP009896 example barcode sheet is variable-length (`4-8 bp`) and contains separate run blocks with reused barcode sequences. The fetcher filters the installed sheet to the requested accession when that run column is present. SRP009896 reads include a leading `N`, so use `DOTMATCH_BARCODE_START=1` with the public example sheet, and use the `k=0` exact-prefix lane with `--barcode-length auto` unless you provide a separate fixed-length sheet. A barcode SOTA claim remains blocked until real repeated rows and comparator evidence pass the gate.

Raw BCL demultiplexing benchmark/report:

```bash
make bench-bcl-small
make fetch-10x-bcl-demo
make bcl-competitor-env
make bcl-linux-env
DOTMATCH_BCL_THREADS=8 make bench-bcl-10x

DOTMATCH_BCL_RUN_FOLDER=/path/to/run \
DOTMATCH_BCL_SAMPLE_SHEET=/path/to/SampleSheet.csv \
make bench-bcl-real
DOTMATCH_BCL_REPEATS=5 make bench-bcl-real-repeated

make bcl-sota-gate
```

`bench-bcl-small` uses a generated classic-BCL run folder and is only a smoke benchmark. Raw-BCL state-of-the-art claims require real classic-BCL and CBCL run folders, BCL Convert/bcl2fastq/CUDA-Demux comparator rows where available, repeated timing, and `bcl-validate` zero-mismatch evidence. `make bcl-sota-gate` intentionally fails until those requirements are met.

For the public 10x tiny-BCL demo row:

```bash
make fetch-10x-bcl-demo
python3 scripts/bench_bcl_demux.py \
  --run-folder examples/bcl_demux/data/cellranger-tiny-bcl-1.2.0 \
  --sample-sheet examples/bcl_demux/data/cellranger-tiny-bcl-samplesheet.normalized.csv \
  --workflow-name public_10x_tiny_bcl \
  --detect-competitors \
  --run-installed-competitors
python3 scripts/generate_bcl_demux_report.py
```

Native Edlib assignment comparison:

```bash
make benchmark-report-native
make bench-small
make bench-paper
make figures
```

`make figures` records repeated-run native benchmark CSVs under `benchmarks/raw/` and SVG/PDF figures under `benchmarks/figures/`. Use `DOTMATCH_NATIVE_REPEATS=N` and `DOTMATCH_NATIVE_REPORT_READS=N` to control runtime. `make bench-paper` enables the larger paper matrix, including `12/16/20/24/32 bp`, `96/737/4096/16384/65536` targets, substitution/indel/no-match/ambiguous modes, peak RSS, and mismatch counts.

External competitor scaffold:

```bash
make competitor-env
python3 scripts/bench_competitors.py \
  --barcodes barcodes.tsv \
  --reads reads.fastq.gz \
  --barcode-start 0 \
  --barcode-length 16 \
  --k 1 \
  --dotmatch ./dotmatch \
  --run-cutadapt \
  --run-bowtie2 \
  --run-guide-counter \
  --out docs/benchmarks/external_competitors.csv
```

Cutadapt, Bowtie2, and guide-counter comparisons are opt-in and run only when those tools are installed and pinned. These are workflow comparators, not the exact oracle; the native Edlib scan remains the exact assignment comparator. guide-counter is a particularly important CRISPR comparator because it already provides gzipped FASTQ guide counting with one mismatch and no indels; DotMatch now exposes `--metric hamming` for that fair semantic lane, while its differentiator is exact `k=1` Levenshtein assignment with indels, ambiguity semantics, target audit, and a general known-target engine.

Latest native report:

- [Scientific claim ledger](docs/scientific-claims.md)
- [Methods and citation template](docs/methods-and-citation.md)
- [GitHub launch checklist](docs/github-launch-checklist.md)
- [GitHub launch kit](docs/launch-kit.md)
- [Changelog](CHANGELOG.md)
- [Release process](docs/release-process.md)
- [Native Edlib benchmark report](docs/benchmarks/native/README.md)
- [Raw native Edlib assignment CSV](docs/benchmarks/native/native_edlib_assignment.csv)
- [Real CRISPR benchmark report](docs/benchmarks/real/README.md)
- [Public CRISPR workflow comparator](docs/benchmarks/public_crispr/README.md)
- [Raw real CRISPR Edlib CSV](benchmarks/raw/real_crispr_edlib.csv)
- [Raw public CRISPR workflow CSV](benchmarks/raw/public_crispr_workflow.csv)
- [Raw repeated public CRISPR CSV](benchmarks/raw/public_crispr_repeated.csv)
- [Raw public CRISPR count agreement CSV](benchmarks/raw/count_agreement_summary.csv)
- [Raw public CRISPR Edlib validation CSV](benchmarks/raw/public_crispr_edlib_validation.csv)
- [Innovation positioning](docs/innovation-positioning.md)
- [Usability comparison](docs/usability-comparison.md)
- [Contributing guide](CONTRIBUTING.md)
- [Support policy](SUPPORT.md)

![Native speedup vs Edlib](docs/benchmarks/native/native_speedup_vs_edlib.svg)

![Native candidates per read](docs/benchmarks/native/native_candidates_per_read.svg)

![Native assignment throughput](docs/benchmarks/native/native_assignment_throughput.svg)

Real public CRISPR guide-counting benchmark:

```bash
DOTMATCH_REAL_READS=25 DOTMATCH_REAL_FETCH_RECORDS=25 make bench-real-report
```

This uses the public MAGeCK/Yusa guide library and real FASTQ reads from `ERR376998`/`ERR376999`, compares DotMatch indexed `k=1` assignment to native Edlib exhaustive scan, and writes [the real-data report](docs/benchmarks/real/README.md).

Python Edlib smoke comparison, useful for Python users but not for headline claims:

```bash
make shared
python3 -m pip install edlib matplotlib pandas numpy
python3 scripts/bench_vs_edlib.py
```

Appendix report:

- [Python binding benchmark report](docs/benchmarks/README.md)
- [Raw Edlib comparison CSV](docs/benchmarks/edlib_python.csv)
- [Raw batch assignment CSV](docs/benchmarks/batch_assignment.csv)

Benchmark output is CSV-shaped so it can be redirected and plotted. Record hardware, compiler, OS, command, and git commit with any published result. Headline graphs must use native C/C++ comparators, not Python binding overhead.

## Publication Target

The publication target is not a broad aligner claim. DotMatch should be published when it can show that, for real FASTQ barcode workloads, it is faster, simpler, and more predictable than common demultiplexing paths while matching a dynamic-programming oracle for assignment semantics.

See [Publication Target](docs/publishing-target.md) for the benchmark plan, competitor set, pass/fail criteria, and claim discipline. See [Publication Roadmap](docs/publication-roadmap.md) for the paper-grade product checklist, required workflows, figures, release gates, and remaining gaps. Cutadapt, Bowtie/Bowtie2-based workflows, native Edlib scans, and simple scripting baselines should be treated as competitors where their semantics fit the workload.

See [Public Schemas](docs/schemas.md) for the stable TSV/JSON contracts emitted by the open core, and [Scientific Claim Ledger](docs/scientific-claims.md) for the current supported/blocked claim status.

Real CRISPR guide-counting example:

```bash
cd examples/crispr_guides
python3 ../../scripts/fetch_mageck_demo.py --small --out data
./run.sh
```

Use the same script with `--subsample 1000` for a small real public FASTQ subset, or without `--small`/`--subsample` to fetch the full public MAGeCK/Yusa demo FASTQ files.

Paper workflow benchmark:

```bash
make bench-public-crispr-small
make bench-public-crispr-competitors
make validate-public-crispr-edlib
make count-agreement
make public-crispr-report
make public-crispr-smoke-gate
# repeated public-data runs; defaults to 10k and 100k records/sample, 5 repeats
make bench-public-crispr-repeated
make public-crispr-claim-gate
# full public FASTQ download/run:
make bench-public-crispr
# optional external workflow comparators:
python3 scripts/run_public_crispr_benchmark.py --run-mageck --run-cutadapt --run-bowtie2
```

`make bench-public-crispr-small` downloads a small real FASTQ subsample. The tiny `--small` fixture is only for deterministic example smoke tests.
Use `DOTMATCH_PUBLIC_READ_SIZES=10000,100000` and `DOTMATCH_PUBLIC_REPEATS=5` to control the repeated benchmark. README and manuscript speed claims should use only repeated rows with zero Edlib validation mismatches and explicitly matched semantics.
`make public-crispr-smoke-gate` validates the benchmark machinery with relaxed thresholds. `make public-crispr-claim-gate` is strict: repeated real-data rows, count agreement, and Edlib validation must all pass before any publication-style claim is promoted.

Current MAGeCK/Yusa artifact status: `make public-crispr-claim-gate` passes on five repeated 10k and 100k-record/sample runs with installed MAGeCK and guide-counter comparators, DotMatch exact-vs-MAGeCK exact count agreement, DotMatch Hamming-vs-guide-counter count correlation, and 1,000-read/sample native Edlib validation with zero mismatches. This supports a narrow real-data CRISPR workflow claim, not a broad state-of-the-art alignment claim.

Two-dataset SOTA evidence is stricter and is designed for atlas/Linux full-data runs:

```bash
make fetch-sanson-crispr
make bench-crispr-sota
make validate-crispr-sota-edlib
make count-agreement-sota
make crispr-sota-report
make crispr-sota-gate

# On atlas:
make atlas-crispr-sota

# To restart a long interrupted repeated-run job without discarding completed rows:
PATH="$PWD/build/guide-counter/bin:$PWD/build/competitor-env/bin:$PATH" \
  python3 scripts/run_crispr_sota_repeated.py --run-mageck --run-guide-counter --full --resume

# To focus only on missing full-FASTQ rows for the strict gate:
PATH="$PWD/build/guide-counter/bin:$PWD/build/competitor-env/bin:$PATH" \
  python3 scripts/run_crispr_sota_repeated.py --run-mageck --run-guide-counter --full-only --resume

# To target one dataset during a long full-row resume:
PATH="$PWD/build/guide-counter/bin:$PWD/build/competitor-env/bin:$PATH" \
  python3 scripts/run_crispr_sota_repeated.py --datasets mageck_yusa --run-mageck --run-guide-counter --full-only --resume

# To collect Sanson/Brunello full evidence one sample at a time on limited disk:
PATH="$PWD/build/guide-counter/bin:$PWD/build/competitor-env/bin:$PATH" \
  python3 scripts/run_crispr_sota_repeated.py --datasets sanson_brunello \
    --sanson-samples RepB --run-mageck --run-guide-counter \
    --dotmatch-threads 8 --full-only --resume

# To resume long native Edlib validation by sample:
DOTMATCH_SOTA_VALIDATION_RECORDS=10000 DOTMATCH_SOTA_VALIDATION_SAMPLE=10000 \
  python3 scripts/validate_crispr_sota_edlib.py --jobs 2 --edlib-threads 4 --resume
```

`make crispr-sota-gate` requires both MAGeCK/Yusa and Sanson/Brunello real-data rows, competitor rows, DotMatch Hamming throughput above guide-counter by the configured threshold, count agreement at publication depth, Levenshtein candidate collapse versus exhaustive target count, full-row evidence at expected dataset read depth unless relaxed, and native Edlib validation. README/manuscript SOTA wording should be generated only from rows that pass this gate.

Current two-dataset strict-gate status: `make crispr-sota-gate` passes on the checked raw artifacts. MAGeCK/Yusa full FASTQ rows are present, Sanson/Brunello full exact, Hamming, Levenshtein, and guide-counter rows are present via complete `full_sample` chunks, the full Hamming lane passes the guide-counter speed check on both datasets, and Levenshtein candidate collapse remains below the exhaustive-scan threshold.
The Sanson/Brunello full fetcher streams source gzip members directly into per-sample FASTQ.gz files, avoiding a second `.source.fastq.gz` cache copy, and validates cached/downloaded full FASTQs before they are used as benchmark inputs. The strict gate accepts either one full all-sample row or complete `full_sample` rows for `plasmid`, `RepA`, `RepB`, and `RepC`; this permits sample-by-sample collection on limited disks.

## C API Shape

Core pairwise functions:

- `qdaln_edit_distance`
- `qdaln_edit_distance_leq`
- `qdaln_edit_distance_myers64`
- `qdaln_edit_distance_dp`

Batch assignment:

- `qdaln_match_many`
- `qdaln_assign_many`
- result status: `QDALN_MATCH_NONE`, `QDALN_MATCH_UNIQUE`, `QDALN_MATCH_AMBIGUOUS`, `QDALN_MATCH_INVALID`
- ambiguity policy: `QDALN_POLICY_BEST` or `QDALN_POLICY_RADIUS`
- edit provenance: `QDALN_EDIT_EXACT`, `QDALN_EDIT_K1_SUB`, `QDALN_EDIT_K1_INS`, `QDALN_EDIT_K1_DEL`, `QDALN_EDIT_K2`, `QDALN_EDIT_OTHER`

For v0.1, `N` is treated as a literal byte, not a wildcard.

## Roadmap

### Milestone 1: credible baseline repo

- [x] C core
- [x] CLI
- [x] static/shared library
- [x] correctness oracle
- [x] fuzz tests
- [x] microbenchmarks
- [x] GitHub Actions CI
- [x] Python Edlib comparison

### Milestone 2: short-DNA matching system

- [x] `distance_leq(a, b, k)` with thresholded early rejection
- [x] batch many-read vs many-target assignment
- [x] unique/ambiguous/no-match semantics
- [x] synthetic barcode benchmark
- [x] Python ctypes bindings
- [x] native C/C++ Edlib assignment benchmark
- [x] FASTQ/FASTQ.gz count-table command
- [x] native C count-table command
- [x] MAGeCK-compatible count output
- [x] target-set audit command
- [x] indexed-vs-exhaustive validation command
- [x] optional native Edlib validation helper
- [x] local/GitHub Python wheel builds with bundled native core
- [ ] PyPI manylinux/musllinux Linux wheels

### Milestone 3: stronger performance story

- [ ] length-specialized kernels: `<=16`, `<=32`, `<=64`, `<=128`
- [ ] target profile precomputation API
- [ ] 2-bit DNA encoding option
- [ ] Apple Silicon/ARM NEON experiment
- [ ] adaptive dispatcher by length/error regime
- [x] native Edlib assignment comparison
- [x] repeated native benchmark statistics
- [x] peak RSS reporting in native benchmark CSV
- [x] full synthetic error-mode matrix for paper runs
- [x] public CRISPR FASTQ workflow benchmark driver
- [x] Cutadapt/Bowtie2 external workflow comparator hooks
- [x] strict public CRISPR claim gate passing on repeated 10k and 100k-record/sample MAGeCK/Yusa rows
- [ ] native comparisons against SeqAn and Parasail where applicable

## Claim Discipline

Do not claim broad state-of-the-art alignment performance from the current benchmarks.

The credible claim shape is narrower:

> Fast exact short-DNA global edit-distance and threshold assignment for barcode/primer-style workloads, with correctness verified against dynamic programming and native Edlib assignment scans.

The current native Edlib benchmark supports workload-level assignment speedups, not a generic pairwise-alignment speed claim. Cutadapt/Bowtie/Bowtie2 FASTQ demultiplexing comparisons and native SeqAn/Parasail comparisons are still required before broader ecosystem claims.

Do not lead with exact `k=0` speedups against Edlib scan. Exact lookup should be compared against hash-table baselines. The first high-impact DotMatch story is `k=1` known-target short-DNA assignment, where the index avoids almost all full alignments while preserving exact assignment semantics.
