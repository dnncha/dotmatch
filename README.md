# DotMatch

[![CI](https://github.com/Dnncha/dotmatch/actions/workflows/ci.yml/badge.svg)](https://github.com/Dnncha/dotmatch/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![Citation](https://img.shields.io/badge/cite-CITATION.cff-green.svg)](CITATION.cff)

DotMatch is a fast, deterministic tool for assigning short sequencing reads to a known set of short DNA targets.

It is built for guides, barcodes, primers, adapters, and other short target lists where the expected sequences are known in advance.

In many CRISPR and barcode workflows, the real question is simple: which target does this read belong to, and how should ties or near-misses be handled? DotMatch stays focused on that question.

It reports each read as one of:

- `unique`
- `ambiguous`
- `none`
- `invalid`

It reports ties explicitly. There is no SAM/BAM output, no CIGAR, and no reference-mapping layer.

## Why Scientists Use DotMatch

Scientists usually reach for DotMatch when they want a known-target assignment tool that is:

- explicit about ambiguity
- fast on short fixed targets such as guides and barcodes
- usable from the command line in real FASTQ workflows
- simple to audit when a library has near-collisions
- narrow enough that the behavior is easy to reason about

DotMatch fits cases where the target space is already known and the priority is clean, reproducible assignment.

## When To Use It

DotMatch is a good fit for:

- CRISPR guide counting
- inline barcode demultiplexing
- barcode-pair or dual-guide fixed-window counting
- primer or adapter prefix assignment
- feature-barcode or guide-capture per-read assignment
- target-library collision audits before enabling one-edit rescue

Use something else when you need:

- genome or transcriptome alignment
- SAM/BAM output
- traceback or CIGAR strings
- general adapter trimming
- cell-level or UMI-level quantification
- variant calling, consensus generation, or general-purpose mapping

## Why Use DotMatch Instead Of Something Else?

DotMatch is shaped around a specific job.

Many existing tools are excellent. They answer different questions: whole-reference alignment, trimming, demultiplexing with broader workflow assumptions, or downstream single-cell quantification. When the job is "assign each read to one of these known short targets and tell me when that assignment is not unique," a smaller tool is often easier to trust.

DotMatch keeps the contract small:

- short known targets in, deterministic assignments out
- exact and bounded-edit-distance matching
- explicit `unique` / `ambiguous` / `none` / `invalid` outcomes
- count, demux, validation, and audit commands built around that assignment model

That narrow scope is deliberate.

## Quick Start

Build the native CLI:

```bash
git clone https://github.com/Dnncha/dotmatch.git
cd dotmatch
make
```

Try the basic pairwise commands:

```bash
./dotmatch dist ACGT AGGT
./dotmatch leq 1 ACGT AGGT
```

Install the Python package from source:

```bash
python3 -m pip install .
python3 -c "import dotmatch; print(dotmatch.distance('ACGT', 'AGGT'))"
```

## A Small Example

This example shows the core behavior in a small, readable form.

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

Expected output:

```text
mode	read_id	read_seq	target_index	target_seq	distance	status	match_count	second_best_distance
assign	r0	ACGT	0	ACGT	0	unique	3	1
assign	r1	ACGC	0	ACGT	1	ambiguous	2	-1
assign	r2	TTTT	-1		-1	none	0	-1
```

Read `r1` is one edit away from more than one target, so DotMatch reports it as `ambiguous`. That behavior is useful in guide-counting and barcode workflows where forced assignments can distort counts.

## Common Workflows

### CRISPR Guide Counting

```bash
./dotmatch count \
  --targets guides.csv \
  --reads sample_R1.fastq.gz \
  --target-start 0 \
  --target-length 20 \
  --k 1 \
  --metric levenshtein \
  --indel-window 1 \
  --out counts.tsv \
  --sample-qc sample_qc.tsv \
  --summary summary.json
```

For CRISPR-style use, DotMatch can produce guide counts, per-sample QC summaries, assignment tables, and MAGeCK-compatible count matrices.

### Inline Barcode Demultiplexing

```bash
./dotmatch demux \
  --barcodes barcodes.tsv \
  --reads pooled.fastq.gz \
  --barcode-start 0 \
  --barcode-length 8 \
  --k 1 \
  --metric hamming \
  --out-dir demuxed \
  --summary demux.qc.json
```

The demultiplexer is built for fixed-position inline barcodes in FASTQ or FASTQ.gz, with explicit handling of ambiguous and unmatched reads.

### Paired Target Counting

```bash
./dotmatch pair-count \
  --left-targets guides_a.tsv \
  --right-targets guides_b.tsv \
  --reads paired_barcodes.fastq.gz \
  --left-start 0 \
  --left-length 20 \
  --right-start 24 \
  --right-length 20 \
  --k 1 \
  --metric hamming \
  --out pair_counts.tsv \
  --summary pair_summary.json
```

`pair-count` assigns two fixed windows independently and counts only reads where both sides are uniquely assigned.

### Target-Library Audit

```bash
./dotmatch audit \
  --targets guides.tsv \
  --k 1 \
  --audit-mode auto \
  --out-dir audit/
```

This is useful before turning on one-edit rescue. It tells you whether the target library itself creates ambiguity risk.

### AssaySpec Workflows

```bash
dotmatch assay init --template crispr --out assay.toml
dotmatch assay check assay.toml
dotmatch assay plan assay.toml
dotmatch assay run assay.toml
```

AssaySpec is a TOML layer for fixed-window `count`, `demux`, and `pair-count` workflows. It validates the spec, prints the native commands, runs target audit before assignment, writes workflow-facing reports, and records provenance in `assay_manifest.json`. See [DotMatch AssaySpec v1](docs/assayspec.md).

## Output Semantics

DotMatch keeps a small, stable set of assignment outcomes:

- `unique`: exactly one acceptable target assignment
- `ambiguous`: more than one acceptable target assignment
- `none`: no target within the requested distance rule
- `invalid`: the read could not be evaluated under the requested extraction rules

Explicit ambiguity handling is central to the project. Ambiguous reads stay visible and should not be treated as confident assignments.

## Alphabet Policy

DotMatch uses literal-byte matching for target assignment. `A`, `C`, `G`, `T`, `N`, and IUPAC ambiguity symbols are ordinary symbols; `N` and IUPAC codes are not expanded as wildcards.

The public policy string is:

```text
literal-byte; A/C/G/T/N/IUPAC symbols are ordinary byte symbols; no wildcard expansion
```

## Installation Status

Today, the checked install path is from source.

```bash
python3 -m pip install .
```

Public package channels for PyPI, Bioconda, GHCR/BioContainers, and Zenodo DOI are pending or tracked separately. Treat them as available only after the corresponding release checks pass.

For the current state, see:

- [Packaging notes](docs/packaging.md)
- [Distribution release record](docs/distribution-release.json)

Post-release public channel checks stay separate from local release checks:

```bash
make distribution-channels
make workflow-adoption-status
```

## Accuracy, Validation, And Scope

DotMatch has repository checks for:

- native C tests against a dynamic-programming oracle
- deterministic FASTQ assignment behavior
- indexed-versus-exhaustive validation
- benchmark and workflow evidence for specific supported lanes

DotMatch has intentionally narrow evidence boundaries. It is well supported as a known-target assignment tool for short-DNA workflows. Broader claims should stay out of the README.

Before tagging a release, run the consolidated local pre-tag gate:

```bash
make pretag-ready
```

If you are writing documentation, a methods section, or release notes, read:

- [Evidence notes](docs/scientific-claims.md)
- [Methods and citation](docs/methods-and-citation.md)

## Python API

```python
import dotmatch

dotmatch.distance("ACGT", "AGGT")
# 1

dotmatch.distance_leq("ACGT", "AGGT", 1)
# True

dotmatch.assign(["ACGT", "ACGC"], ["ACGT", "AGGT", "ACGA"], k=1)
```

The Python package uses the same native core through `ctypes`, so the CLI and Python surfaces share the same assignment behavior.

## Where To Go Next

- [CRISPR count first run](docs/tutorials/crispr-count-first-run.md)
- [DotMatch AssaySpec v1](docs/assayspec.md)
- [Evidence notes](docs/scientific-claims.md)
- [Public schemas](docs/schemas.md)
- [Changelog](CHANGELOG.md)
- [Release process](docs/release-process.md)
- [Contributing](CONTRIBUTING.md)

## Summary

DotMatch assigns short reads to a known library of short targets, quickly and with explicit ambiguity handling. For that class of problem, it gives scientists a smaller and clearer surface to trust.
