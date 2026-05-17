# DotMatch

[![CI](https://github.com/dnncha/dotmatch/actions/workflows/ci.yml/badge.svg)](https://github.com/dnncha/dotmatch/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![Citation](https://img.shields.io/badge/cite-CITATION.cff-green.svg)](CITATION.cff)

DotMatch assigns short sequencing reads to a known set of short DNA targets.
It is intended for CRISPR guide counting, inline barcode demultiplexing,
fixed-window feature/barcode assignment, primer or adapter-prefix checks, and
similar workflows where the reference is a target library rather than a genome.

For each read, DotMatch reports one of four outcomes:

- `unique`: exactly one target is the best match under the configured edit radius;
- `ambiguous`: more than one target is compatible with the read;
- `none`: no target is within the configured radius;
- `invalid`: the configured sequence window cannot be extracted.

Ambiguous reads are exposed explicitly instead of being assigned silently. That
is the main design choice: DotMatch is for deterministic known-target
assignment, not probabilistic read mapping.

## When To Use DotMatch

DotMatch is a good fit when you have a table of expected short sequences and
want reproducible per-read assignment or count tables.

Common uses include:

- CRISPR pooled-screen guide counting with MAGeCK-compatible output;
- fixed-position barcode demultiplexing from FASTQ/FASTQ.gz;
- per-read assignment of 10x guide-capture or feature-barcode windows;
- primer-start, amplicon-panel, adapter-prefix, or whitelist-style assays;
- target-library audits before allowing one-edit correction;
- validating an indexed assignment run against an exhaustive scan or Edlib.

DotMatch is not a genome aligner. It does not produce SAM/BAM, CIGAR strings,
variant calls, cell/UMI quantification, expression matrices, or guide-level
statistics. It works on extracted short windows and known target lists.

## Installation

DotMatch currently supports source builds and local Python package installs on
Linux and macOS. You need a C compiler, `make`, Python 3.9 or newer for the
Python package, and zlib for FASTQ.gz support.

```bash
git clone https://github.com/dnncha/dotmatch.git
cd dotmatch
make

./dotmatch --version
./dotmatch dist ACGT AGGT
./dotmatch leq 1 ACGT AGGT
```

Python install from a checkout:

```bash
python3 -m pip install .
python3 -c "import dotmatch; print(dotmatch.distance('ACGT', 'AGGT'))"
```

Docker build from the repository:

```bash
docker build -t dotmatch:dev .
docker run --rm -v "$PWD:/work" dotmatch:dev dist ACGT AGGT
```

Package-channel status for PyPI, Bioconda, containers, and release archives is
tracked in [Packaging Notes](docs/packaging.md).

## Quick Example

The core operation is many-read versus many-target assignment. Target files and
read files can be simple TSVs with `id<TAB>sequence`.

```bash
cat > targets.tsv <<'EOF'
bc0	ACGT
bc1	AGGT
bc2	ACGA
EOF

cat > reads.tsv <<'EOF'
r0	ACGT
r1	ACGC
r2	TTTT
EOF

./dotmatch assign 1 targets.tsv reads.tsv
```

Expected output:

```text
mode	read_id	read_seq	target_index	target_seq	distance	status	match_count	second_best_distance
assign	r0	ACGT	0	ACGT	0	unique	3	1
assign	r1	ACGC	0	ACGT	1	ambiguous	2	-1
assign	r2	TTTT	-1		-1	none	0	-1
```

`r1` is deliberately ambiguous: it is within one edit of more than one target,
so DotMatch reports the ambiguity instead of choosing a target.

## CRISPR Guide Counting

For pooled CRISPR screens, `crispr-count` wraps the FASTQ counting engine and
writes a MAGeCK-style count matrix.

```bash
cat > samples.tsv <<'EOF'
sample_id	fastq
plasmid	plasmid_R1.fastq.gz
treatment	treatment_R1.fastq.gz
EOF

./dotmatch crispr-count \
  --library guides.csv \
  --samples samples.tsv \
  --guide-start 23 \
  --guide-length 20 \
  --k 1 \
  --metric hamming \
  --out counts.mageck.tsv \
  --summary qc.json \
  --ambiguous discard
```

Use `--metric hamming` for one-mismatch/no-indel guide-counter-style counting.
Use `--metric levenshtein --indel-window 1` when one-base insertions and
deletions around the guide window should be considered. Ambiguous reads are not
added to guide counts unless you explicitly request diagnostic reporting.

A small worked example is available in
[examples/crispr_guides](examples/crispr_guides/README.md), and a step-by-step
fixture walkthrough is in
[docs/tutorials/crispr-count-first-run.md](docs/tutorials/crispr-count-first-run.md).

## General FASTQ Counting

The lower-level `count` command works with arbitrary known targets and one or
more FASTQ/FASTQ.gz inputs.

```bash
./dotmatch count \
  --targets targets.tsv \
  --reads sample_R1.fastq.gz \
  --sample-label sample_1 \
  --target-start 0 \
  --target-length 20 \
  --k 1 \
  --metric levenshtein \
  --indel-window 1 \
  --ambiguity-policy best \
  --out counts.tsv \
  --target-counts-long target_counts.long.tsv \
  --sample-qc sample_qc.tsv \
  --assignments assignments.tsv \
  --summary summary.json
```

The count table separates exact matches, one-substitution corrections,
one-insertion corrections, one-deletion corrections, and other accepted
corrections. `sample_qc.tsv` records assignment rate, rescue rate, ambiguous and
unmatched fractions, target coverage, zero-count targets, Gini index, and the
number of candidate targets checked after indexing.

Output schemas are documented in [Public Schemas](docs/schemas.md).

## Barcode Demultiplexing

For fixed-position inline barcodes, `demux` writes one FASTQ per uniquely
assigned barcode and can optionally retain ambiguous and unmatched reads.

```bash
./dotmatch demux \
  --barcodes barcodes.tsv \
  --reads pooled.fastq.gz \
  --barcode-start 0 \
  --barcode-length 8 \
  --k 1 \
  --metric hamming \
  --max-correction-qual 20 \
  --out-dir demuxed \
  --summary demux.qc.json \
  --assignments demux.assignments.tsv \
  --ambiguous-out ambiguous.fastq \
  --unmatched-out unmatched.fastq
```

Use `--barcode-length auto` when the barcode sheet contains multiple lengths.
Prefix-overlapping exact matches are reported as ambiguous rather than resolved
by length.

DotMatch also includes an early classic per-cycle BCL demultiplexing command for
small RunInfo/SampleSheet/BCL workflows. CBCL and NovaSeq-style inputs are not
part of the current BCL scope.

## Target Library Audit

Before enabling one-edit correction, audit the target set for collisions and
near neighbors.

```bash
./dotmatch audit \
  --targets guides.tsv \
  --k 1 \
  --audit-mode auto \
  --out-dir audit/
```

The audit output includes duplicate targets, nearby target pairs, collision
clusters, per-target safety, and example query variants that would be ambiguous
at `k=1`.

## Python API

The Python package loads the native library through `ctypes`.

```python
import dotmatch

dotmatch.distance("ACGT", "AGGT")
# 1

dotmatch.distance_leq("ACGT", "AGGT", 1)
# True

matcher = dotmatch.Matcher(["ACGT", "AGGT", "ACGA"])
results, stats = matcher.assign_with_stats(["ACGT", "ACGC"], k=1)
```

When working from a source checkout, build the shared library first:

```bash
make shared
DOTMATCH_LIB=$PWD/libdotmatch.dylib PYTHONPATH=$PWD/python python3
```

On Linux, use `libdotmatch.so` instead of `libdotmatch.dylib`.

The historical `quickdna` Python package, `quickdna` console script, and `qda`
native CLI target remain as compatibility aliases. New workflows should use
`dotmatch`.

## Matching Semantics

DotMatch uses literal-byte DNA matching. `A`, `C`, `G`, `T`, `N`, and IUPAC
ambiguity symbols are ordinary byte symbols; `N` and IUPAC codes are not
expanded as wildcards.

Supported assignment modes include:

- exact matching (`k=0`);
- Hamming matching for fixed-length one-substitution workflows;
- global Levenshtein matching for substitutions, insertions, and deletions;
- fixed-window `k=2` correction through the exhaustive assignment path;
- explicit ambiguity policies for best-target and whole-radius assignment.

The public policy string reported by the C and Python APIs is:

```text
literal-byte; A/C/G/T/N/IUPAC symbols are ordinary byte symbols; no wildcard expansion
```

## Validation And Benchmarks

The repository includes native C tests, CLI fixture tests, Python tests,
deterministic fuzz checks against a dynamic-programming oracle, and optional
Edlib validation for assignment runs.

Useful local checks:

```bash
make test
make cli-test
make python-test
make python-package-test
```

Benchmark and evidence reports:

- [Benchmark overview](docs/benchmarks/README.md)
- [Native Edlib assignment report](docs/benchmarks/native/README.md)
- [Public CRISPR guide-counting report](docs/benchmarks/public_crispr/README.md)
- [Barcode demultiplexing report](docs/benchmarks/barcode_demux/README.md)
- [Feature-barcode assignment report](docs/benchmarks/feature_barcode/README.md)
- [CRISPR guide-capture assignment report](docs/benchmarks/perturb_seq/README.md)
- [Amplicon/panel primer-start report](docs/benchmarks/amplicon_panel/README.md)
- [Oligo/adapter prefix-assignment report](docs/benchmarks/oligo_adapter/README.md)

For the current limits of public claims, see
[Evidence Notes](docs/scientific-claims.md). For methods text and citation
language, see [Methods and Citation](docs/methods-and-citation.md).

## Development

```bash
make
make test
make cli-test
make coverage
```

Workflow-manager examples are included for Galaxy, Nextflow, nf-core-style
modules, Snakemake, and MultiQC custom content under
[examples/workflows](examples/workflows/).

Contributions are welcome. Please read [CONTRIBUTING.md](CONTRIBUTING.md),
[SUPPORT.md](SUPPORT.md), and [SECURITY.md](SECURITY.md) before opening issues
or pull requests.

## Citation

If DotMatch is useful in your work, cite the software release using
[CITATION.cff](CITATION.cff). Suggested methods wording is provided in
[docs/methods-and-citation.md](docs/methods-and-citation.md).

## License

DotMatch is released under the [Apache License 2.0](LICENSE).
