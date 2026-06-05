# Getting Started

This page is the fastest useful path from installation to a checked DotMatch
run. It assumes you have FASTQ or FASTQ.gz reads and a table of expected short
DNA targets such as guides, sample barcodes, feature barcodes, primers, or
adapter prefixes.

## Install

For the current PyPI release:

```bash
python3 -m pip install dotmatch==0.1.7
dotmatch --version
```

For the currently published Bioconda package:

```bash
conda create -n dotmatch -c conda-forge -c bioconda dotmatch=0.1.4
conda activate dotmatch
dotmatch --version
```

From a source checkout:

```bash
git clone https://github.com/dnncha/dotmatch.git
cd dotmatch
make
python3 -m pip install .
dotmatch --version
```

The source build needs a C compiler, `make`, Python 3.9 or newer, and zlib for
FASTQ.gz support.

## Prepare Targets

Targets are ordinary tabular records with an identifier and sequence. Keep the
file small and explicit: one expected guide, barcode, primer, feature tag, or
panel target per row.

```text
target_id	sequence
guide_001	ACGTACGTACGTACGTACGT
guide_002	ACGTACGTACGTACGTAGGT
guide_003	TGCATGCATGCATGCATGCA
```

Before allowing error correction, audit the target table:

```bash
dotmatch audit \
  --targets guides.tsv \
  --k 1 \
  --audit-mode auto \
  --out-dir audit/
```

Review duplicate sequences, near-neighbor targets, and ambiguous example
queries before using `--k 1` or higher in production. A target set that is not
safe for correction should usually be counted exactly or redesigned.

## Count Known Targets From FASTQ

Use `dotmatch count` when reads contain one fixed target window.

```bash
dotmatch count \
  --targets guides.tsv \
  --reads sample_R1.fastq.gz \
  --sample-label sample_1 \
  --target-start 23 \
  --target-length 20 \
  --k 1 \
  --metric hamming \
  --ambiguity-policy radius \
  --out counts.tsv \
  --target-counts-long target_counts.long.tsv \
  --sample-qc sample_qc.tsv \
  --assignments assignments.tsv \
  --summary summary.json
```

Use Hamming distance when all targets and read windows have the same length and
only substitutions should be corrected. Use Levenshtein distance when one-base
insertions or deletions should be considered:

```bash
dotmatch count \
  --targets targets.tsv \
  --reads sample_R1.fastq.gz \
  --target-start 0 \
  --target-length 20 \
  --k 1 \
  --metric levenshtein \
  --indel-window 1 \
  --ambiguity-policy radius \
  --out counts.tsv \
  --sample-qc sample_qc.tsv \
  --summary summary.json
```

The default `radius` ambiguity policy is conservative: a read is counted only
when exactly one target lies inside the configured radius. Use `best` only when
you deliberately need best-distance compatibility with an existing workflow.

## Demultiplex Inline Barcodes

For fixed-position inline barcodes:

```bash
dotmatch demux \
  --barcodes barcodes.tsv \
  --reads pooled.fastq.gz \
  --barcode-start 0 \
  --barcode-length 8 \
  --k 1 \
  --metric hamming \
  --ambiguity-policy radius \
  --max-correction-qual 20 \
  --out-dir demuxed \
  --summary demux.summary.json \
  --assignments demux.assignments.tsv \
  --ambiguous-out ambiguous.fastq \
  --unmatched-out unmatched.fastq
```

Open the summary and assignment tables before trusting the split FASTQs. High
unmatched, ambiguous, or invalid rates usually mean the barcode start, barcode
length, sample sheet, or correction policy needs review.

## Diagnose a Barcode Run

`barcode autopsy` is the quickest way to inspect common failure modes:

```bash
dotmatch barcode autopsy \
  --barcodes barcodes.tsv \
  --reads pooled.fastq.gz \
  --scan-starts 0:12 \
  --k-values 0,1 \
  --out-dir autopsy/
```

Start with `autopsy/report.html`, then use the TSV and JSON files for pipeline
records, lab handoff, or methods review.

## Read the Outputs

The most important files are:

- `counts.tsv` or `target_counts.long.tsv`: counts for uniquely assigned reads.
- `sample_qc.tsv`: assignment rate, rescue rate, ambiguous reads, unmatched
  reads, invalid windows, target coverage, and representation metrics.
- `assignments.tsv`: per-read assignment states when requested.
- `summary.json`: run configuration, assignment policy, and provenance.
- HTML reports: human-readable review pages for assay, barcode, panel, and QC
  workflows.

Ambiguous reads are intentionally visible and are not added silently to target
counts under the default policy. This is the central trust contract of
DotMatch.

## Next Steps

- Use [AssaySpec](assayspec.md) for declarative, reproducible workflow runs.
- Use [CRISPR Count QC](crispr-qc.md) before downstream screen statistics.
- Use [Barcode Panel Design](barcode-panel-design.md) when creating or checking
  barcode panels.
- Use [Public Schemas](schemas.md) when integrating DotMatch with Snakemake,
  Nextflow, MultiQC, notebooks, or LIMS exports.
- Use [Methods and Citation](methods-and-citation.md) and `dotmatch citation`
  when recording the software version in methods sections, reports, or workflow
  provenance. Release citation metadata is kept in `CITATION.cff`.
