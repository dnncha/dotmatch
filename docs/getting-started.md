# Getting started

This guide walks from installation to a verified DotMatch run. It assumes you have FASTQ or FASTQ.gz reads and a table of expected short DNA targets such as CRISPR guides, sample barcodes, feature barcodes, primers, or adapter prefixes.

## Install

From PyPI:

```bash
python3 -m pip install dotmatch
dotmatch --version
```

From Bioconda:

```bash
conda create -n dotmatch -c conda-forge -c bioconda dotmatch
conda activate dotmatch
dotmatch --version
```

From a source checkout (requires a C compiler, `make`, Python 3.9+, and zlib for FASTQ.gz support):

```bash
git clone https://github.com/dnncha/dotmatch.git
cd dotmatch
make
python3 -m pip install .
dotmatch --version
```

PyPI normally receives new releases first. Bioconda can lag a release tag while its recipe builds, so check `dotmatch --version` when an exact release matters.

## Verify your install offline

Run the packaged synthetic check without using study data or making network requests:

```bash
dotmatch demo --out-dir first-run/
```

Open `first-run/index.html` to inspect the generated review bundle.

## Count a CRISPR screen

`crispr quickstart` turns one or more FASTQ paths and a guide table into a self-contained AssaySpec project. It infers the guide window, writes an inference report, and leaves the project in `draft` status for review:

```bash
dotmatch crispr quickstart \
  --library guides.csv \
  --fastq 'fastqs/*.fastq.gz' \
  --out crispr_screen/
```

Inspect `crispr_screen/inference_report.json` and `crispr_screen/assay.toml` to confirm the read window, orientation, library, and sample files. Change `status = "draft"` to `status = "ready"` in `assay.toml`, then run the count and build a raw-FASTQ-free handoff bundle:

```bash
dotmatch assay start crispr_screen/assay.toml

# Package config, QC, methods, and checksums without raw FASTQs:
dotmatch assay handoff crispr_screen/assay.toml
```

To skip the draft pause on an already-known layout, pass `--accept-inference` to `crispr quickstart` (or `--no-run` for a review-only setup). A non-zero exit code after execution means a reliability gate needs inspection: open `assay_out/reliability_report.html` and apply `assay_out/assay_fixes.tsv` before using counts downstream.

When starting from an entire directory of FASTQs, `dotmatch assay new crispr --library guides.tsv --reads-dir fastqs/ --out crispr_screen/` (or `dotmatch crispr new`) scaffolds the same project layout with a `./run.sh` wrapper.

## Prepare and audit a target table

Targets are tab-separated records with an identifier and a sequence: one expected guide, barcode, primer, feature tag, or panel target per row.

```text
target_id	sequence
guide_001	ACGTACGTACGTACGTACGT
guide_002	ACGTACGTACGTACGTAGGT
guide_003	TGCATGCATGCATGCATGCA
```

Before enabling mismatch correction (`--k 1` or higher), audit the target set for duplicates and near neighbors:

```bash
dotmatch audit \
  --targets guides.tsv \
  --k 1 \
  --audit-mode auto \
  --out-dir audit/
```

Review duplicate sequences, near-neighbor pairs, and ambiguous example queries before running correction in production. Count unsafe target sets exactly (`--k 0`) or redesign the panel.

## Count known targets directly from FASTQ

Use `dotmatch count` when reads contain one fixed target window:

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

Use `--metric hamming` when all targets and read windows share the same length and only substitutions should be corrected. Use `--metric levenshtein --indel-window 1` when one-base insertions or deletions should be considered:

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

The default `radius` ambiguity policy counts a read only when a single target lies inside the configured distance radius. Pass `--ambiguity-policy best` only when you intentionally want nearest-target compatibility with another tool.

## Demultiplex inline barcodes

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

Inspect `demux.summary.json` and `demux.assignments.tsv` before using the split FASTQs. High unmatched, ambiguous, or invalid rates usually point to a wrong barcode offset, length mismatch, sample-sheet error, or overly permissive correction radius.

## Diagnose a barcode run

`dotmatch barcode autopsy` diagnoses common demultiplexing failures in one pass:

```bash
dotmatch barcode autopsy \
  --barcodes barcodes.tsv \
  --reads pooled.fastq.gz \
  --scan-starts 0:12 \
  --k-values 0,1 \
  --out-dir autopsy/
```

Open `autopsy/report.html` first, and archive the accompanying TSV and JSON tables for pipeline logs or lab handoff.

## Read the outputs

- `counts.tsv` or `target_counts.long.tsv`: counts for uniquely assigned reads.
- `sample_qc.tsv`: assignment rate, rescue rate, ambiguous reads, unmatched reads, invalid windows, target coverage, and representation metrics.
- `assignments.tsv`: per-read assignment states when requested.
- `summary.json`: run configuration, assignment policy, and provenance.
- HTML reports: self-contained review pages for assay, barcode, panel, and QC workflows.

Ambiguous reads stay separate and never contribute silently to target counts under the default policy.

## Next steps

- [AssaySpec](assayspec.md): full reference for `assay new`, `start`, `check`, and `run`.
- [Lab evaluation and handoff](lab-evaluation.md): create portable review bundles without copying raw reads.
- [Command reference](command-reference.md): complete CLI command map.
- [CRISPR count QC](crispr-qc.md): checks to review before downstream screen statistics.
- [Barcode panel design](barcode-panel-design.md): design, check, and export barcode sets.
- [Public schemas](schemas.md): output file contracts for Snakemake, Nextflow, MultiQC, and LIMS exports.
- [Methods and citation](methods-and-citation.md): release citation guidance, `CITATION.cff`, and `dotmatch citation`.
