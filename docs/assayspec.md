# DotMatch AssaySpec v1

AssaySpec is a TOML workflow layer for fixed-window known-target assays. It
does not replace the native matching/counting code; it validates a declarative
spec, writes reproducible run files, prints the native commands it will execute,
and records an `assay_manifest.json` beside the outputs.

## Commands

```bash
dotmatch assay init --template crispr --out assay.toml
dotmatch assay infer --mode count --assay-type crispr --targets guides.csv --reads sample.fastq.gz --out assay.toml --report inference_report.json
dotmatch assay check assay.toml
dotmatch assay plan assay.toml
dotmatch assay run assay.toml
dotmatch assay autopsy assay.toml --out-dir autopsy/
```

Templates are:

- `crispr`
- `feature-barcode`
- `inline-barcode-count`
- `inline-barcode-demux`
- `amplicon-panel`
- `oligo-adapter`
- `pair-count`

`plan` is a dry run: it prints deterministic native commands and does not create
the output directory. `run` creates the output directory, writes generated files,
runs target audit first, runs the compiled native workflow, and records command
exit codes and warnings in `assay_manifest.json`.

`infer` samples FASTQ reads, scores fixed-window candidates against the supplied
target table, writes a candidate AssaySpec, and writes `inference_report.json`
plus `inference_candidates.tsv`. Low-confidence inference writes
`status = "draft"`. `run` refuses draft specs until a user reviews the report
and changes the status to `ready`.

`autopsy` diagnoses suspicious runs by wrapping native target audit and
`inspect-unmatched`. It writes `autopsy_summary.json`, `findings.tsv`, and
`top_unmatched.*.tsv` files. `run` also triggers autopsy automatically when
sample QC crosses conservative thresholds.

## Count Example

```toml
schema_version = 1
status = "ready"
mode = "count"
assay_type = "crispr"
targets = "guides.csv"

[[samples]]
id = "control"
fastq = "control.fastq.gz"

[[samples]]
id = "treated"
fastq = "treated.fastq.gz"

[run]
out_dir = "dotmatch_assay_out"
threads = 1

[extract]
start = 23
length = 19

[assignment]
k = 1
metric = "hamming"
ambiguous = "discard"

[outputs]
format = "mageck"
assignments = true
ambiguous = true
unmatched = true
```

Count mode writes `counts.mageck.tsv` for CRISPR/MAGeCK output or `counts.tsv`
for DotMatch output, plus `target_counts.long.tsv`, `sample_qc.tsv`,
`summary.json`, `report.html`, `audit/`, and optional row-level diagnostics.

## Demux And Pair Modes

Demux mode uses `mode = "demux"`, `barcodes`, `reads`, `[extract]`, and writes
`demuxed/`, `summary.json`, optional `assignments.tsv`, `ambiguous.fastq`, and
`unmatched.fastq`.

Pair mode uses `mode = "pair-count"`, `left_targets`, `right_targets`, `reads`,
`[left]`, and `[right]`. It writes `pair_counts.tsv`, `pair_summary.json`, and
optional `pair_assignments.tsv`.

## Safety Policy

AssaySpec always runs native target audit before assignment. If the audit says
the target set is unsafe at the configured `k`, `dotmatch assay run` records a
warning and continues. It never changes `k`, target sequences, or ambiguity
policy automatically; DotMatch's explicit `unique`/`ambiguous`/`none` semantics
remain the authority.

Automatic autopsy triggers when any sample has assignment rate below `0.80`,
ambiguous rate above `0.05`, no-match rate above `0.15`, or invalid rate above
`0.02`. These thresholds are recorded in `assay_manifest.json`.
