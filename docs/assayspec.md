# AssaySpec workflows

AssaySpec is a TOML workflow layer for fixed-window known-target assays. It
does not replace the native matching/counting code; it validates a declarative
spec, writes reproducible run files, prints the native commands it will execute,
and records an `assay_manifest.json` beside the outputs.

AssaySpec is part of the Python package surface. PyPI/source installs and the
Python-enabled Bioconda recipe expose `dotmatch assay ...` through the standard
`dotmatch` console script.

## Commands

```bash
dotmatch assay new crispr --library guides.csv --reads-dir fastqs/ --out crispr_screen/
dotmatch assay start assay.toml
dotmatch assay check assay.toml
dotmatch assay optimize assay.toml
dotmatch assay plan assay.toml
dotmatch assay run assay.toml
dotmatch assay init --template crispr --out assay.toml
dotmatch assay infer --mode count --assay-type crispr --targets guides.csv --reads sample.fastq.gz --out assay.toml --report inference_report.json
dotmatch assay autopsy assay.toml --out-dir autopsy/
dotmatch assay handoff assay.toml
```

`start` is the default production entrypoint: it runs `check`, then `run`, and
prints the reliability verdict plus paths to `reliability_report.html` and
`assay_fixes.tsv`. Exit codes follow the verdict: `0` for `passed`, `1` for
`needs_review`, and `2` for `failed` or `blocked`.

Templates are:

- `crispr`
- `feature-barcode`
- `inline-barcode-count`
- `inline-barcode-demux`
- `amplicon-panel`
- `oligo-adapter`
- `pair-count`

`optimize` reads the spec and writes a benchmark-informed CPU/GPU backend
recommendation without changing the assignment authority. `plan` is a dry run:
it prints deterministic native commands and does not create the output
directory. The plan also lists the reliability artifacts that a run or check
will write. `check` validates the spec and writes preflight reliability
artifacts without running read assignment. `run` creates the output directory,
writes generated files, runs target audit first, runs the compiled native
workflow, and records command exit codes and warnings in `assay_manifest.json`.
It also writes `assay_report.html` as the primary workflow report and
`assay_manifest.summary.tsv` for workflow systems and MultiQC custom content.
Production runs also write `methods.md`, `CITATION.bib`, and
`software_versions.yml` for lab notebooks, workflow submissions, and report
handoff.

`new` scaffolds a reviewable project directory from a target or barcode table and
a directory of FASTQ files. It stages inputs under `inputs/` and `reads/`,
infers the extract window from the first sample, writes `assay.toml`,
`inference_report.json`, `samples.generated.tsv`, `README.md`, and `run.sh`.
Low-confidence inference writes `status = "draft"` until a scientist promotes
the spec to `ready`.

`infer` samples FASTQ reads, scores fixed-window candidates against the supplied
target table, writes a candidate AssaySpec, and writes `inference_report.json`
plus `inference_candidates.tsv`. Low-confidence inference writes
`status = "draft"`. Production runs refuse draft specs by default until a user
reviews the report and changes the status to `ready`; this can be explicitly
relaxed with `[reliability] fail_on_draft_inference = false`.

`autopsy` helps diagnose suspicious runs by wrapping native target audit and
`inspect-unmatched`. It writes `autopsy_summary.json`, `findings.tsv`, and
`top_unmatched.*.tsv` files. `run` also triggers autopsy automatically when
sample QC crosses conservative thresholds.

`handoff` creates `<out_dir>/handoff/` after a completed run. It copies a
review-sized set of reports, configuration, QC, count, methods, and citation
files; writes `SHA256SUMS`; and records SHA-256 hashes for every declared input.
It does not copy FASTQ files, so raw reads remain in the originating controlled
workspace. Use `--out-dir PATH` to write the bundle to an empty directory.

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
ambiguity_policy = "radius"
ambiguous = "discard"

[reliability]
profile = "production"
fail_on_unsafe_targets = true
fail_on_draft_inference = true
min_assignment_rate = 0.80
max_ambiguous_rate = 0.05
max_unmatched_rate = 0.15
max_invalid_rate = 0.02
require_public_evidence_boundary = true

[backend]
mode = "auto"
allow_gpu = true

[outputs]
format = "mageck"
assignments = true
ambiguous = true
unmatched = true
```

Count mode writes `counts.mageck.tsv` for CRISPR/MAGeCK output or `counts.tsv`
for DotMatch output, plus `target_counts.long.tsv`, `sample_qc.tsv`,
`summary.json`, native `report.html`, `assay_report.html`,
`assay_manifest.json`, `assay_manifest.summary.tsv`, `methods.md`,
`CITATION.bib`, `software_versions.yml`, `audit/`, and optional row-level
diagnostics. CRISPR count runs also write `crispr_qc.json`,
`crispr_qc.summary.tsv`, and `crispr_qc.html`.

## Reliability Artifacts

AssaySpec writes closed-loop reliability artifacts for `check` and `run`:

- `reliability_summary.json`
- `reliability_findings.tsv`
- `reliability_report.html`
- `reliability_manifest.summary.tsv`

`reliability_summary.json` records the profile, thresholds, backend authority,
GPU eligibility, evidence boundary, artifact paths, and normalized findings.
Findings use `info`, `warning`, `error`, and `blocked` severity. `check` records
read-dependent QC as unavailable because no assignment has run yet. `run`
aggregates target audit results, sample QC thresholds, autopsy findings, command
failures, and assay evidence metadata.

`profile = "production"` fails fast after target audit if
`fail_on_unsafe_targets = true` and the target set is unsafe at the configured
radius. Runtime threshold failures are written after assignment and leave output
artifacts in place for diagnosis. `profile = "exploratory"` records the same
conditions as findings without using unsafe preflight status to stop the run.

`[backend] mode = "auto"` keeps CPU assignment as the production authority and
records whether the assay is eligible for the experimental Metal GPU path.
`gpu-metal-experimental` remains advisory unless an assay-specific real-workload
gate validates it; DotMatch does not silently switch production assignment to
GPU.

## Demux And Pair Modes

Demux mode uses `mode = "demux"`, `barcodes`, `reads`, `[extract]`, and writes
`demuxed/`, `summary.json`, optional `assignments.tsv`, `ambiguous.fastq`, and
`unmatched.fastq`.

Pair mode uses `mode = "pair-count"`, `left_targets`, `right_targets`, `[left]`, and `[right]`.
It accepts one of two input layouts:

```toml
# Both target windows occur in one read.
reads = "reads.fastq.gz"

# Or, remove reads and use synchronized mates.
left_reads = "sample_R1.fastq.gz"
right_reads = "sample_R2.fastq.gz"
```

For paired inputs, left extraction coordinates apply to `left_reads` and right
coordinates apply to `right_reads`. The two files must have the same number of
complete records in matching order. DotMatch compares canonical read IDs and
stops on the first mismatch. It writes `pair_counts.tsv`, `pair_summary.json`, and
optional `pair_assignments.tsv`. The summary records the input layout and
side-specific unmatched and invalid outcomes; methods and assay reports list
both FASTQ inputs.

## Safety Policy

AssaySpec always runs native target audit before assignment. If the audit says
the target set is unsafe at the configured `k`, the reliability profile decides
whether to stop before assignment or record the risk and continue. It never
changes `k`, target sequences, or ambiguity policy automatically; DotMatch's
explicit `unique`/`ambiguous`/`none` semantics remain the authority.

Templates and inferred specs default to `ambiguity_policy = "radius"`, which
keeps any read with more than one target inside the configured radius out of
forced assignments. Use `ambiguity_policy = "best"` only when best-distance
compatibility is deliberate.

Automatic autopsy uses the configured reliability thresholds for assignment,
ambiguous, no-match, and invalid rates. These thresholds are recorded in
`assay_manifest.json`.
