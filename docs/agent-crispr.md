# Run CRISPR guide counting with a local agent

This route prepares and executes deterministic fixed-window counting against a
known CRISPR guide library. It preserves ambiguous reads, writes
MAGeCK-compatible counts, evaluates guide-library and representation gates, and
can create a raw-data-free hashed handoff.

## Exact inputs

- `targets`: local TSV or CSV guide library with stable identifiers and guide sequences;
- `reads_dir`: local directory containing only the intended FASTQ or FASTQ.gz sample files;
- `output_dir`: absent or empty local directory;
- optional `threads`, `max_reads`, and `max_start` integers.

With `threads > 1`, the agent keeps aggregate assignment and ambiguity counts
in `summary.json` and `sample_qc.tsv` but disables incompatible ordered
row-level assignment, ambiguous-read, and unmatched-read files. Use one thread
when those per-read diagnostic files are required.

No field accepts a shell command. Paths are resolved locally. A non-empty
output directory, a symlinked output directory, missing input, unsafe target
library, or truncated FASTQ stops the workflow.

## Copyable start

Create `crispr-request.json`:

```json
{
  "intent": "crispr-guide-counting",
  "targets": "/absolute/path/guides.tsv",
  "reads_dir": "/absolute/path/fastqs",
  "output_dir": "/absolute/path/dotmatch-crispr-run",
  "threads": 4
}
```

Then run:

```bash
dotmatch agent invoke prepare_assay --input crispr-request.json
```

Use the returned `spec.path` in the next structured request:

```json
{"spec": "/absolute/path/dotmatch-crispr-run/assay.toml"}
```

```bash
dotmatch agent invoke inspect_assay --input assay-spec.json
dotmatch agent invoke run_assay --input assay-spec.json
dotmatch agent invoke review_assay --input assay-spec.json
```

For a passed revision, provide a new empty handoff directory:

```json
{
  "spec": "/absolute/path/dotmatch-crispr-run/assay.toml",
  "output_dir": "/absolute/path/dotmatch-crispr-handoff"
}
```

```bash
dotmatch agent invoke handoff_assay --input handoff.json
```

## Outputs

- candidate AssaySpec, inference report, candidates, and generated sample sheet;
- target collision audit and preflight reliability artifacts;
- MAGeCK-compatible count matrix, sample QC, CRISPR QC, normalized spec, methods, citation, and software-version records;
- one stable agent envelope with status/exit mapping, spec revision and hash, artifact hashes, normalized findings, resource use, and next actions;
- optional handoff manifest, review files, and SHA-256 list. Raw FASTQ is not copied.

## Automatic correction boundary

The run may create at most three numbered candidate specs for evidence-backed
changes to extraction start/length, read orientation, assignment metric,
correction radius, or CPU fallback. It never edits the original spec, changes
target sequences, counts ambiguous reads, relaxes QC thresholds, changes the
reliability profile, or performs downstream CRISPR statistics. Repeated states
and unresolved reliability blocks stop rather than widening the assay.

## Reference contract and scope

This route provides known-guide counting only: no downstream screen statistics.
Inspect the [checked contract fixture](https://dnncha.github.io/dotmatch/agent-reference-crispr.json)
when integrating the tool envelope. Its intentionally failed verdict exercises
unsafe and low-assignment states; it is not biological validation.
