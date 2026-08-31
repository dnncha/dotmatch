# Run Perturb-seq direct-guide capture with a local agent

This route assigns a reviewed fixed read window to a finite direct-guide list.
It is the same deterministic assignment contract used by CRISPR counting, with
the narrower scientific boundary required for Perturb-seq guide-capture reads.

## Exact inputs

- `targets`: local TSV or CSV containing the known guide-barcode identifiers and sequences;
- `reads_dir`: local directory containing only the intended guide-capture FASTQ or FASTQ.gz files;
- `output_dir`: absent or empty local directory;
- optional `threads`, `max_reads`, and `max_start` integers.

The target table must be known before assignment. DotMatch does not discover
new guide sequences from the reads.

## Copyable start

Create `perturb-seq-request.json`:

```json
{
  "intent": "perturb-seq-guide-capture",
  "targets": "/absolute/path/direct-guides.tsv",
  "reads_dir": "/absolute/path/guide-capture-fastqs",
  "output_dir": "/absolute/path/dotmatch-perturb-seq-run",
  "threads": 4
}
```

```bash
dotmatch agent invoke prepare_assay --input perturb-seq-request.json
```

Pass the returned spec through `inspect_assay`, `run_assay`,
`review_assay`, and—only when the reliability boundary permits it—
`handoff_assay`. The JSON shapes are identical to the
[CRISPR agent route](agent-crispr.md).

## Outputs

- inferred AssaySpec and inference evidence;
- per-guide unique counts and optional per-read assignments;
- separate ambiguous, unmatched, and invalid outcomes in QC and findings;
- reliability, provenance, artifact hashes, resource use, and raw-data-free handoff records.

## Evidence boundary

The GSE146194 reference evaluates a frozen 32-guide direct-capture rule on
48,000 held-out reads. DotMatch and matched independent per-read oracles have
zero differences at `k=0` and `k=1` for that recorded dataset, target list,
window, orientation, metric, and ambiguity rule.

This does not establish cell-barcode correction, UMI deduplication,
guide-per-cell calls, expression quantification, indirect-capture designs,
perturbation effects, or biological conclusions. Those steps and claims must
remain outside the DotMatch handoff.
