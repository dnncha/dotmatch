# Sequence duplication audit

A thousand identical guide reads may represent a thousand observations of the same guide. Removing them because their sequences match would change the abundance being measured. This audit measures repetition while preserving every input read.

```sh
dotmatch duplication --reads sample.fastq.gz
dotmatch duplication --reads R1.fastq.gz --reads2 R2.fastq.gz --json > duplication.json
```

The command is available through the installed Python `dotmatch` entry point (or `python -m dotmatch.cli`), not the standalone native executable. No new third-party dependency is required.

For a quick look, add `--max-records 100000`. This examines an input prefix, not a random sample; it must not be extrapolated to the full library. For a full audit, omit that option. Use `--temp-dir /path/to/scratch` to select scratch storage and `--cache-mb 64` to configure SQLite's page cache. The cache setting is not a hard process memory bound. Disk usage grows with distinct sequences and database overhead; no fixed input-to-disk multiplier is promised. Scratch files are removed on success and handled errors. Abrupt process termination can leave scratch files behind.

Four-line FASTQ and gzip FASTQ are supported. Sequences must contain A/C/G/T/N, case is normalized, N is literal, and quality values do not define identity. Different lengths and orientations remain distinct. Paired identity requires both full sequences to match; pair identifiers (ignoring terminal /1 and /2) and record counts must agree. Pair checks cover only the audited prefix when a limit is set. Wrapped FASTQ is rejected.

## Agent contract

`--json` writes one JSON object to stdout. Successful audits exit 0; input or storage failures exit 2 with `status: "error"` and a message. Argument syntax errors use argparse's stderr and exit 2. Agents must check both exit status and `status` before reading metrics.

- `schema_version`: `dotmatch.duplication.v1`.
- `scope`: `full` or `prefix`; `unit`: `reads` or `read_pairs`.
- `records`: observations examined; `distinct_sequences`: exact sequence groups (full ordered mate tuples for paired data).
- `repeated_records`: records minus distinct groups, i.e. observations beyond the first in each group.
- `repeated_fraction`: repeated records divided by records; null for empty input.
- `singleton_sequences`: groups observed once; `largest_sequence_group`: largest multiplicity, or zero for empty input.
- `pcr_duplicate_fraction`: always null; this audit cannot establish PCR origin.
- `recommended_action`: `preserve_read_counts`; `input_modified`: false.

Do not interpret a high repetition fraction as a failed assay or automatically deduplicate a subsequent count job. The report intentionally sets no universal pass/fail threshold. No sequence or read identifier is included in the output.

Python callers can use `from dotmatch.duplication import audit_duplication`; this returns the same report dictionary and raises exceptions on failure.

## Research decisions

[Fastq-dupaway (Sigorskikh et al., 2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12749013/) motivates moving exact sequence grouping onto disk when RAM is constrained. Dotmatch uses SQLite full-string keys rather than implementing that paper's external sort. This is an independent audit implementation, not a reproduction of its performance results. The paper explicitly lacks ground-truth false-positive/negative evaluation and discusses loss of biological signal in low-diversity data. Its approximate tail-Hamming and prefix-based removal modes were not adopted: their grouping semantics would be inappropriate defaults for known-target abundance.

[UMI-tools (Smith et al., 2017)](https://genome.cshlp.org/content/27/3/491) supports using molecular identifiers with assay context and error modeling when molecular counts are required. A generic UMI deduplicator is not added here: the correct grouping context depends on the assay. This audit does not infer molecules from sequence identity.

Validation covers exact identity, paired identity, input preservation, malformed and mismatched data, gzip, empty inputs, prefix scope, order invariance, quality independence, and temporary database cleanup. It does not establish large-dataset speed or memory performance.
