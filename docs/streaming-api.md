# Python API

`dotmatch.stream_assign` is the Python API for notebook and workflow code that
needs row-level FASTQ assignment without loading a full run into memory. It uses
the same indexed native matcher and the same `unique`, `ambiguous`, `none`, and
`invalid` outcomes as the CLI.

```python
import dotmatch

rows = dotmatch.stream_assign(
    "reads.fastq.gz",
    "guides.tsv",
    target_start=23,
    target_length=20,
    k=1,
    policy="radius",
)

summary = dotmatch.write_assignments_tsv(rows, "assignments.tsv")
print(summary["assignment_rate"])
```

Targets can be a TSV/CSV path, a pandas or polars DataFrame accepted by
`targets_from_dataframe`, a list of `(target_id, sequence)` pairs, or a list of
sequences. FASTQ input may be plain text or `.gz`.

Each yielded `StreamAssignment` contains:

```text
read_id
observed_seq
target_index
target_name
target_seq
best_distance
second_best_distance
match_count
status
status_name
```

Only `unique` rows carry `target_name` and `target_seq`. Ambiguous reads are not
silently assigned.

## Helpers

- `dotmatch.load_targets(path)` loads TSV/CSV target tables into `(id, seq)`
  pairs.
- `dotmatch.iter_fastq(path)` yields validated FASTQ records from plain or
  gzipped files.
- `dotmatch.assignment_summary(rows)` returns counts and assignment,
  ambiguous, no-match, and invalid rates.
- `dotmatch.write_assignments_tsv(rows, path)` writes the stable assignment TSV
  shape and returns `assignment_summary`.

Use the CLI for production count matrices and reports. Use `stream_assign` when
Python code needs to inspect reads, join assignments with external metadata, or
feed a workflow-specific table while preserving DotMatch's ambiguity contract.
