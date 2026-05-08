# DotMatch Public Schemas

These are the open file contracts for DotMatch Core. They are intentionally plain TSV/JSON so workflow systems, MultiQC custom content, notebooks, and future workbench layers can consume them without linking to the C library.

## `target_counts.long.tsv`

One row per sample and target.

```text
sample_id
target_id
group
sequence
exact_count
k1_sub_count
k1_ins_count
k1_del_count
other_count
total_count
ambiguous_nearby
```

Rules:

- counts include only uniquely assigned reads;
- ambiguous reads are never added to a target count;
- `ambiguous_nearby=1` means another target can create ambiguity within the configured radius.

## `sample_qc.tsv`

One row per sample.

```text
sample_id
fastq
total_reads
valid_extracted_reads
assigned_reads
exact_reads
k1_rescued_reads
k1_sub_reads
k1_ins_reads
k1_del_reads
ambiguous_reads
no_match_reads
invalid_reads
assignment_rate
exact_rate
rescue_rate
ambiguous_rate
no_match_rate
targets_observed
zero_count_targets
gini_index
top_1pct_read_fraction
candidates_verified
```

Rules:

- rates are fractions from `0.0` to `1.0`;
- `valid_extracted_reads = total_reads - invalid_reads`;
- `assigned_corrected` is the preferred total for uniquely assigned non-exact reads;
- `k1_rescued_reads` is retained for compatibility and equals `assigned_corrected`,
  including in Levenshtein `k=2` runs.

## `pair_counts.tsv`

Nonzero paired/combinatorial target counts from `dotmatch pair-count`.

```text
left_id
right_id
count
```

Only reads with uniquely assigned left and right windows contribute to `count`.

## `pair_assignments.tsv`

Optional row-level diagnostics from `dotmatch pair-count --assignments`.

```text
read_id
left_observed
left_index
left_id
left_status
left_distance
right_observed
right_index
right_id
right_status
right_distance
pair_status
```

`pair_status` is `unique` only when both windows are uniquely assigned. If
either side is ambiguous, unmatched, or invalid, the read is excluded from
`pair_counts.tsv`.

## `pair_summary.json`

Top-level fields:

```text
workflow
k
metric
alphabet_policy
left_start
left_length
right_start
right_length
n_left_targets
n_right_targets
total_reads
assigned_pairs
pair_ambiguous
left_unmatched
right_unmatched
invalid
candidates_considered
candidates_verified
```

Rules:

- rates are fractions from `0.0` to `1.0`;
- `valid_extracted_reads = total_reads - invalid_reads`;
- `assigned_corrected` is the preferred total for uniquely assigned non-exact reads;
- `k1_rescued_reads` is retained for compatibility and equals `assigned_corrected`,
  including in Levenshtein `k=2` runs.

## `audit_summary.tsv`

Key-value summary of target-library safety.

```text
metric
value
```

Required metrics:

```text
audit_mode
targets
unique_sequences
duplicate_sequences
min_edit_distance
safe_at_k0
safe_at_k1
safe_at_k2
pairs_distance_0
pairs_distance_1
pairs_distance_2
pairs_within_requested_k
risk_pairs_for_k1
risk_pairs_for_k2
ambiguous_query_variants_k1
recommended_k
```

`audit_mode=exact` computes exhaustive pairwise distances. `audit_mode=fast` computes `k=1` safety through one-edit variant indexing and may report `not_computed` for `k=2` metrics.

## `audit_summary.json`

JSON equivalent of the audit summary for workflow engines and dashboards.

Fields:

```text
audit_mode
k
targets
unique_sequences
duplicate_sequences
min_edit_distance
safe_at_k0
safe_at_k1
safe_at_k2
pairs_distance_0
pairs_distance_1
pairs_distance_2
pairs_within_requested_k
risk_pairs_for_k1
risk_pairs_for_k2
ambiguous_query_variants_k1
recommended_k
```

Rules:

- safety fields are booleans when computed;
- `safe_at_k2` and `risk_pairs_for_k2` are `null` in fast audit mode;
- `min_edit_distance` is numeric in exact mode and may be the string `">=3"` in fast mode.

## `collision_pairs.tsv`

One row per target pair with collision risk.

```text
target_a
target_b
sequence_a
sequence_b
distance
risk_at_k1
risk_at_k2
example_ambiguous_query
```

## `target_safety.tsv`

One row per target.

```text
target_id
sequence
nearest_target
nearest_distance
safe_at_k1
safe_at_k2
num_nearby_k1_risk_targets
```

## `ambiguous_variants.tsv`

One row per query variant that would be within one edit of multiple targets.

```text
query_variant
targets_within_k1
```

This file answers the practical question behind one-edit rescue: which observed sequences would be ambiguous under exact `k=1` Levenshtein semantics?

## `top_unmatched.tsv`

One row per frequent unassigned extracted sequence.

```text
sequence
count
length
nearest_target
nearest_distance
nearest_edit_class
possible_reason
reverse_complement
revcomp_nearest_target
revcomp_nearest_distance
offset_hint
adapter_hint
```

Current reason labels:

```text
near_known_target_above_k
reverse_complement_candidate
offset_shift_candidate
adapter_or_primer_candidate
low_quality_candidate
contains_N
wrong_length
unknown
```

## `summary.json`

Run-level machine-readable summary. Top-level fields:

```text
k
metric
ambiguity_policy
alphabet_policy
max_correction_qual
indel_window
target_start
auto_offset
target_length
n_targets
samples
```

For count and demux summaries, `k=2` is currently a Levenshtein-only fixed-window
mode. Hamming summaries remain limited to `k=0` and `k=1`.

`alphabet_policy` records the assignment alphabet contract reported by
`qdaln_alphabet_policy()`: `N` and IUPAC ambiguity symbols are literal byte
symbols, not wildcard expansions. Demultiplexing summaries include the same
field.

`max_correction_qual` is either `null` or the Sanger Phred threshold supplied
with `--max-correction-qual`. When set, one-edit substitution and read-insertion
rescues require the observed edited base to have quality at or below this
threshold; exact matches and read-deletion rescues are not rejected by this
gate. Demultiplexing summaries include the same field.

Each sample object includes:

```text
sample
selected_target_start
total_reads
assigned_unique
assigned_exact
assigned_corrected
k1_rescued_reads
percent_rescued_by_k1
ambiguous
percent_ambiguous
unmatched
percent_unmatched
invalid
library_covered_targets
library_coverage_fraction
top_target_id
top_target_count
candidates_considered
candidates_verified
```
