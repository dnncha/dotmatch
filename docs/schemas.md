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
- `k1_rescued_reads` includes all uniquely assigned non-exact reads.

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
indel_window
target_start
auto_offset
target_length
n_targets
samples
```

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
