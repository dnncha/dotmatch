# Data contract

## Manifest

Generate a complete example with `editwitness demo`; inspect
`editwitness schema manifest` or the committed schema snapshots under
`src/editwitness/schemas/`. Runtime validation is authoritative.

JSON Schema describes structure. It does **not** encode every cross-field rule:
reference bounds, primer matches, edit nonoverlap, known IDs, or work budgets.
Always run `editwitness validate` or construct `Manifest` through Pydantic.

| Field | Meaning |
|---|---|
| `schema_version` | `1.0` and `1.1` inputs are accepted; the current default is `1.1`. |
| `observation_model` | `original-sites-presence-v1` if omitted; select `exact-local-sites-presence-v2` explicitly for edited-sequence rematching. |
| `generation` | Optional declared deletion-generation provenance; not an authenticated attestation. |
| `coordinate_system` | `0-based-half-open` (default; no automatic coordinate conversion). |
| `reference` | Name, uppercase A/C/G/T local sequence, optional assembly/contig/genomic start, synthetic flag. |
| `alleles` | Unique IDs, optional descriptions, sorted local replacements. |
| `hypotheses` | Unique IDs and exactly two known allele IDs. Repetition represents two modeled copies. |
| `expected_hypothesis` | ID of one declared hypothesis. |
| `assays` | One or more existing assays. |
| `candidates` | Optional follow-up assays with positive integer cost units. |
| `deletion_scan` | Optional single-deletion endpoint grid. |

IDs use letters, numbers, underscores, dots and hyphens; they start with a letter
or digit and have at most 80 characters. Assay IDs are unique across both existing
and candidate collections. Unknown fields are rejected. Integers cannot be
strings, booleans, or floats. NaN, Infinity and duplicate JSON keys are rejected.

Reference metadata does not trigger a lookup or convert coordinates. FASTA import
uppercases A/C/G/T; JSON sequence fields are strictly uppercase. Ambiguous N/IUPAC
bases are unsupported in this model rather than guessed or silently discarded.

## Coordinates without ambiguity

For reference `ACGTACGT`, `[2,5)` is `GTA`.

```json
{"start": 2, "end": 5, "sequence": ""}
```

This deletes those three bases. A substitution at local index 2 is `[2,3)` with a
one-base alternate. An insertion immediately before local base 2 is `[2,2)` with
a nonempty alternate. All replacements refer to the unchanged input reference.

A 1-based genomic/VCF position is not a local position. For a supplied window
starting at zero-based genomic coordinate `g`, a simple 1-based genomic position
`p` maps to local base `p - 1 - g`. That formula does not normalize VCF alleles:
VCF anchoring, strand orientation, left normalization, and complex alleles require
a dedicated importer, which is not shipped here.

## Primer orientation and reads

`left_primer` and `right_primer` are intervals on the supplied forward-oriented
reference, with a nonempty gap between them. Both actual primer oligos are
specified 5′→3′. `left_oligo` equals the left interval sequence;
`right_oligo` equals the **reverse complement** of the right interval sequence.
Optional oligos are checked exactly against their intervals.

Product size includes both primer intervals. The observed insert excludes both.
`min_product_bp` and `max_product_bp` are inclusive, declared inclusion bounds;
they are not automatically inferred from a polymerase or instrument.

`paired_end` requires `read_bases`: usable **post-primer-trim insert bases per
end**, not nominal sequencer cycle count. `full_insert` forbids `read_bases`.

## Output families

Schema 1.1 `editwitness.analysis` contains full allele edit definitions, a distinct-alternative count, generation provenance, the reference digest, assays, per-allele
observations (including every retained exact product), per-hypothesis signatures, canonical genotype representatives, comparisons, witnesses, candidate-panel
plan, assumptions, notices, model version and checksums.

`editwitness.deletion_scan` contains the explicit grid, denominator, counts, up to
20 blind examples, limitations and checksums.

`editwitness.summary` is a compact projection of an analysis. It does not contain
the full evidence and cannot be verified as a full result. Its `analysis_sha256`
identifies the complete analysis from which it was derived.

`editwitness.witness` returns one named equivalent alternative. With
`--include-sequences`, it includes the relevant allele observations across
existing and candidate assays. It is a focused evidence view, not a separately
sealed analysis.

`manifest_validation` and `integrity` outputs establish structural consistency or
content integrity only. Their names must not be repurposed as scientific
validation status.

## Files and limits

Input JSON and FASTA are bounded at 8 MiB; loaded full result JSON at 64 MiB.
Analysis also has sequence/work budgets described in the architecture document.
A valid but too-large design must be deliberately split, never silently sampled.

Output defaults to stdout. Files are created atomically and are not overwritten
unless `--force` is supplied. Input replacement is refused even with `--force`.
Parent directories must already exist. Temporary file permissions default to
owner-only; operating-system behavior and filesystem support still apply.

Each output file is atomic independently. An analysis JSON plus HTML report is
not a multi-file transaction: a late filesystem error can leave one complete
file and no other file. Use result checksums and exit status; do not treat a
partial output set as completed work.

## Product collections and aliases in 1.1

For exact-local observations, `products` is authoritative. Each product has
`start`, `end`, `orientation`, `product_length`, ordered `reads`, and `signal_id`.
Coordinates are on the edited allele. A multi-product observation intentionally
has no singular `signal_id` or `product_length`; null is not absence of signal.
The engine unions signal identifiers across **all** retained products.

`representative_hypothesis` identifies the canonical representative of a local
genotype. `same_local_genotype_as_expected` distinguishes a renamed/reordered
expected state from a distinct but observationally confounded alternative.
`distinct_alternatives` counts canonical alternatives, not input labels. The
`witness` command accepts an alias but identifies the requested and representative
IDs explicitly.

`editwitness.model_comparison` is a summary of analyses under two models, with
analysis digests, shared/changed witnesses and panel changes. It is not itself a
sealed full analysis. `expand-deletions` emits an ordinary full manifest with
explicit alleles/hypotheses and declared provenance, not an opaque generator job.

See [migration notes](migration-0.2.md) for result compatibility. Schema 1.0 or
unknown-version results are rejected by the current verifier rather than migrated
under an old digest. Use the package that produced an older result to verify it.
