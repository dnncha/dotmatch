# Data contract (schema 1.1)

Use `editwitness schema manifest` or `schema analysis` for the complete generated
JSON Schema. Schema snapshots ship inside the wheel and are checked for drift.
Cross-field/biological invariants also require `validate`; JSON Schema alone is
not a replacement for the runtime validator.

## Input

All edit and primer coordinates are **zero-based, half-open local reference
coordinates**, not genomic positions unless the supplied local sequence happens
to begin at genomic zero. Oligos are 5′-to-3′ strings; the right primer is the
reverse complement of its annotated plus-strand reference interval. Full-insert
and paired-end readouts require different explicit configuration.

Every hypothesis names exactly two declared alleles. IDs are unique within their
namespace; referential integrity and nonoverlapping edit invariants are checked.
Extra fields and invalid primitive types are rejected. Public APIs revalidate
Pydantic instances rather than trusting `model_copy` or `model_construct`.

Schema 1.0 and 1.1 inputs are supported. An omitted observation model means legacy
v1 for backward compatibility. New inputs should set schema 1.1 and explicitly
select `exact-local-sequence-presence-v2` or `original-sites-presence-v1`.

Optional `generation` describes the declared source of expanded alternatives;
it is provenance metadata, not an authenticated proof of how a file was created.

## Output

`model_version` identifies the actual observation function. `schema_version`
identifies serialization. `package_version` identifies executable behavior.
`allele_evidence` retains edit definitions, final sequence length and SHA-256 for
all alleles, including those with no sequence signal.

In v2 an allele/assay may produce multiple products and multiple distinct signals:

- `products` preserves final-allele plus-strand site coordinates, orientation,
  product length and measured reads for every eligible heteroprimer product.
- `signal_ids` is the authoritative set of distinct measured sequence signals.
- Legacy `signal_id` and `reads` describe a singleton distinct signal only;
  **null `signal_id` does not imply absence when `signal_ids` has multiple items**.

A genuinely absent signal differs from an empty insert sequence. Product lengths,
allele IDs, hashes of whole genomes, and binding-site coordinates must never be
used as extra observed data in paired-end equivalence comparisons.

`same_local_genomic_state_as_expected` marks representation/name aliases, excluded
from `witnesses`. `plan.unresolved_hypotheses` is never silently discarded.
`dominated_candidates` explains safe optimizer preprocessing.

## CLI process contract

| Exit | Meaning |
|---|---|
| 0 | Requested computation completed; no statement of biological safety. |
| 2 | Invalid input, unsupported configuration, or resource limit exceeded. |
| 3 | I/O failure. |
| 4 | Analysis completed with ambiguity and `--fail-on-ambiguity` was requested. |
| 5 | Checksum or same-version replay mismatch. |

JSON goes to stdout unless `--output` is supplied. Errors are JSON on stderr.
`--compact` emits a summary, not a complete replayable analysis. Output paths are
preflighted, inputs are not overwritten, and existing files need explicit
`--force`. Hashes detect content changes, not maliciously re-signed data.

## Replay

`verify result.json` checks the checksum of the actual archived JSON payload,
including old 0.1.0a1 artifacts. `verify result.json --manifest input.json` also
reruns the analysis and requires the exact originating package version. This
avoids claiming that a new schema's additional fields reproduce an old byte
representation. [Migration notes](migration-0.2.md).

## Additional bounded-work checks in 0.2.0a2

The observation engine permits at most 100,000 total hypothesis-to-signal
references in one analysis. Deletion generation permits at most 200 million
reconstructed bases across its valid grid cells, including duplicate states.
Both limits are disclosed by `capabilities`; exceeding either returns a
structured input error without partial evidence. The scientific observation
models and schema identifiers are unchanged by these execution guards.

`self-test` emits `editwitness.software_self_test` schema 1.0. Its `passed` field
and process status (0 or 6) concern the installed software only.
