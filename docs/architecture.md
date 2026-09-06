# Architecture

## One deterministic engine, multiple interfaces

```text
Local JSON / local FASTA
           |
Revalidated strict contracts ----------------> JSON Schema
           |
Optional bounded deletion expansion ----------> explicit, inspectable manifest
           |
Reconstruct local allele sequences once
           |
Explicit versioned observation model
  original sites v1  |  exact rematched local products v2
           |
Allele signal sets -> diploid signal unions
           |
Canonical local genotype groups -> concrete counterexamples
           |
Candidate coverage -> exact or explicitly heuristic panel selection
           |
Versioned, hashed Analysis
   Python objects | full JSON | compact agent summary | script-free HTML

Separate: endpoint grid -> original-site geometry scan -> hashed ScanResult
```

The report never decides a scientific conclusion. The CLI owns no duplicate
model. A model comparison runs the same engine twice with explicit model choices.
The geometry scan is separate because structural eligibility and observational
equivalence answer different questions.

## Modules

| Module | Responsibility |
|---|---|
| `models.py` | Strict immutable input contracts, cross-field invariants and conservative budgets. Public entry points revalidate unchecked Pydantic copies. |
| `sequence.py` | Reference-coordinate edits and sequence/interval primitives. |
| `observations.py` | Historical original-site response function. |
| `exact.py` | Exact edited-sequence matching, both heteroprimer orientations, bounded multi-product enumeration and signal projection. |
| `generate.py` | Explicit reference-deletion alternatives, canonical local-genotype deduplication, provenance, no silent subsampling. |
| `engine.py` | Reconstruction cache, hypotheses, alias grouping, witnesses and qualified conclusions. |
| `compare.py` | Paired response-model analysis and a compact sensitivity comparison. |
| `planner.py` | Exact bounded minimum-cost search; deterministic rational weighted-cover fallback. |
| `scan.py` | Streaming original-site deletion geometry, not sequence/genotype inference. |
| `io.py` | Bounded UTF-8 JSON, duplicate-key rejection, canonical digests, protected atomic writes. |
| `fasta.py` | Single-record local reference import and exact primer-site initialization. |
| `cli.py` | Stable machine-oriented commands, structured errors, compact responses. |
| `report.py` | Escaped, offline evidence presentation with explicit preview limits. |

## Deliberate choices

**Python, with one direct runtime dependency.** Pydantic provides input validation
and schemas. Computation does not require an LLM, service, API key or network.
The sequence/model code is inspectable and portable. Add a native extension only
when measured workloads justify it, preserving oracle agreement and wheels.

**Two scientific models, not a silent algorithm swap.** Existing omitted-model
manifests retain v1. New CLI examples explicitly select v2. Package, schema and
model versions are separate. Schema 1.1 results expose products, canonical
hypothesis representatives and generation provenance; old evidence requires its
producing package for integrity/replay verification.

**No probabilistic score without a probabilistic model.** The signal is qualitative
sequence presence, not genomic dosage or read fraction. Adding a confidence
number would require a specified statistical model and independent calibration.

**Agent use is ordinary deterministic software use.** A tool discovery command,
schemas, typed Python functions and consistent JSON errors are the primary
interface. No arbitrary executable plugins or remote reference lookup can be
requested by a manifest. An optional future MCP wrapper must stay thin and keep
scientific semantics in this package.

## Determinism and provenance

Canonical JSON sorts keys and fixes separators and ASCII escaping. The input hash
covers normalized validated values, including defaults. The result hash excludes
only itself. Timestamps, random IDs and environment-specific paths do not affect
analysis output. Input array ordering can change the digest: this is not a
semantic hash over all equivalent designs.

Genotype comparison is different: unordered pairs of reconstructed sequences are
canonicalized with multiplicity. Aliases retain provenance but cannot inflate
counterexamples or planning gains. `verify --manifest` replays the calculation;
a checksum alone is neither an authenticated signature nor proof of correctness.

## Bounded work

The manifest caps 20 kb of reference, 128 alleles, 1,000 hypotheses, 16 baseline
assays and 24 candidates. Reference/reconstructed work across assays is limited
to 20 million bases. Exact products are generated using sorted hit lists and
binary-search bounds, with hard limits of 512 matches per search, 128 products
per allele/assay and 20 million total product bases. No partial result survives a
resource-limit error.

Deletion generation rejects over 5,000 grid cells or an expanded manifest that
would exceed any contract bound. The separate scanner accepts 500,000 cells and
stores at most 20 examples. Exact panel planning uses at most 262,144 subsets;
the larger-candidate path declares itself heuristic.

HTML previews at most 100 hypotheses, 50 witnesses, 8 products per observation
and 2,400 bases per displayed read. Limits are labeled; full bounded JSON remains
authoritative. UI limits never change analysis or selection.

## Extension rules

New assays need explicit response-model semantics and independent tests. New
biological claims additionally need external evidence. Raw-data adapters must
preserve upstream files, coordinate definitions and source hashes. Unknown
readouts must be rejected, not interpreted as full insert. Batch execution should
use separate outputs and retain model versions and limits per sample.

No cloud platform, database, GPU layer, job service or model training pipeline is
required for this question. Introduce infrastructure only in response to a
measured need from actual users.
