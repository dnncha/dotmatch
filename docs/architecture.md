# Architecture

## One scientific engine, several interfaces

```text
Local JSON / local FASTA
          │
 strict immutable contracts ───────────► JSON Schema
          │
 sequence transforms + original-site observation model
          │
 per-allele signals → per-hypothesis signal sets
          │
 explicit equivalence witnesses → candidate coverage → panel planner
          │
 sealed, versioned Analysis
          ├── Python objects
          ├── full JSON / compact agent summary
          └── script-free HTML evidence report

Separate path: bounded deletion grid → streaming interval scan → sealed ScanResult
```

The report never decides a scientific result. The CLI never owns a second
implementation of the observation model. The scan and the equivalence analysis
are separate because they answer different questions.

## Module responsibilities

| Module | Responsibility |
|---|---|
| `models.py` | Strict types, immutability, local coordinates, cross-field invariants, budgets. |
| `sequence.py` | Small sequence and interval primitives. |
| `observations.py` | One explicit original-site observation function. |
| `engine.py` | Signal caching, hypothesis comparisons, witnesses, qualifications. |
| `planner.py` | Exact bounded subset search; deterministic greedy fallback. |
| `scan.py` | Streaming single-deletion geometry, not genotype inference. |
| `io.py` | Bounded parsing, duplicate-key rejection, canonical hashes, atomic writes. |
| `fasta.py` | Single-record local import and exact primer-site initialization. |
| `cli.py` | Machine-oriented commands, exit codes, compact responses. |
| `report.py` | Presentation of computed evidence, escaped and offline. |

## Decisions made deliberately

**Python first, not Rust-first packaging.** The priority is an inspectable model,
a portable wheel, and straightforward adoption by existing scientific Python
workflows. The runtime has one direct dependency, Pydantic. The low-level
sequence helpers have no third-party dependencies. A native extension should be
added only after a reproducible profile identifies a real bottleneck; it must
pass the independent oracle tests unchanged.

**No LLM in the calculation.** Agents discover the CLI/schema and invoke the same
pure functions as a human or workflow runner. No model provider, key, hosted
backend, MCP server, or subscription is required. An eventual MCP adapter should
be a thin, separately packaged wrapper, not a second scientific implementation.

**Explicit finite hypotheses.** We do not hide an undocumented biological prior
inside a score. Named alleles and alternatives let users inspect exactly what
was tested and let agents preserve the distinction between input assumptions
and computed results.

**Qualitative sequence presence.** This initial model avoids interpreting PCR
read fractions as genomic dosage. A quantitative model would require evidence
for amplification, sampling, and error assumptions. It cannot be introduced by
merely adding a floating-point confidence field.

**Version the science.** Package, schema, and response-model versions are distinct.
The output includes all three. Altering eligibility or observation semantics
requires a new model identifier and migration notes. An optimization must not
silently alter the observation function.

## Determinism and provenance

Canonical JSON sorts keys and uses fixed separators and ASCII escaping. The
manifest digest is over validated, normalized values including defaults. Results
include package/model versions and a digest over all result fields except the
digest itself. Timestamps, environment-specific paths, and random IDs are not
included. Reordering input arrays can change the digest: this is not a semantic
normal form over all equivalent designs.

`verify --manifest` recomputes the analysis or scan. A digest is an integrity
check, not an author signature, biological validation, or protection against an
attacker who replaces both content and digest.

## Bounded computation and presentation

Reference: 20,000 bases. At most 128 alleles, 1,000 diploid hypotheses, 16 existing
assays, and 24 candidates. The sum of reconstructed/reference work across
alleles and assays is conservatively limited to 20 million bases, including
large inserted sequences. Deletion scans reject grids over 500,000 endpoint
pairs before allocating per-grid work.

The planner's exact path considers at most 262,144 subsets; larger candidate
sets use a clearly labeled heuristic. The scan keeps aggregate counts and at
most 20 blind examples, not every simulated sequence.

HTML shows at most 50 counterexamples and 2,400 bases per displayed read,
explicitly labeling truncated previews. Full JSON retains the complete bounded
analysis. Report limits cannot change witnesses, panel choices, or JSON evidence.

## Extension rules

A future response model should return explicit observations and qualifications
through a versioned boundary. Do not overload `full_insert` with an assay that
observes something else. New genomic event classes need normalized semantics and
an independent sequence-level oracle. Raw-read adapters must preserve upstream
provenance and distinguish unobserved outcomes from zero-frequency outcomes.

Do not add arbitrary executable plugin loading from manifests. Do not fetch
references implicitly. Batch and workflow integrations should compose the local
CLI and use distinct output paths per manifest.

## Deferred on purpose

No cloud, database, jobs service, authentication system, dashboard framework,
GPU layer, distributed scheduler, or model-training pipeline. None is necessary
to make the first scientific question useful. Add infrastructure only after
independent users encounter a concrete limit.
