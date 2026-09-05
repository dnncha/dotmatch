# Validation: what is established and what is not

## Current status

This alpha has deterministic software tests and deliberately synthetic examples.
It has no prospective laboratory validation, no measured sensitivity or
specificity, no externally adjudicated assay-design benchmark, and no clinical
validation. No published data were relabeled as a successful package benchmark.

For exact executed checks, see [BUILD_STATUS.md](../BUILD_STATUS.md). A configured
CI job, a proposed test, or a packaged schema is not a passing result.

## Software verification strategy

**Independent sequence oracle.** Exhaustive small deletion, insertion, substitution
and replacement cases are reconstructed by an independent labeled-base oracle.
It determines surviving original primer sites using base provenance, not the
production disruption/mapping helpers. Tests compare full-insert and paired-end
observations, boundaries, empty inserts, and product-size eligibility. Randomized
compound examples use a fixed seed.

**Panel optimality.** Small randomized set-cover instances are checked against an
independent combinations-based enumeration. Ties, no candidates, unrecoverable
alternatives, and the deliberately nonoptimal larger-candidate path are tested.

**Scan agreement.** Interval-only deletion counts are compared with the sequence
observation path across small exhaustive grids, with length and size filters.
This tests internal model consistency, not empirical deletion frequency.

**Contracts and execution.** Tests cover strict types, bounds, duplicate fields
and IDs, immutable inputs, both schema and semantic validation, altered result
hashes, replay, malformed input, atomic writes, overwrite refusal, structured
exit codes, packaged examples and escaped reports. The core is tested with
socket creation blocked.

These tests substantially constrain implementation mistakes. They do not make an
idealized PCR model realistic merely by being numerous.

## Independent biological validation plan

### Gate A: adjudicated examples

Obtain public or explicitly permitted examples with a reference assembly, exact
primer sequences and coordinates, true edit structures, read configuration,
quality-control measurements, and an independent assessment of the relevant
alleles. Record accession, license, checksum, provenance and exclusions before
analysis. An article abstract alone is not sufficient benchmark metadata.

Include known hidden-deletion examples **and negative controls** where the
measurement truly separates the competing states. Have an independent scientist
review each model encoding without seeing the package prediction.

### Gate B: model adequacy

Challenge assumptions deliberately: boundary changes, repetitive sites,
replacement-induced rescued sites, primer mismatches, size-dependent bias,
unequal allele amplification, incomplete reads, sampling, mixtures, and
copy-number changes. Record failures instead of shrinking the benchmark to
well-behaved examples. Decide which require a new model and which remain excluded.

### Gate C: prospective decision usefulness

At least two independent facilities should use the package on a second design
without developer-led encoding. Measure whether it changes a validation decision,
how much manual input is needed, whether its warnings are actionable, and whether
an independent measurement supports the specific ambiguity it highlighted.

No arbitrary "accuracy" score should combine incomplete metadata, synthetic
cases and real biological examples. Report denominators, excluded cases,
confidence/uncertainty where justified, and observed model failures separately.

### Gate D: broader release

A methods claim requires an independently reviewed model, audited benchmark
provenance, reproducible scripts, negative controls, and evidence that users can
interpret reports without mistaking them for safety certificates. Publication
and integration work should follow that evidence, not replace it.
