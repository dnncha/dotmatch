# Architecture

## Small deterministic core, multiple interfaces

```text
Local JSON / FASTA
  -> strict immutable manifest and cross-field validation
  -> optional bounded deletion expansion + provenance
  -> cached reconstruction of each final allele
  -> selected versioned observation function
  -> per-hypothesis signal sets and sequence-state identity
  -> ambiguity witnesses + candidate coverage
  -> bounded exact / greedy panel selection
  -> typed result + checksum
  -> full JSON / compact JSON / static HTML / Python API
```

`models.py` owns immutable Pydantic contracts; public `analyze` and generation
entry points revalidate even constructed/copied models. `sequence.py` owns edit
application. `exact.py` rematches final sequences and enumerates products in
both orientations. `observations.py` retains the separately named legacy
function. `engine.py` composes observations without filesystem or network access.
`planner.py` treats panel selection as a finite set-cover objective; `design.py`
generates finite deletion challenges. `scan.py` is a distinct geometry scanner.

`io.py` handles bounded file input, content hashing and atomic non-destructive
writes. `cli.py` owns process exit semantics. `report.py` is a script-free escaped
HTML view of the same result, not a second scientific calculation.

## Performance without hidden approximation

Each final sequence is reconstructed once per analysis. Sorted primer hits and
binary search limit product pairing to nonoverlapping, size-eligible ranges;
paired-end mode slices only sequenced ends, not long unobserved gaps. Shared
budgets limit aggregate product evidence, rather than multiplying a permissive
per-call budget by the number of alleles. Panel dominance removes redundant
choices before combinatorial optimization.

Limits are part of the interface. The engine never selects an arbitrary first
primer match or quietly truncates a counterexample set. Do not substitute a
native implementation until an independently checked workload demonstrates a
benefit. Preserve an independent reference oracle.

## Trust boundaries

No network, LLM, telemetry, executable manifest content, shell expansion, or
remote sequence retrieval occurs during analysis. User sequences remain local
unless the user explicitly shares outputs. Reports and witness JSON can contain
sensitive sequence information and must be treated accordingly.

Schemas, package versions and biological model identifiers are separate. Adding
a probabilistic or empirical response function requires a new model identity,
calibration evidence and migration documentation. Existing callers are future
read-only adapters, not a reason to silently reinterpret their outputs.

Publishing is a separate authenticated operator script. It checks source and
artifact inventories, exact commit CI, visibility and identity. It never uses
the analysis package as a credential store. See [release process](releasing.md).
