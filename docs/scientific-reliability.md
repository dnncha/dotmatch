# Scientific Reliability Contract

DotMatch makes a narrow, testable claim: known short DNA sequences are assigned
deterministically from explicitly extracted read windows, and counts never
silently absorb ambiguous, unmatched, invalid, or lossy input. Reliability work
protects that contract; it does not by itself prove biological correctness.

## What is guaranteed

- **Closed assignment states.** Every read ends as `unique`, `ambiguous`, `none`,
  or `invalid`. `MatchStatus` is an `IntEnum` and stream summaries are typed
  `TypedDict`s, so the states and counters are checked at the API boundary.
- **Oracle equivalence.** The native scan, indexed, Hamming, exact, and policy
  paths are compared with a scalar reference implementation over generated
  targets and reads.
- **Count conservation.** `total_reads == unique + ambiguous + none + invalid`
  for every streamed summary; rates are derived from those counters.
- **Unique-only counting.** Ambiguous, unmatched, and invalid observations never
  enter count matrices. This is enforced by the feature-matrix and AnnData
  adapters and pinned by example tests.
- **Lossless counts.** Raw counts are parsed as exact integers, never floats;
  missing values, duplicates, negative values, fractional values, overflow, and
  unselected-but-invalid sample cells fail explicitly. There is no zero-fill
  fallback.
- **Input integrity.** FASTQ structure, printable ASCII, sequence/quality length
  agreement, repeated-identifier agreement, gzip/plain equivalence, and target
  ID uniqueness are checked before assignment.
- **Determinism.** Seeded panel design and feature-matrix artifacts are
  reproducible for the same inputs; assignment TSVs are byte-identical.
- **Boundary fuzzing.** Generated valid and corrupted FASTQ, count TSV, and
  AssaySpec TOML inputs must either parse into typed records or raise a
  documented input error, never an unexpected exception.
- **Mutation sensitivity.** Pinned cases prove the oracle and count tests
  distinguish radius-vs-best, exact-vs-k1, unique-only-vs-count-all, and
  lossless-vs-float count semantics.
- **Typed public surface.** The entire `dotmatch` and `quickdna` package is
  checked by mypy strict and pyright at the Python 3.9 language level.

## Enforced by

```bash
make scientific-tests   # oracle, conservation, lossless counts, determinism
make quality            # typecheck + lint + scientific tests + full tests + 3.9
make python-test        # full Python suite
```

Key files:

- `python/tests/test_scientific_invariants.py`: scalar oracle equivalence,
  stream conservation, count round-trips above float precision, gzip/plain
  equivalence, deterministic panel design.
- `python/tests/test_scientific_fuzz.py`: FASTQ/count/AssaySpec fuzzing,
  deterministic assignment/panel/feature-matrix artifacts, and
  mutation-sensitivity cases.
- `python/tests/test_native_oracle_audit.py`: independent distance and policy
  oracle at length boundaries.
- `python/tests/test_input_validation_boundaries.py`: FASTQ and count parsing
  boundaries.
- `python/tests/test_feature_matrix.py`, `test_sparse_integrity_audit.py`,
  `test_sensitivity*.py`: count-matrix and artifact integrity.
- `docs/scientific-claims.md`: the evidence boundary for public performance and
  comparison claims.

## What these tests do not prove

They do not prove biological accuracy, guide efficiency, cell calling, UMI
molecule counts, or equivalence to a comparator. Those need the workflow-specific
evidence gates and recorded datasets described in `docs/scientific-claims.md`.
Passing the scientific invariants means the implementation preserves its stated
semantics, not that the assay or biology is correct.

When behaviour changes, update the scalar oracle and add a focused invariant or
regression to `test_scientific_invariants.py` before changing the implementation.
A change that cannot be expressed as a check is not ready to ship.
