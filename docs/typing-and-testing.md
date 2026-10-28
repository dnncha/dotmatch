# Typed Python Surface and Test Flow

DotMatch ships `py.typed`, so editors and type checkers see the public API.
The Python runtime floor remains 3.9; the type checker targets 3.10 because
current mypy no longer supports 3.9 as a check target.

## The type vocabulary

The public vocabulary lives in `python/dotmatch/_types.py`:

- `MatchStatus` is an `IntEnum` with four closed values: `INVALID`, `NONE`,
  `UNIQUE`, `AMBIGUOUS`. The historical `MATCH_*` names remain aliases, and
  every value is still an `int`.
- `Metric`, `AssignmentPolicy`, `PanelMetric`, `DistanceMetric`, and the assay
  enums are `Literal` aliases, so a typo fails static analysis instead of
  reaching a runtime string comparison.
- `AssignmentCounts` and `AssignmentSummary` are `TypedDict`s with integer
  counters and float rates, so a misspelled key or a rate used as an integer
  fails static analysis.
- `TargetRow` is a `NamedTuple` of `target_id` and `sequence`.
- `PanelCertificate` is the typed shape of `panel check` output.
- `require_int`, `require_float`, and `require_str` are the parse boundary for
  TOML, CLI, JSON, and other unstructured input.

Internal code trusts these types. New external input gets parsed once at the
boundary rather than cast at each use site.

## Static analysis

```bash
make typecheck   # mypy --strict plus pyright
make lint        # ruff: pyflakes, import order, error-level pycodestyle
make quality     # typecheck + lint + Python tests
```

`python/tests/typecheck_api.py` is analyzed by both checkers but never
collected by pytest. It uses `assert_type` for positive contracts and
`# type: ignore[...]` for negative contracts; `warn_unused_ignores` proves the
checker really rejects the bad call.

## Example and snapshot tests

```bash
make examples-test
UPDATE_EXPECT=1 make examples-test   # re-record intentional changes
```

`python/tests/test_examples_typed.py` calls the public API the way a user does
and pins deterministic output with the Rust `expect-test` style helper in
`python/tests/expect.py`. Snapshots live in `python/tests/expected/`.

Public docstrings carry executable `>>>` examples. `test_doctests.py` runs
them for `dotmatch._types`, `core`, and `panel`.

`python/tests/test_invariants.py` adds Hypothesis properties for distance
symmetry and thresholds, equal-length Hamming/Levenshtein bounds, exact vs
Hamming k=0 agreement, radius-vs-best policy direction, and stream batch-size
invariance.

`examples/python/typed_quickstart.py` is a standalone, typed example and is
executed by the example tests.


## Current coverage

The strict gate covers the whole `dotmatch` and `quickdna` packages, plus the
type tests, snapshots, doctests, typed example, and binding benchmark. `pyright`
checks at Python 3.9; mypy checks at 3.10 because current mypy no longer
supports a 3.9 check target. Runtime imports and the example suites run on 3.9.
