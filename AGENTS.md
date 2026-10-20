# Working on DotMatch

DotMatch is a scientific C engine with Python APIs and a Next.js website.
Read [docs/architecture.md](docs/architecture.md) for ownership and checks.

- Matching algorithms and memory ownership live in `src/qdalign.c`; public C ABI
  lives in `include/qdalign.h`. Keep ctypes declarations in
  `python/dotmatch/_bindings.py` synchronized with that header.
- Reuse `fastq_io.py`, `target_io.py`, `count_io.py`, and `_validation.py` for
  parsing and validation. Preserve IDs, exact integer counts, literal sequence
  symbols, ambiguity statuses, and streaming behavior.
- Implement Python behavior in `dotmatch`. `quickdna` and `assaycode` are
  compatibility surfaces; preserve their public imports and entry points.
- Library code must not import CLI orchestration or repository scripts. Keep
  optional scientific dependencies lazy. `make architecture-check` enforces
  these boundaries without importing the native engine.
- Use Python 3.9-compatible runtime syntax/behavior. Tooling runs on 3.11+;
  the compatibility matrix verifies 3.9 through 3.14.
- Activate a venv, install `requirements-quality.txt` plus `pytest build`, and
  run `make check` before submitting Python/native changes. Run
  `make python-package-test` for ABI, packaging, or import changes. Run
  `make native-analyze`, `make asan`, and `make coverage` for native changes.
- For website changes run `npm ci`, `npm run lint`, `npm run typecheck`,
  `npm run check:site`, and `npm run build`.
- Fix failing checks at their source. Keep any unavoidable tool suppression
  local, explain it, and add a regression test for the behavior it protects.
- Evidence, public-claim, and release requirements remain in
  [CONTRIBUTING.md](CONTRIBUTING.md). Run only the gates relevant to the change;
  routine code work does not require downloading public sequencing datasets.

Keep changes focused. Avoid broad rewrites, regenerated benchmark artifacts,
or compatibility changes merely to satisfy a style preference.
