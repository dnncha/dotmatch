# Architecture and development checks

DotMatch keeps its matching semantics in the native engine and exposes the same
behavior through Python, command-line tools, and compatibility packages.

| Responsibility | Owner | Constraint |
| --- | --- | --- |
| Matching, indexes, allocation | `src/qdalign.c`, `include/qdalign.h` | Public C ABI; deterministic statuses and distances |
| Native CLI and FASTQ workflows | `src/qda.c`, `src/target_table.h` | Preserve native/Python input parity |
| ctypes ABI and library discovery | `python/dotmatch/_bindings.py` | Sole shared-library loader; trusted explicit/package/source paths |
| Matching API and index lifetime | `python/dotmatch/core.py` | Uses bindings; no CLI or repository-tool dependencies |
| Strict shared input parsing | `fastq_io.py`, `target_io.py`, `count_io.py`, `_validation.py` | Standard library plus the declared validator dependency |
| Optional scientific conversions | `dataframes.py`, `tl.py`, `multiqc.py` | Load scientific stacks on use; MultiQC owns its plugin imports |
| Python CLI routing | `entrypoint.py`, `cli.py`, `first_run.py` | Compose library operations |
| Compatibility | `python/quickdna`, `python/assaycode` | Depend on DotMatch; never the reverse |
| Evidence and repository automation | `scripts/` | May use the runtime; runtime must not import scripts |
| Website | `app/`, `lib/` | Next.js/TypeScript; separate lint, type, test, and build gates |

`scripts/check_architecture.py` statically checks internal import direction,
low-level dependencies (including function-body imports), eager optional
dependencies, eager import cycles, and ownership of ctypes library loading.
The tests deliberately introduce forbidden imports and cycles. This is a check
of ordinary Python imports and loader calls, not a sandbox for dynamic code.

## Local workflow

Use Python 3.11 or newer for development tooling. The runtime still supports
3.9, and CI exercises Python 3.9–3.14.

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-quality.txt pytest build
make check
```

`make quality` runs correctness lint across the Python runtime, tests, and
scripts; formatting and strict typing on the boundary modules and architecture
checker; architecture checks; and native compiler warnings as errors. It does
not need the optional scientific stack, downloaded datasets, or a native build.
`make check` follows it with the native, CLI, and Python test suites.

Strict typing currently covers the native loader, validators, shared parsers,
optional-module loader, native CLI bridge, and architecture checker, as listed
in `tool.mypy.files`. The larger analytical modules are not yet fully typed.
Add new boundary modules to that list and `PYTHON_FORMAT_PATHS`; expand coverage
when extracting responsibilities from the larger modules. Ruff correctness
rules cover existing and new files in the checked directories automatically.

`make python-format` formats the checked boundary files. Other legacy modules
are intentionally not reformatted wholesale. Optional local hooks use the same
installed versions and Make targets as CI:

```bash
pre-commit install
pre-commit run --all-files
```

For native changes, `make native-analyze` uses GCC's analyzer in a dedicated
output directory. `make asan` and `make coverage` use the existing sanitizer and
coverage checks. These targets rebuild native outputs; run them sequentially
with package tests. On macOS the compiler gate also checks the Metal source.

For web changes:

```bash
npm ci
npm run lint
npm run typecheck
npm run check:site
npm run build
```

ESLint checks TypeScript and React hooks with zero warnings allowed. The
`@typescript/native` alias retains the TypeScript 7 checker; `typescript` 6
provides the compiler API required by the current ESLint parser and Next.js.
The typecheck command selects the native checker explicitly so npm executable
link order cannot silently choose a different compiler.

## CI behavior

The fast quality job must pass before native/platform, Python compatibility,
website, and reviewer-reproducibility jobs run. Native static analysis and
sanitizer/coverage checks run independently on Linux. CI uses read-only
repository permissions, cancels superseded PR runs, limits job duration, and
uploads useful test/analysis evidence on failure. Existing scientific evidence,
package integrity, CodeQL, and release gates remain authoritative.

The stable `CI required` job fails if any required CI job fails, is cancelled,
or is skipped unexpectedly. Repository administrators can require this one
check in branch protection. Adding a workflow does not itself change repository
branch-protection settings.
