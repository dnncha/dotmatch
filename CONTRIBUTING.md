# Contributing

A useful contribution changes a real scientific or operational decision while
keeping the evidence inspectable. Small, reviewable changes are preferred.

## Local setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python -m pytest -q
python scripts/generate_schemas.py --check
python scripts/check_style.py
python -m mypy src/editwitness
python -m build
```

On Windows activate the environment with `.venv\Scripts\Activate.ps1` instead.
For branch coverage including CLI subprocesses:

```bash
python scripts/check_coverage.py
```

Subprocess coverage is slower than the plain suite. Keep its results separate
from benchmark timing. Developer tools are not runtime dependencies.

## Before changing scientific behavior

Read `AGENTS.md` and the scientific-model document. Provide a minimal synthetic
counterexample, an independently derived expected result, and a regression test.
Explain whether the change fixes implementation of the current model or changes
the model itself. The latter needs a new model version and migration notes.

When updating contracts, regenerate JSON Schema snapshots and examples, then
check both structural and semantic validation. Preserve strict errors rather
than silently coercing questionable inputs.

## Data and disclosure

Do not attach patient-derived or confidential sequences to public issues.
Use a synthetic reproducer or an appropriately licensed public dataset with
provenance. Identify synthetic data explicitly. Cite primary work accurately and
credit adjacent tools instead of claiming their established ideas as new.

See `roadmap.json` for bounded tasks with acceptance criteria. Update task status
only with executable evidence or an actual external review outcome.
