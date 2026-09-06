# EditWitness 0.2.0a1 — build and publication status

**6 September 2026.** Locally tested research alpha. **Not newly published to
GitHub or PyPI.** No remote repository was modified during this audit.

## Executed checks

| Check | Actual result |
|---|---|
| Original 0.1.0a1 baseline | 572 tests passed locally before changes. |
| Final audited suite | **617 passed, 0 failed**, including CLI subprocesses. |
| Runtime statement coverage | **1,005 / 1,008 — 99.70%**. |
| Runtime branch coverage | **310 / 316 — 98.10%**. |
| Combined coverage | 99.32%; configured 95% gate passed. |
| Wheel and source-distribution builds | Built through installed setuptools PEP 517 backend. |
| Extracted source-distribution tests | **617 passed** outside the development checkout. |
| Extracted source inventory | Verified; legacy JSON fixture and workflow/release files included. |
| Installed-wheel smoke | Passed outside the checkout, using already installed dependencies. |
| Wheel functional paths | Full insert, paired ends, HTML, replay, model comparison and deletion generation. |
| Browser checks | 1,440px desktop and 390px phone; no page-level overflow, scripts, external requests or page errors. |
| Schema snapshots | Match the runtime contracts. |
| Source hygiene | 37 Python files pass Python 3.11 syntax and text hygiene. |
| Relative documentation links | No broken targets found. |
| Publication dry run | Passed source inventory, version and plan checks. |

Environment: Python 3.13.5, Pydantic 2.13.4,
pytest 9.0.2, coverage 7.13.3;
Linux-6.18.35-x86_64-with-glibc2.41. The final instrumented suite completed in 42.16 seconds.
Timings are not a biological or comparative performance benchmark.

Source/test/script fingerprint: `1adf1c4bf36a16b483358609dd7d7a25803449e01b1b9de39d756d3b551cda3d`.
The machine-readable record and per-file coverage are in `docs/verification.json`.
Documentation records were finalized after testing without changing runtime or
test code. Historical alpha CI evidence is preserved in `docs/history/`, not
claimed for this candidate.

## Unexecuted gates

Strict **mypy and Ruff were not available locally** and dependency downloads were
unavailable. They remain required CI jobs. The new six-combination Linux/macOS/
Windows and Python 3.11/3.13 matrix has **not run remotely**. Local wheel smoke
reused existing dependencies; it is not a clean dependency-resolution test. The
configured CI smoke uses a fresh environment and normal dependency installation.

No independent scientific review or adjudicated biological benchmark has been
performed. Test coverage is software evidence, not PCR sensitivity, clinical
validity, scientific novelty, or proof that real alternatives are exhaustive.

## Publication attempt and remaining action

The standalone `dnncha/editwitness` lookup returned Not Found. Current-session
GitHub connector discovery exposed only reads. Running the supplied public
publisher exited 2 with:

```text
gh is required and unavailable here; no remote writes were attempted
```

No standalone repo, new GitHub prerelease or PyPI project was created. DotMatch
was untouched. From this extracted source on an already authenticated local
GitHub CLI installation:

```bash
python scripts/publish_github.py --public --release --dry-run
python scripts/publish_github.py --public --release
```

The publisher creates independent source history, then requires exact-commit
push CI before attaching the checked CI-built wheel and sdist to an alpha
prerelease. Resolve failed gates; do not disable them or claim publication before
it succeeds. See `docs/releasing.md`, the audit ledger and `roadmap.json`.
