# Build and validation status — 0.1.0a1

## Executed locally

- Python 3.13.5 on Linux x86_64; Pydantic 2.13.4.
- Built a pure-Python wheel and a source distribution. Installed the wheel into
  a separate directory with preinstalled runtime dependencies, then ran packaged
  examples, analysis and replay from outside the source tree. Fresh network
  dependency resolution was not tested.
- Full suite: **572 tests passed**. The suite includes parameterized/exhaustive
  synthetic cases, not 572 independent biological experiments.
- Covered full suite: **572 tests passed**, with 732/733 executable statements and 218/222 branches covered
  (99.48% combined statement/branch coverage). Coverage includes CLI
  subprocesses; it measures software exercise, not scientific validation.
- Structural schema snapshots match the runtime schemas; semantic validation is
  tested separately.
- Dependency-free source-hygiene check passed. This is not a substitute for a
  static type checker.
- Deterministic synthetic full-insert and paired-end examples reproduced, with
  checksum verification and replay.
- Chromium render inspected at 1440 px desktop and 390 px mobile widths, with no
  document-level horizontal overflow and no script elements.
- Synthetic warm benchmark: a 450,000-endpoint grid yielded 325,250 valid
  deletions; five-run median 0.135 seconds in this environment. No empirical
  assay-performance claim or competitor comparison is implied.

Machine-readable verification details are in `docs/verification.json`.
No unexecuted check is marked as passed.

## Not established

No biological validation, measured sensitivity/specificity, external scientific
review, prospective facility adoption, PyPI/Bioconda publication, hosted docs or
DOI. Strict mypy/Ruff execution was unavailable in the local build environment;
strict mypy and a Linux/macOS/Windows test matrix are configured in GitHub Actions
but must be checked against actual run results. No cross-platform pass is implied
by configuration alone.

See `docs/validation.md` and `roadmap.json` for acceptance criteria and remaining
work. Run `python scripts/check_coverage.py` to reproduce branch-aware coverage.
