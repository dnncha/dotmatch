# Closed-Loop Assay Reliability Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the first shippable closed-loop assay reliability slice to AssaySpec check, plan, and run.

**Architecture:** Keep AssaySpec as the public contract. Add reliability config parsing, backend eligibility, normalized findings, and reliability artifact writers inside the Python AssaySpec workflow layer, because the existing code already owns plan/run/report generation. The reliability engine consumes existing audit, sample QC, CRISPR QC, autopsy, manifest, and assay evidence metadata rather than reimplementing assignment.

**Tech Stack:** Python 3.9+, TOML via `tomllib`/`tomli`, CSV/JSON/HTML standard library, existing native DotMatch CLI.

---

## File Structure

- Modify `python/dotmatch/assayspec.py`: parse `[reliability]` and `[backend]`, add reliability artifacts, evaluate findings, write JSON/TSV/HTML/summary, include reliability in manifest/report, and enforce production preflight blocking.
- Modify `python/tests/test_assayspec.py`: add tests for check artifacts, plan output, run artifacts, exploratory unsafe behavior, production unsafe blocking, backend validation, and report escaping.
- Modify `docs/assayspec.md`: document reliability configuration and generated artifacts.
- Create no new production module for the MVP. A separate module can follow once behavior stabilizes; keeping the first slice in `assayspec.py` avoids circular imports with current private helpers.

## Task 1: Reliability Config And Artifact Surface

**Files:**
- Modify: `python/tests/test_assayspec.py`
- Modify: `python/dotmatch/assayspec.py`

- [ ] **Step 1: Write failing tests for reliability config, plan artifacts, and plan output**

Add tests that expect `compile_assay_plan()` to include reliability artifacts and `dotmatch assay plan` to mention the reliability report without creating output directories.

- [ ] **Step 2: Run the tests and verify RED**

Run:

```bash
DOTMATCH_LIB=$PWD/libdotmatch.dylib PYTHONPATH=$PWD/python python3 -m pytest \
  python/tests/test_assayspec.py::test_load_count_spec_and_compile_deterministic_plan \
  python/tests/test_assayspec.py::test_assay_plan_prints_native_commands_without_creating_outputs -q
```

Expected: fails because reliability artifacts and plan output do not exist.

- [ ] **Step 3: Implement config defaults and artifact registration**

Add constants and helpers in `assayspec.py`:

- profiles: `production`, `exploratory`;
- backend modes: `auto`, `cpu`, `gpu-metal-experimental`;
- default thresholds matching `AUTOPSY_THRESHOLDS`;
- `AssaySpec.reliability`;
- `AssaySpec.backend`;
- validation for `[reliability]` and `[backend]`;
- reliability artifact paths in `compile_assay_plan()`;
- a comment block in `format_plan()`.

- [ ] **Step 4: Run the same tests and verify GREEN**

Run the command from Step 2. Expected: pass.

## Task 2: Preflight Reliability Artifacts For `assay check`

**Files:**
- Modify: `python/tests/test_assayspec.py`
- Modify: `python/dotmatch/assayspec.py`

- [ ] **Step 1: Write failing test for `assay check` artifacts**

Add a test that runs `dotmatch assay check`, then asserts:

- `reliability_summary.json` exists;
- `reliability_findings.tsv` exists;
- `reliability_report.html` exists;
- `reliability_manifest.summary.tsv` exists;
- read-dependent findings are recorded as unavailable;
- the report contains escaped, evidence-bounded content.

- [ ] **Step 2: Run the test and verify RED**

Run:

```bash
DOTMATCH_LIB=$PWD/libdotmatch.dylib PYTHONPATH=$PWD/python python3 -m pytest \
  python/tests/test_assayspec.py::test_assay_check_writes_preflight_reliability_artifacts -q
```

Expected: fails because `assay check` currently only prints `ok`.

- [ ] **Step 3: Implement preflight reliability writer**

Add helpers:

- `_write_preflight_reliability(plan)`;
- `_build_reliability_summary(plan, manifest=None, stage="preflight")`;
- `_write_reliability_artifacts(plan, summary)`;
- `_write_reliability_findings(path, findings)`;
- `_write_reliability_manifest_summary(path, summary)`;
- `_write_reliability_report(path, summary)`.

`assay check` should compile the plan, write the preflight artifacts, print `ok`,
and return zero.

- [ ] **Step 4: Run the test and verify GREEN**

Run the command from Step 2. Expected: pass.

## Task 3: Runtime Reliability Aggregation

**Files:**
- Modify: `python/tests/test_assayspec.py`
- Modify: `python/dotmatch/assayspec.py`

- [ ] **Step 1: Write failing test for run reliability artifacts**

Extend the existing count-run test to assert:

- manifest artifact map includes all reliability artifacts;
- `reliability_summary.json` records `overall_status = "failed"` for the existing low-assignment fixture;
- threshold findings include `assignment_rate`;
- backend summary records CPU authority and GPU eligibility/skipping;
- HTML report links to the reliability report.

- [ ] **Step 2: Run the test and verify RED**

Run:

```bash
DOTMATCH_LIB=$PWD/libdotmatch.dylib PYTHONPATH=$PWD/python python3 -m pytest \
  python/tests/test_assayspec.py::test_assay_run_count_reproduces_existing_crispr_fixture -q
```

Expected: fails because runtime reliability artifacts are not written.

- [ ] **Step 3: Implement runtime aggregation**

After all native steps and optional autopsy complete, build reliability from:

- `sample_qc.tsv` threshold rows;
- audit `audit_summary.json` safety;
- autopsy artifact presence and findings;
- manifest command exit codes;
- `docs/assay-evidence.json` boundary;
- backend eligibility from assay mode, metric, `k`, extract length, and `allow_gpu`.

Write reliability artifacts before the manifest so the manifest can include
their paths and the assay report can link to them.

- [ ] **Step 4: Run the test and verify GREEN**

Run the command from Step 2. Expected: pass.

## Task 4: Production Blocking And Exploratory Recording

**Files:**
- Modify: `python/tests/test_assayspec.py`
- Modify: `python/dotmatch/assayspec.py`

- [ ] **Step 1: Write failing tests for unsafe target behavior**

Create two small count specs with colliding targets:

- production profile with `fail_on_unsafe_targets = true` should stop after audit, return non-zero, write reliability artifacts, and not write counts;
- exploratory profile should complete assignment, write reliability artifacts, and record unsafe target finding without failing only because of preflight.

- [ ] **Step 2: Run tests and verify RED**

Run:

```bash
DOTMATCH_LIB=$PWD/libdotmatch.dylib PYTHONPATH=$PWD/python python3 -m pytest \
  python/tests/test_assayspec.py::test_assay_run_production_blocks_unsafe_targets_before_assignment \
  python/tests/test_assayspec.py::test_assay_run_exploratory_records_unsafe_targets_without_preflight_block -q
```

Expected: fails because unsafe audit currently warns and continues.

- [ ] **Step 3: Implement post-audit preflight blocking**

After each audit step, inspect `audit_summary.json`. If the configured radius is
unsafe and the profile is production with `fail_on_unsafe_targets = true`, write
reliability artifacts and manifest, print a clear error, and return exit code 2
before assignment commands run.

- [ ] **Step 4: Run tests and verify GREEN**

Run the command from Step 2. Expected: pass.

## Task 5: Docs And Full Verification

**Files:**
- Modify: `docs/assayspec.md`

- [ ] **Step 1: Update docs**

Document `[reliability]`, `[backend]`, generated artifacts, production versus
exploratory profiles, and CPU/GPU authority boundaries.

- [ ] **Step 2: Run targeted tests**

Run:

```bash
python3 -m py_compile python/dotmatch/assayspec.py
DOTMATCH_LIB=$PWD/libdotmatch.dylib PYTHONPATH=$PWD/python python3 -m pytest python/tests/test_assayspec.py -q
```

Expected: all pass.

- [ ] **Step 3: Run release verification**

Run:

```bash
make test cli-test python-test release-ready
```

Expected: all pass, including release gates.

- [ ] **Step 4: Commit**

Commit the implementation with:

```bash
git add python/dotmatch/assayspec.py python/tests/test_assayspec.py docs/assayspec.md docs/superpowers/plans/2026-05-21-closed-loop-assay-reliability.md
git commit -m "feat: add assay reliability artifacts"
```
