# ML Backend Optimizer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a benchmark-informed AssaySpec backend optimizer that recommends fast CPU/GPU execution plans while preserving deterministic CPU assignment authority.

**Architecture:** The optimizer lives in `python/dotmatch/assayspec.py` beside existing backend summaries. It extracts workload features from an AssaySpec, scores CPU and GPU candidates with conservative benchmark priors, writes `backend_optimization.json`, and surfaces the result in reliability summaries.

**Tech Stack:** Python 3.9+, existing AssaySpec CLI, JSON artifacts, pytest.

---

### Task 1: Optimizer Data Model And Summary

**Files:**
- Modify: `python/dotmatch/assayspec.py`
- Test: `python/tests/test_assayspec.py`

- [ ] **Step 1: Write failing tests**

Add tests that import `optimize_assay_backend`, load the existing CRISPR count fixture, and assert:

```python
def test_backend_optimizer_recommends_gpu_candidate_for_public_crispr(tmp_path: Path) -> None:
    from dotmatch.assayspec import load_assay_spec, optimize_assay_backend

    assay = load_assay_spec(_write_count_spec(tmp_path))
    plan = optimize_assay_backend(assay)

    assert plan["authority"] == "cpu"
    assert plan["selected_backend"] == "cpu"
    assert plan["candidate_backend"] == "gpu-metal-experimental"
    assert plan["recommendation"] == "gpu_candidate_requires_cpu_validation"
    assert plan["expected_speedup_band"] == "1.5-3x"
    assert "public_gpu_gate_validated" in plan["reason_codes"]
    assert "cpu_count_checksum_required" in plan["accuracy_gates"]
```

- [ ] **Step 2: Verify the test fails**

Run:

```bash
python -m pytest python/tests/test_assayspec.py::test_backend_optimizer_recommends_gpu_candidate_for_public_crispr -q
```

Expected: fail because `optimize_assay_backend` does not exist.

- [ ] **Step 3: Implement minimal optimizer**

Add `optimize_assay_backend(assay: AssaySpec) -> dict[str, Any]`, plus small helpers for target features and speedup band. The first implementation should return CPU authority and a GPU candidate for eligible public CRISPR Hamming `k=1` specs.

- [ ] **Step 4: Verify the test passes**

Run:

```bash
python -m pytest python/tests/test_assayspec.py::test_backend_optimizer_recommends_gpu_candidate_for_public_crispr -q
```

Expected: pass.

### Task 2: Gating Cases

**Files:**
- Modify: `python/dotmatch/assayspec.py`
- Test: `python/tests/test_assayspec.py`

- [ ] **Step 1: Write failing tests**

Add tests for demux and Levenshtein:

```python
def test_backend_optimizer_gates_compute_compatible_demux_without_public_gpu_gate(tmp_path: Path) -> None:
    from dotmatch.assayspec import load_assay_spec, optimize_assay_backend

    assay = load_assay_spec(_write_demux_spec(tmp_path))
    plan = optimize_assay_backend(assay)

    assert plan["selected_backend"] == "cpu"
    assert plan["candidate_backend"] == "gpu-metal-experimental"
    assert plan["recommendation"] == "gpu_candidate_gated"
    assert "compute_compatible_no_public_gpu_gate" in plan["reason_codes"]


def test_backend_optimizer_requires_cpu_for_levenshtein(tmp_path: Path) -> None:
    from dotmatch.assayspec import load_assay_spec, optimize_assay_backend

    spec = _write_count_spec(tmp_path)
    spec.write_text(spec.read_text(encoding="utf-8").replace('metric = "hamming"', 'metric = "levenshtein"'), encoding="utf-8")
    assay = load_assay_spec(spec)
    plan = optimize_assay_backend(assay)

    assert plan["candidate_backend"] == "cpu"
    assert plan["recommendation"] == "cpu_required"
    assert "metric_not_gpu_supported" in plan["reason_codes"]
```

- [ ] **Step 2: Verify the tests fail**

Run the two tests by name with `python -m pytest ... -q`. Expected: fail before gating logic exists.

- [ ] **Step 3: Implement gating reason codes**

Extend optimizer helpers so unsupported metrics, edit radii, variable target length, non-A/C/G/T targets, pair-count mode, and missing public GPU gates produce explicit reason codes.

- [ ] **Step 4: Verify the tests pass**

Run the two targeted tests. Expected: pass.

### Task 3: CLI Artifact

**Files:**
- Modify: `python/dotmatch/assayspec.py`
- Test: `python/tests/test_assayspec.py`

- [ ] **Step 1: Write failing CLI test**

Add:

```python
def test_assay_optimize_writes_backend_optimization_artifact(tmp_path: Path) -> None:
    spec = _write_count_spec(tmp_path)

    rc = _run_cli(["assay", "optimize", str(spec)])

    assert rc.returncode == 0, rc.stderr
    assert "gpu-metal-experimental" in rc.stdout
    artifact = tmp_path / "assay_out" / "backend_optimization.json"
    data = json.loads(artifact.read_text(encoding="utf-8"))
    assert data["authority"] == "cpu"
    assert data["candidate_backend"] == "gpu-metal-experimental"
```

- [ ] **Step 2: Verify the CLI test fails**

Run the test by name. Expected: fail because `optimize` is not a recognized subcommand.

- [ ] **Step 3: Add `dotmatch assay optimize`**

Add an `optimize` subparser, compile the plan, call `optimize_assay_backend`, write `backend_optimization.json`, and print a concise summary.

- [ ] **Step 4: Verify the CLI test passes**

Run the targeted test. Expected: pass.

### Task 4: Reliability Integration

**Files:**
- Modify: `python/dotmatch/assayspec.py`
- Test: `python/tests/test_assayspec.py`

- [ ] **Step 1: Write failing reliability test**

Extend `test_assay_check_writes_preflight_reliability_artifacts` to assert:

```python
assert summary["backend_optimizer"]["authority"] == "cpu"
assert summary["backend_optimizer"]["candidate_backend"] == "gpu-metal-experimental"
```

- [ ] **Step 2: Verify the reliability test fails**

Run the test by name. Expected: fail because `backend_optimizer` is missing.

- [ ] **Step 3: Add optimizer metadata to reliability summary**

Add `backend_optimizer: optimize_assay_backend(plan.spec)` to `_build_reliability_summary`.

- [ ] **Step 4: Verify the reliability test passes**

Run the targeted test. Expected: pass.

### Task 5: Full Verification And Commit

**Files:**
- Modify: all files from previous tasks.

- [ ] **Step 1: Run focused tests**

```bash
python -m pytest python/tests/test_assayspec.py -q
```

Expected: all tests in `test_assayspec.py` pass.

- [ ] **Step 2: Run formatting/whitespace check**

```bash
git diff --check
```

Expected: no output.

- [ ] **Step 3: Stage only optimizer files**

```bash
git add python/dotmatch/assayspec.py python/tests/test_assayspec.py docs/superpowers/specs/2026-05-23-ml-backend-optimizer-design.md docs/superpowers/plans/2026-05-23-ml-backend-optimizer.md
```

- [ ] **Step 4: Commit and push**

```bash
git commit -m "feat: add assay backend optimizer"
git push origin codex/evidence-gallery
```
