# Performance Domination Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the next high-impact CPU, GPU, and ML optimization lanes without changing scientific assignment semantics.

**Architecture:** Split work into four disjoint lanes: core library Hamming radius expansion in `src/qdalign.c`, CLI counting throughput in `src/qda.c`, benchmark-informed strategy routing in Python, and Metal benchmark/prototype improvements in `tools/`. Every lane must keep CPU deterministic assignment as the authority and prove equality with existing scan/oracle paths.

**Tech Stack:** C11, pthreads, Objective-C++/Metal benchmarks, Python 3 tests, existing Makefile test targets.

---

### Task 1: Core CPU Hamming k=2/k=3 Seed Index

**Files:**
- Modify: `include/qdalign.h`
- Modify: `src/qdalign.c`
- Modify: `tests/test_qdalign.c`
- Modify: `tests/test_qdalign_threshold_alloc.c`

- [x] **Step 1: Write failing tests for indexed Hamming k=2 and k=3**

Add tests comparing `qdaln_index_assign_hamming_stats` against `qdaln_match_many` for exact, duplicate, ambiguous, no-match, non-ACGT fallback, and lengths 4, 20, and 32.

- [x] **Step 2: Verify red**

Run: `make test`
Expected: fail because indexed Hamming currently rejects `k > 1`.

- [x] **Step 3: Implement segmented seed tables**

Extend the existing seed index to support `k=1..3` by storing `k + 1` deterministic seed segments per encodable target and verifying candidates with packed Hamming distance. Keep scan fallback for unsupported lengths/alphabet.

- [x] **Step 4: Verify green**

Run: `make test`
Expected: all qdalign tests pass.

### Task 2: CLI Counting Sparse Merge

**Files:**
- Modify: `src/qda.c`
- Modify: `tests/test_cli_fastq.sh`

- [x] **Step 1: Add regression coverage**

Add a threaded large-target hamming count case that compares serial and threaded count TSV output and validates summary JSON.

- [x] **Step 2: Implement dirty-slot local count merging**

Track local count slots touched by each worker and merge only those dirty slots instead of sweeping `threads * targets * 5` counters.

- [x] **Step 3: Verify**

Run: `make cli-test`
Expected: all CLI smoke tests pass.

### Task 3: ML/Stat Strategy Router

**Files:**
- Modify: `python/dotmatch/assayspec.py`
- Modify: `python/tests/test_assayspec.py`

- [x] **Step 1: Add optimizer tests**

Add tests that `optimize_assay_backend` reports recommended CPU strategy, thread hints, benchmark prior count, and route reasons for Hamming, Levenshtein, GPU-eligible, and GPU-ineligible assays.

- [x] **Step 2: Implement benchmark-informed route metadata**

Extend optimizer output beyond static GPU candidate fields to include CPU strategy, threading hint, diagnostics constraints, GPU crossover notes, and confidence from benchmark prior coverage.

- [x] **Step 3: Verify**

Run: `PYTHONPATH=python python3 -m pytest python/tests/test_assayspec.py -q`
Expected: pass.

### Task 4: GPU Seed-Index Benchmark Prototype

**Files:**
- Modify: `tools/bench_gpu_metal.mm`
- Modify: `tools/bench_gpu_crispr_metal.mm`
- Modify: `scripts/bench_gpu.py`
- Modify: `scripts/bench_gpu_crispr.py`

- [x] **Step 1: Add benchmark mode and equality checks**

Add a benchmark row identifying whether the GPU path is brute-force scan or seed-index prototype. Keep CPU checksum/mismatch validation mandatory.

- [x] **Step 2: Prototype GPU seed-index candidate reduction**

For fixed A/C/G/T Hamming k=1, add Metal-side seed candidate lookup or a clearly isolated CPU-built seed candidate buffer consumed by Metal. The prototype must preserve exact match_count and ambiguity semantics for benchmark cases.

- [x] **Step 3: Verify**

Run: `make build/bench_gpu_metal`
Expected: build succeeds on macOS with Metal available. If Metal execution is unavailable, the binary must report unavailable rather than fail.

### Task 5: Integration

**Files:**
- Modify: benchmark docs only if fresh benchmark evidence is generated.

- [x] **Step 1: Run full verification**

Run: `make test`, `make cli-test`, and targeted Python tests for touched modules.

- [x] **Step 2: Benchmark accepted lanes**

Run focused before/after microbenchmarks for any CPU/GPU lane that lands.

- [x] **Step 3: Commit and push scoped changes**

Stage only files changed by accepted lanes. Do not stage unrelated dirty benchmark/docs/example files already present in the worktree.
