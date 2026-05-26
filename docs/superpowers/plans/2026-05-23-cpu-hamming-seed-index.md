# CPU Hamming Seed Index Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a memory-efficient CPU seed lookup path for core indexed Hamming assignment so each k=1 read uses two seed probes instead of enumerating every possible one-base mutation.

**Architecture:** Extend `qdaln_index` with a two-seed lookup table for encodable DNA targets of length 2-32. Keep existing exact lookup and fallback mutation enumeration intact, and dispatch the new seed path only for `qdaln_index_assign_hamming_stats(..., k=1, ...)`.

**Tech Stack:** C11, existing `qdalign.c` index internals, Makefile native tests.

---

### Task 1: Seed Index Data Structure

**Files:**
- Modify: `src/qdalign.c`
- Test: `tests/test_qdalign_threshold_alloc.c`

- [ ] **Step 1: Write the failing test**

Add an internal test that builds a mixed target panel and asserts the index owns a populated two-seed Hamming table.

- [ ] **Step 2: Run test to verify it fails**

Run: `make build/test_qdalign_threshold_alloc`
Expected: compile failure because `hamming_seed_ready`, `n_hamming_seeds`, and `hamming_seed_hash_cap` do not exist yet.

- [ ] **Step 3: Add seed fields and builder**

Add `qdaln_hamming_seed_entry`, seed table fields to `qdaln_index`, a hash/insert helper, and populate two seed entries per encodable target with length 2-32.

- [ ] **Step 4: Run test to verify it passes**

Run: `make build/test_qdalign_threshold_alloc && ./build/test_qdalign_threshold_alloc`
Expected: pass.

### Task 2: Seed-Based Hamming Assignment

**Files:**
- Modify: `src/qdalign.c`
- Test: `tests/test_qdalign.c`

- [ ] **Step 1: Route k=1 Hamming assignment through the seed table**

Add a seed visitor that deduplicates candidates with `candidate_seen`, verifies candidate Hamming distance with the packed distance helper, and preserves `match_count`, `best_distance`, `second_best_distance`, and ambiguity semantics.

- [ ] **Step 2: Run native regression tests**

Run: `make test`
Expected: both native test binaries pass.

### Task 3: Verification and Benchmark

**Files:**
- Modify: no additional files

- [ ] **Step 1: Run CLI regression tests**

Run: `make cli-test`
Expected: FASTQ and CRISPR CLI smoke tests pass.

- [ ] **Step 2: Run a focused before/after microbenchmark**

Compare `HEAD~1` against the working tree for `qdaln_index_assign_hamming_stats` at guide-like lengths and report reads/sec plus checksum agreement.

- [ ] **Step 3: Commit and push**

Stage only `src/qdalign.c`, `tests/test_qdalign.c`, `tests/test_qdalign_threshold_alloc.c`, and this plan if changed. Commit as `perf: add core hamming seed index`, then push `codex/evidence-gallery`.
