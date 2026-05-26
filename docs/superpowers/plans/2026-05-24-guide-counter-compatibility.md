# GuideCounter Compatibility Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a GuideCounter-compatible count entrypoint that accepts GuideCounter-style flags and emits GuideCounter-style counts, extended counts, and stats outputs.

**Architecture:** Implement a thin `guide-counter count` compatibility wrapper in the native CLI. The wrapper delegates assignment/counting to the existing `run_count` engine using GuideCounter semantics: Hamming, best-distance one-mismatch by default, exact mode with `--exact-match`, automatic offset detection by default, and sample-name inference. After the count engine writes intermediate DotMatch/MAGeCK outputs, the wrapper rewrites the requested output prefix into GuideCounter-compatible files.

**Tech Stack:** C11 CLI, existing FASTQ counting engine, shell CLI regression tests.

---

### Task 1: GuideCounter CLI Red Test

**Files:**
- Modify: `tests/test_cli_fastq.sh`

- [x] **Step 1: Add a failing GuideCounter compatibility fixture**

Add a fixture invoking `dotmatch guide-counter count --input ... --library ... --output ...` with inferred sample names, exact and one-mismatch reads, automatic offset detection, and guide annotations.

- [x] **Step 2: Verify red**

Run: `make cli-test`
Expected: fail because the command is not implemented.

### Task 2: Compatibility Wrapper

**Files:**
- Modify: `src/qda.c`

- [x] **Step 1: Parse GuideCounter flags**

Support `guide-counter count`, `guide-counter-count`, and `guide-count` with:
`--input/-i`, `--samples/-s`, `--library/-l`, `--output/-o`, `--exact-match/-x`, `--offset-sample-size/-N`, `--offset-min-fraction/-f`, `--essential-genes/-e`, `--nonessential-genes/-n`, `--control-guides/-c`, and `--control-pattern/-C`.

- [x] **Step 2: Delegate to `run_count`**

Build an internal `count` argv with `--format mageck`, `--metric hamming`, `--ambiguity-policy best`, `--auto-offset 499`, and the requested offset sampling options.

- [x] **Step 3: Emit GuideCounter output files**

Create `{output}.counts.txt`, `{output}.extended-counts.txt`, and `{output}.stats.txt`.
The extended file inserts `guide_type` after `gene`; stats include GuideCounter-compatible mapped/mean/zero-read fields.

### Task 3: Verification

**Files:**
- Modify: `tests/test_cli_fastq.sh`
- Modify: `src/qda.c`

- [x] **Step 1: Verify green**

Run: `make cli-test`, `make test`, and `git diff --check` for touched files.
