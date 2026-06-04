# Stale Branch Integration Audit

Date: 2026-06-04

## Summary

The old `codex/performance-hotpath` and `codex/evidence-gallery` branches contain useful ideas, but they are not safe integration branches. Both are based on old release/doc states and direct merges would revert current Read the Docs, JOSS, release metadata, dependency, and readiness work.

Keep mining individual ideas only when current `origin/main` lacks an equivalent implementation. Do not merge either branch wholesale.

## `codex/performance-hotpath`

Useful themes:

- same-length low-edit fast paths;
- hamming seed/candidate indexes;
- stricter CLI input validation and allocation checks;
- release/readiness gate hardening.

Current status:

- Superseded by current `origin/main` through newer changes including `6cdf82b` (`perf: fast-path same-length low-edit comparisons`), `d297245` (`Validate BCL numeric inputs`), `9698140` (`Harden CLI numeric parsing and BCL writes`), and later release-readiness work.
- Direct diff from current main would delete newer docs and JOSS files and revert release metadata.

Decision:

- Do not merge.
- Local branch can be pruned after this audit.

## `codex/evidence-gallery`

Useful themes:

- AssaySpec reliability artifacts;
- backend optimizer and CPU/GPU advisory model;
- GPU/Metal benchmark lane;
- GuideCounter-compatible counting surface;
- richer evidence/gallery release materials.

Current status:

- Substantial evolved versions are already on current main: `python/dotmatch/assayspec.py` includes reliability and backend optimizer support, `src/qda.c` includes GuideCounter-compatible and hamming strategy paths, and public evidence/readiness documentation has advanced past this branch.
- Direct diff from current main would delete newer Read the Docs, JOSS, release, workflow, and benchmark gate files.
- Remote branch `origin/codex/evidence-gallery` still exists, so the reference remains available even if the local branch is pruned.

Decision:

- Do not merge.
- Keep remote reference as archive unless explicitly closing/deleting remote branches.
- Local branch can be pruned after this audit.

## Remaining Actionable Work

- Consolidated current, safe work lives on `codex/consolidate-repo-review`.
- If future performance work is desired, start from current `origin/main` and benchmark targeted changes against the existing hamming seed/direct lookup paths.
- If future evidence UX work is desired, start from current AssaySpec reliability outputs and add narrow tests before changing public claims or benchmark artifacts.
