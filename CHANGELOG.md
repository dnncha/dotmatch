# Changelog

## 0.2.0a2 — release hardening

- Require an explicit fixed allele for heterozygous deletion generation; allele
  ordering no longer silently chooses the preserved haplotype.
- Bound repeated deletion reconstruction and total hypothesis signal references;
  exhausted budgets raise errors rather than silently truncating evidence.
- Add `self-test`: two bundled synthetic scenarios, content integrity and replay,
  JSON output, and a distinct failure exit code. No network or test dependencies.
- Verify installed wheels and independently extracted source archives in CI.
- Expand the CI matrix and test the declared minimum Pydantic dependency.
- Prepare namespaced public prerelease distribution without changing DotMatch main.
  Actual publication and exact-commit CI must be checked in the release record.


## 0.2.0a1 — 2026-09-06

Sequence-aware exact local model; representation-invariant evidence; identical-state
exclusion; complete missing-allele evidence; bounded deletion generation; explicit
readout initialization; public API revalidation; safe panel dominance; legacy
artifact integrity; strengthened packaging and CI-gated prerelease publisher.
See `docs/audit-0.2.0a1.md` and `docs/migration-0.2.md`.


## 0.1.0a1 — 2026-09-05

Initial research alpha. Added a strict local manifest, original-primer-site
observation model, full-insert and paired-end sequence-presence comparisons,
explicit alternative-hypothesis witnesses, bounded exact/greedy candidate-panel
selection, streaming deletion geometry scans, FASTA initialization, reproducible
JSON and replay verification, a script-free HTML report, portable agent guidance,
synthetic examples, independent-oracle tests, packaging and a prioritized roadmap.

No raw-read calling, empirical performance claims, copy-number response model,
clinical validation, or automatic package-index publication is included.
