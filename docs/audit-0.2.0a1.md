# Audit of EditWitness 0.1.0a1 and changes in 0.2.0a1

**Review date: 6 September 2026.** This is an implementation self-audit, not an
independent scientific endorsement. The original 572-test baseline passed locally;
the audit therefore targeted assumptions and missing adversarial cases rather
than treating the existing test count as proof of correctness.

The two most consequential findings were edit-representation dependence and
identical genomic states being reported as alternatives. The new exact-local
model and state-identity checks address those issues. It still describes a
bounded idealized measurement, not observed PCR performance.

## EW-A01 · HIGH — Edit notation changed the scientific conclusion

**Finding.** A broad replacement retaining the original primer DNA was classified as disrupted by the original-site model, while a compact edit yielding the same final sequence was observable.

**Change.** Added explicitly versioned exact final-sequence rematching with both orientations, multisite products and representation-invariance tests. Legacy remains named and warned, not silently redefined.

**Evidence.** `tests/test_exact_model.py::test_sequence_invariance_rescues_unchanged_site_in_broad_replacement`

## EW-A02 · HIGH — Aliases could inflate genomic counterexamples

**Finding.** Renamed hypotheses or alternate allele IDs encoding the same final diploid sequence state could be counted as distinct alternatives.

**Change.** Compare actual unordered final sequence pairs; label identical local states and exclude them from witnesses in both models.

**Evidence.** `tests/test_exact_model.py::test_same_state_with_different_allele_ids_is_not_an_alternative`

## EW-A03 · HIGH — Readout initialization silently assumed complete observation

**Finding.** A new CLI manifest could default to full-insert readout without a deliberate user declaration of actual sequence coverage.

**Change.** Require mutually exclusive --full-insert or --read-bases in CLI initialization. Model configuration remains visible.

**Evidence.** `tests/test_audit_regressions.py::test_cli_init_never_silently_assumes_full_insert_observation`

## EW-A04 · HIGH — Public API could trust bypassed model validation

**Finding.** Frozen Pydantic instances can still be constructed/copied without validation; the public engine accepted those instances.

**Change.** Enable instance revalidation and validate at public computational boundaries, including nested copied models.

**Evidence.** `tests/test_exact_model.py::test_public_engine_revalidates_model_copy_and_construct`

## EW-A05 · MEDIUM — The missing allele was not fully inspectable in result evidence

**Finding.** An allele with no signal could be named without its complete edit definition in a result/report. Sequence inspection emphasized observed reads rather than the hidden DNA.

**Change.** Store every allele edit definition, final length and SHA; emit final sequences through focused witness inspection. Add candidate evidence to HTML.

**Evidence.** `tests/test_audit_regressions.py::test_complete_witness_includes_reconstructable_final_alleles`

## EW-A06 · MEDIUM — Geometry scan was not a hypothesis-generating workflow

**Finding.** The alpha required manual encoding of useful alternatives; structural scan counters did not establish read/genotype equivalence.

**Change.** Add bounded reference-deletion generation, fixed-allele declaration, sequence deduplication, explicit provenance and fail-closed caps; keep geometry scanner separate.

**Evidence.** `tests/test_audit_regressions.py::test_cli_generated_challenges_are_actual_analysis_inputs`

## EW-A07 · MEDIUM — Redundant candidate assays could needlessly disable exact optimization

**Finding.** Candidate-count threshold switched to greedy even when many candidates were safely dominated duplicates.

**Change.** Perform deterministic objective-preserving dominance reduction and expose removed candidates. Retain impossible-to-separate cases.

**Evidence.** `tests/test_audit_regressions.py::test_equivalent_candidate_dominance_enables_exact_search`

## EW-A08 · MEDIUM — New schema defaults could invalidate old artifact checksums

**Finding.** Hashing the newly validated model rather than the archived payload would add defaults and change the digest.

**Change.** Verify actual archived JSON before migration defaults; old artifacts pass integrity. Require originating executable for exact replay.

**Evidence.** `tests/test_audit_regressions.py::test_old_result_integrity_survives_new_optional_schema_fields`

## EW-A09 · HIGH — Release status and source staging were not actual independent publication

**Finding.** The earlier release was staged in an unrelated repository, and current-session GitHub operations cannot create or push a standalone repository.

**Change.** Keep DotMatch untouched. Provide a fresh-history publisher with identity, visibility, source and exact-commit CI gates. Record publication as blocked, not completed.

**Evidence.** `tests/test_release_safety.py; docs/releasing.md`

## EW-A10 · MEDIUM — Source-distribution contents omitted new regression fixtures

**Finding.** Packaging patterns included Python tests but not their JSON fixtures, workflow files or full documentation/evidence inventory.

**Change.** Explicitly include JSON fixtures, documentation JSON, workflows, line-ending rules and release notes. Test the source archive independently.

**Evidence.** `MANIFEST.in; source-distribution verification recorded in BUILD_STATUS.md`

## Evidence scope and next release gates

See `BUILD_STATUS.md` and `docs/verification.json` for checks actually completed
in the current environment. The historical alpha's remote six-platform CI,
coverage and typing results are preserved under `docs/history/`, not recycled as
proof about this candidate. Local package builds are not a public release.

Public repository/prerelease creation is blocked by unavailable authenticated
writes in this session. The local publication helper is tested with offline
preflights/mocks, not a live GitHub creation test. Strict typing and Ruff must be
run by the configured remote CI before releasing its artifacts.

Independent scientific review, a provenance-complete biological benchmark and a
real caller adapter remain uncompleted. No clinical safety, accuracy, sensitivity,
independent users, laboratory reviewers or PyPI publication is asserted.
