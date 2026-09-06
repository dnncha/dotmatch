# Public-release hardening audit — 0.2.0a2

6 September 2026. Implementation self-review, not independent scientific review.
The separately archived 0.2.0a1 audit remains part of the evidence record.

## EW-A11: allele ordering silently chose the preserved haplotype

For a heterozygous expectation, deletion generation defaulted to the first allele
identifier. Reversing a biologically unordered pair could therefore change the
challenge set without the user making a different experimental assumption.

**Fix:** require an explicit fixed allele whenever the two expected final DNA
sequences differ. Sequence-identical aliases remain acceptable and choose a
stable identifier. Regression tests cover both input orders, explicit choices,
invalid types and sequence-identical aliases.

## EW-A12: bounded product enumeration did not bound repeated evidence expansion

The exact model limited products, but each product signal could be referenced by
many hypotheses. Highly repetitive, otherwise valid inputs could amplify the
output and computation far beyond the original product budget. Deletion
sequence deduplication similarly did not bound repeated reconstruction work.

**Fix:** cumulative limits of 100,000 hypothesis-signal references and 200 million
reconstructed deletion bases. Budget exhaustion raises a structured input error;
no silent sampling, partial result, or reassuring zero is emitted. Boundary tests
check the exact accepted threshold and the first rejected request.

## EW-A13: prepared distributions did not demonstrate a working installation

A build and source-tree tests do not establish that the installed artifact can
load packaged fixtures, schemas or its actual runtime modules. The previous
candidate had no new remote type-check or cross-platform evidence.

**Fix:** add a dependency-free-after-install `self-test` command, exercise it in
the isolated wheel smoke, and run the tests and source inventory from an
independently extracted source distribution. CI gates all release artifacts on
tests, typing, lint, schema checks and coverage. Actual run outcomes belong in
the linked CI record, not in claims made before execution.

The test module `tests/test_release_hardening.py` adds 13 regression cases. The
local total after the initial implementation was 630 passing tests. Refer to
`BUILD_STATUS.md` and the exact-commit release evidence for subsequent checks.

## Still unproven

No biological sensitivity estimate, independent laboratory review, complete
repair-outcome catalogue, validated copy-number model, clinical use, or adoption
claim is supported by this pass. Public software distribution does not change
those scientific limits.
