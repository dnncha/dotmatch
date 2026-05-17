# Barcode Troubleshooting Failure Fixtures

These tiny fixtures document the failure classes the barcode troubleshooting report is
expected to explain. They are intentionally synthetic, because the purpose is
diagnostic coverage rather than public-data performance evidence.

The public evidence lanes live in `docs/barcode-science-readiness.json`.
These fixtures cover the failure vocabulary:

- wrong offset;
- duplicate barcode;
- unsafe one-edit collision;
- ambiguous read;
- unmatched low-complexity read;
- low-quality correction rejected;
- invalid extraction window;
- reverse-complement candidate.

The expected classifications are in `expected_findings.tsv`.
