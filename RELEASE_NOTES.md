# EditWitness 0.2.0a2 — public research alpha

**Know what your CRISPR assay can—and cannot—see.**

This is an independent EditWitness package temporarily hosted under a namespaced
prerelease in `dnncha/dotmatch`. It is not a DotMatch release, dependency or change
to DotMatch main. Do not merge the isolated EditWitness branch into DotMatch.

## What it does

Reconstruct explicit local genomic alternatives, model exact inward-facing
heteroprimer products in both orientations, and show alternatives indistinguishable
from the expected edit under the declared readout. Compare follow-up assay panels
while preserving unresolved cases. Full-insert and paired-end observations remain
separate; copy number and unsequenced gaps are not inferred.

## Improvements since the initial alpha

Exact final-sequence rematching replaces representation-dependent reasoning when
explicitly selected. Complete alternative allele definitions, bounded deletion
hypothesis generation, model comparisons, strict input checks, offline reports,
JSON schemas, integrity checks and replay are included.

This hardening release also requires explicit haplotype selection for heterozygous
expectations, bounds repeated evidence expansion, and adds `editwitness self-test`.
The public artifacts are gated on the linked CI run, including installed-wheel
and extracted-source tests. `release-evidence.json` identifies the tested commit.

## Install and check

Download the matching wheel and verify its hash against `SHA256SUMS`, then:

```bash
python -m pip install ./editwitness-0.2.0a2-py3-none-any.whl
editwitness self-test
editwitness demo -o demo.json
editwitness analyze demo.json -o analysis.json --html report.html
editwitness verify analysis.json --manifest demo.json
```

Python 3.11 or later is required. The first dependency installation may use the
network; the computation and bundled self-test do not. Apache-2.0 licensed.

## Evidence and limits

The example HTML is synthetic. The software is not empirically biologically
validated, does not diagnose defects, does not estimate actual PCR performance
or editing-outcome probabilities, and does not certify a clone as safe. A passing
software self-test concerns the executable, not any real experiment.

No PyPI release is claimed. Standalone repository migration, independent
scientific review and an adjudicated biological benchmark remain on the embedded
roadmap. See the clean source archive for model details, audits and agent guidance.
