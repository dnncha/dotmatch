# EditWitness

**Know what your CRISPR assay can—and cannot—see.**

EditWitness finds explicit genomic alternatives that produce the same modeled
observations as an intended edit. It then asks which of your proposed follow-up
assays can separate those alternatives.

Keep your aligner and editing-analysis pipeline. Add a check on what their
measurements are capable of establishing.

> **Research alpha · 0.1.0a1.** This is working, software-tested scientific
> software with an intentionally bounded, idealized model. It has **not** been
> empirically validated. It does not diagnose an editing defect, estimate its
> probability, or certify a clone as correctly edited or safe.

## The distinction that matters

An amplicon containing only the intended sequence may be compatible with two
very different states:

```text
intended + intended                 intended + allele invisible to this assay
          │                                        │
          └──────────── the same sequence-presence observation ────────────┘
```

In the model, detecting the intended sequence does not establish the state of a
second allele that emits no signal. A witness is the **specific alternative**,
the reason it remains indistinguishable, and the candidate measurements that
would distinguish it under stated assumptions.

This is an observability check, not a new variant caller.

## Try the synthetic example

Requires Python 3.11 or newer. From this source directory:

```bash
python -m venv .venv
# macOS / Linux:
source .venv/bin/activate
# Windows PowerShell instead: .venv\Scripts\Activate.ps1
python -m pip install .

editwitness demo -o demo.json
editwitness analyze demo.json -o analysis.json --html report.html
editwitness analyze demo.json --compact --pretty
editwitness verify analysis.json --manifest demo.json
```

Open `report.html` locally. It needs no server, JavaScript, account, API key, or
external assets. The example is synthetic, not a published experiment.

**Expected result:** the existing inner amplicon leaves two declared alternatives
indistinguishable from the intended biallelic hypothesis. The planner selects the
outer amplicon, which separates the primer-deletion alternative. Deletion of the
entire local window remains unresolved by every supplied candidate.

```json
{
  "conclusion": "ambiguity_demonstrated",
  "equivalent_alternatives": ["hidden_primer_deletion", "hidden_window_deletion"],
  "selected_assays": ["outer"],
  "unresolved_hypotheses": ["hidden_window_deletion"]
}
```

This excerpt flattens selected fields for readability. The actual versioned
output nests panel information under `plan`.

**This release is not published on PyPI.** Install the source or the accompanying
wheel. Do not assume that an unrelated package with a matching name is ours.

## What works today

| Capability | What it actually does |
|---|---|
| Explicit counterexamples | Compares declared two-allele clonal hypotheses against the expected hypothesis. |
| Sequence-aware readout | Reconstructs local replacements and models full primer-trimmed inserts or ordered paired-end reads. |
| Read-gap ambiguity | Does not pretend that paired-end reads observe the unsequenced middle or product length. |
| Follow-up panel selection | Finds a minimum-cost panel for up to 18 useful candidates; uses an explicitly nonoptimal greedy method above that limit. |
| Local deletion scan | Streams a bounded deletion grid using interval geometry without reconstructing every sequence. |
| Human-readable evidence | Generates a self-contained HTML report with coordinate diagrams and inspectable sequence evidence. |
| Reproducible artifacts | Emits canonical JSON, normalized-input checksums, versioned assumptions, and replay verification. |
| Agent interface | Provides schemas, capability discovery, structured errors, stable exit codes, and compact summaries. |

### Start from a real local reference

`init` accepts one local FASTA record, exact 5′→3′ primer oligos, and one intended
substitution. It resolves unique inward-orientation matches within that window.
It **does not** establish genome-wide primer specificity.

```bash
editwitness init --fasta locus.fasta \
  --left-primer YOUR_FORWARD_OLIGO --right-primer YOUR_REVERSE_OLIGO \
  --edit-position 450 --alternate A -o design.json
```

Replace the illustrative values with your own data. `450` is a **local,
zero-based** position, not a VCF position. The alternate must differ from the
reference. Inspect and extend the generated hypotheses and assay assumptions
before drawing conclusions. JSON supports more general, nonoverlapping local
replacements; see the [data contract](docs/data-contract.md).

### Explore a read-gap blind spot

```bash
editwitness demo --paired-end -o paired.json
editwitness analyze paired.json --compact --pretty
editwitness witness paired.json --hypothesis interior_deletion --include-sequences
```

The paired-end example leaves four alternatives indistinguishable. The outer
full-insert candidate separates three; whole-window loss still remains unresolved.

### Audit a deletion grid

```bash
editwitness scan demo.json -o deletion-scan.json
editwitness verify deletion-scan.json --manifest demo.json
```

These are counts on a **declared grid**, not event probabilities, assay
sensitivity, or a biological frequency distribution. The scan is separate from
the hypothesis-comparison engine and operates on the reference, not the intended
allele.

## For Python users

```python
from editwitness import analyze, load_manifest

manifest = load_manifest("design.json")
result = analyze(manifest)
for witness in result.witnesses:
    print(witness.hypothesis_id, witness.resolving_candidate_assays)
print(result.plan.unresolved_hypotheses)
```

The public entry points share the CLI's model. No LLM participates in the
scientific calculation.

## For agents and workflow authors

```bash
editwitness capabilities
editwitness schema manifest
editwitness validate design.json
editwitness analyze design.json --compact
```

Stdout is JSON; errors are JSON on stderr. Successful execution can legitimately
report ambiguity. To use ambiguity as a workflow gate, explicitly add
`--fail-on-ambiguity` (exit code 4).

Read the [agent guide](docs/agent-guide.md), the portable
[agent skill](skills/editwitness/SKILL.md), and [AGENTS.md](AGENTS.md). An agent
must never translate exit code 0, a valid manifest, or an empty list of declared
counterexamples into “safe” or “biallelic confirmed.”

## The boundaries are part of the product

The current model requires **pristine original annotated primer sites** and
perfect detection of eligible products. It does not model mismatch tolerance,
new or rescued binding sites, stochastic PCR, amplification efficiency, read
counts, copy-number measurements, mosaicism, aneuploidy, or complex structural
rearrangements. It consumes a design manifest, not FASTQ/BAM data.

Changing a supplied hypothesis set can change the answer. “Distinguishable within
the declared model” does not mean “all biological alternatives were excluded.”

Read the [scientific model](docs/scientific-model.md) before applying it to an
experiment. Every full analysis carries these qualifications with the evidence.

## Documentation

- [Scientific model and references](docs/scientific-model.md)
- [Coordinates, inputs, outputs, and limits](docs/data-contract.md)
- [Architecture and extension decisions](docs/architecture.md)
- [Agent integration and exit codes](docs/agent-guide.md)
- [Validation status and independent-validation plan](docs/validation.md)
- [Prioritized roadmap](docs/roadmap.md) and [machine-readable tasks](roadmap.json)
- [Development](CONTRIBUTING.md), [release procedure](docs/releasing.md), and [security](SECURITY.md)

## Development

```bash
python -m pip install -e '.[dev]'
python -m pytest -q
python scripts/generate_schemas.py --check
python scripts/check_style.py
python -m mypy src/editwitness
python -m build
```

Only the checks listed in [BUILD_STATUS.md](BUILD_STATUS.md) were actually run
for the delivered artifact. CI configuration is not evidence of a passing CI run.

## License and attribution

Apache-2.0. Copyright 2026 Donncha O'Toole and EditWitness contributors.
See [LICENSE](LICENSE) and [CITATION.cff](CITATION.cff). No paper, DOI, clinical
validation, or institutional endorsement is claimed for this alpha.
