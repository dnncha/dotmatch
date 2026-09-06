# EditWitness

**Know what your CRISPR assay can—and cannot—see.**

EditWitness finds distinct, explicitly described genomic states that produce the
same modeled assay observations. It shows the counterexample, then compares
additional assays that could distinguish it. Keep your existing caller: this is
an assay-design and evidence-inspection tool, not a replacement aligner.

**Research alpha · Apache-2.0 · Python 3.11+ · Local execution · No API key**

Version **0.2.0a2** includes exact primer rematching against final edited DNA,
bounded deletion-hypothesis generation, complete allele evidence, and an offline
installation self-test. Heterozygous challenge generation requires an explicit
choice of the allele to preserve. It is
software-tested, **not empirically validated for biological or clinical use**.
See [build evidence](BUILD_STATUS.md), [the model audit](docs/audit-0.2.0a1.md),
and [release hardening](docs/audit-0.2.0a2.md).
A prepared release is not a published one: current distribution status is
recorded in `BUILD_STATUS.md`. Do not assume a PyPI release exists.

## The problem in one example

An amplicon shows only the intended sequence. Does that establish that both
alleles carry the edit? Not from sequence presence alone: an intended allele
paired with an unobserved deleted allele can give the same modeled result.

The bundled synthetic example makes that ambiguity inspectable:

| Assay panel | What the model can distinguish |
|---|---|
| Original inner assay | The expected edit remains equivalent to two distinct alternatives. |
| Add the supplied outer assay | The primer-deletion alternative becomes distinguishable. |
| Both supplied assays | The whole-window-deletion alternative remains unresolved. |

The output does **not** say that a deletion occurred, estimate its frequency,
or declare a clone safe. It names the alternatives and assumptions behind the
comparison. [Read the scientific model](docs/scientific-model.md).

## Install the research alpha

The public release location is temporarily
[`dnncha/dotmatch`, tag `editwitness-v0.2.0a2`](https://github.com/dnncha/dotmatch/releases/tag/editwitness-v0.2.0a2).
This is **an independent package**, not a DotMatch release or dependency. Its
source lives on an isolated branch; do not merge that branch into DotMatch.
Repository creation is the remaining administrative step for `dnncha/editwitness`.
There is no PyPI release. Use only the versioned asset after the prerelease is
published; a source branch or a queued CI job is not release evidence.

```bash
python -m pip install "https://github.com/dnncha/dotmatch/releases/download/editwitness-v0.2.0a2/editwitness-0.2.0a2-py3-none-any.whl"
editwitness self-test
```

The release includes `SHA256SUMS`, the clean source distribution, a synthetic
HTML report and execution provenance tied to the tested source commit. The
self-test checks two synthetic cases and deterministic replay offline. It is
**not** a validation of a real assay.

## Try the example

Alternatively, from a downloaded source distribution, in its root directory:

```bash
python -m pip install .
editwitness self-test
editwitness demo -o demo.json
editwitness analyze demo.json -o analysis.json --html report.html
editwitness verify analysis.json --manifest demo.json
```

Open `report.html`. It is a self-contained report with no scripts, CDN,
telemetry, or external font requests. The JSON retains the full evidence;
the HTML deliberately limits long sequence previews.

For paired-end sequencing, the tool observes only the declared sequenced ends,
not the unsequenced gap or an implicitly known product length:

```bash
editwitness demo --paired-end -o paired.json
editwitness analyze paired.json -o paired-analysis.json --html paired.html
```

For a real design, start with `editwitness init --help`. Initialization requires
an explicit choice of `--full-insert` or `--read-bases N`. All coordinates are
local, **zero-based, half-open**, on the supplied reference. Generated templates
select `exact-local-sequence-presence-v2` explicitly.

## What is implemented

**Sequence-aware counterexamples.** Exact local matching finds preserved,
recreated and new sites on final edited DNA. It enumerates inward-facing
heteroprimer products in both orientations. Multiple products remain explicit;
identical final diploid sequence states are not counted as alternative genomes.

**Useful, bounded hypothesis generation.** Supply an explicit `deletion_scan`
grid, then generate sequence-deduplicated diploid alternatives:

```bash
editwitness expand-deletions design.json -o expanded.json
editwitness analyze expanded.json -o expanded-analysis.json --html expanded.html
```

Generation records its input hash, grid, fixed allele and counts. When the
expected alleles encode different final sequences, specify `--fixed-allele ID`;
reordering the allele list must not silently select a different challenge. It fails rather
than silently sampling if the configured limit is exceeded. These are challenges
to the assay design, **not predicted CRISPR outcome frequencies**.

**Follow-up panel selection.** Given candidate assays and integer costs, select
measurements that separate the expectation from the currently equivalent
alternatives. Small nondominated panels use exhaustive optimization; larger
ones use a declared greedy method. Impossible-to-separate alternatives remain
unresolved. No claim of globally optimal primer design is made.

**Assumption sensitivity.** `editwitness compare-models design.json` shows which
counterexamples depend on the old original-site model versus the exact local
sequence model. This comparison is not experimental adjudication.

## For agents and workflows

```bash
editwitness self-test
editwitness capabilities
editwitness schema manifest
editwitness validate design.json
editwitness analyze design.json --compact
editwitness witness design.json --hypothesis hidden_primer_deletion --include-sequences
```

JSON is the primary interface. Diagnostics go to stderr; successful commands do
not print prose into machine output. Versioned schemas, explicit resource
limits, content hashes and exact-version replay make runs inspectable.
An exit code of zero means the calculation completed, not that an edit is safe.

```python
from editwitness import analyze, expand_deletions, load_manifest

manifest = load_manifest("design.json")
result = analyze(manifest)
for witness in result.witnesses:
    print(witness.hypothesis_id)
```

[Agent guide](docs/agent-guide.md) · [Data contract](docs/data-contract.md) ·
[Architecture](docs/architecture.md) · [Migration from 0.1](docs/migration-0.2.md)

## Limits worth reading

The exact model assumes that every eligible exact local product is detected.
It does not predict PCR efficiency, mismatch tolerance, sequencing errors,
read counts, allele competition, dosage, mosaicism, or genome-wide specificity.
Single-primer products and sites outside the supplied sequence window are not
modeled. Identical oligos are rejected because their read orientation is
ambiguous under this contract. The set of supplied hypotheses is finite, not an
exhaustive account of biology.

`scan` remains a separately labeled **original-site geometry scanner**. Its
counts are not exact-model read equivalence, biological risk or sensitivity.
Use `expand-deletions` plus `analyze` for generated sequence-model comparisons.

## Contributing and validation

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
python -m mypy src/editwitness
python -m ruff check src/editwitness --select E4,E7,E9,F
python scripts/generate_schemas.py --check
python scripts/release_manifest.py --check
```

The tests include an independent substring-based observation oracle, equivalent
edit representations, multisite/reverse-orientation products, old artifact
integrity, and panel-selection comparisons against independent enumeration.
These establish software behavior—not biological accuracy.

[Validation plan](docs/validation.md) · [Contribution rules](CONTRIBUTING.md) ·
[Machine-readable roadmap](roadmap.json) · [Next-session handoff](docs/continuation.md)

Independent genome-engineering review and adjudicated experimental benchmarks
are the next scientific gates. No reviewers, users, benchmark accuracy or
clinical validity are claimed without actual evidence. Citation metadata is in
[`CITATION.cff`](CITATION.cff); there is no registered publication DOI yet.
