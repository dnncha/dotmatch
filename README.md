# EditWitness

**What could your CRISPR validation assay miss?**

EditWitness compares explicit genomic alternatives against a sequencing assay's measurement design. It returns concrete counterexamples: different local genotypes that would produce the same modeled observations. It then identifies which supplied follow-up assays separate those alternatives—and which remain unresolved.

Keep your existing caller. Use EditWitness to inspect what the measurement can establish.

**Research alpha · Apache-2.0 · local-first · Python 3.11+**

This is an assay-design and model-inspection tool, not a variant caller, a PCR simulator, or a certificate that an edited clone is safe. Software tests are not experimental validation.

## Try the central example

Install from the public source branch. The repository location is temporary: EditWitness is an independent project hosted on an isolated branch of `dnncha/dotmatch`; **do not merge it into DotMatch**. Standalone repository and PyPI publication are not complete.

```bash
python -m pip install "git+https://github.com/dnncha/dotmatch.git@editwitness/public"
editwitness demo -o experiment.json
editwitness analyze experiment.json -o evidence.json --html report.html
editwitness verify evidence.json --manifest experiment.json
```

For reproducible work, replace the branch with the full audited commit SHA. The public build records that SHA with its distribution checksums.

The synthetic example compares an intended biallelic edit with alternatives including a deletion that removes a primer site. The baseline assay cannot distinguish two alternatives from the expectation. An outer assay separates one; loss of the whole modeled window remains unresolved by every supplied sequence-presence assay.

**That is an ambiguity in the measurement, not evidence that either deletion occurred.**

## What is different about this package?

| Question | EditWitness output |
|---|---|
| Could a different local genotype produce the same readout? | An explicit alternative, its edits, and per-assay sequence evidence. |
| Does the answer depend on how primer binding is modeled? | A comparison of original-site eligibility and exact local sequence rematching. |
| Would short paired-end reads distinguish an interior change? | A comparison using only the configured read ends—not the unsequenced gap or latent product length. |
| Which additional assays help? | A minimum-cost panel for small candidate sets, a labeled heuristic for larger sets, and an explicit unresolved list. |
| Did duplicate allele names inflate the evidence? | Sequence-equivalent diploid genotypes are grouped; aliases do not create extra witnesses. |

## Two explicit observation models

`editwitness demo` and `editwitness init` select **`exact-local-sites-presence-v2`**. It reconstructs each edited allele, rematches the two primer oligos in both inward-facing orientations, and retains every exact local product within the declared size bounds. Replacement representations producing the same allele sequence therefore produce the same modeled sequence signals.

The historical **`original-sites-presence-v1`** remains available for comparison. It requires pristine annotated sites and can give representation-dependent answers. Older manifests that omit `observation_model` retain that historical behavior; they are not silently migrated.

```bash
editwitness compare-models experiment.json --pretty
```

Neither model predicts amplification efficiency, mismatch tolerance, sampling failures, or genome-wide specificity. A modeled difference is not automatically detectable experimentally. See the [scientific model](docs/scientific-model.md).

## Start from your own local reference

Primer oligos are written 5′ to 3′. Positions are local, **zero-based and half-open**. FASTA import accepts one A/C/G/T reference record of at most 20,000 bases. The initial substitution must lie between the annotated primer sites.

```bash
editwitness init \
  --fasta locus.fasta \
  --left-primer "$FORWARD_OLIGO" \
  --right-primer "$REVERSE_OLIGO" \
  --edit-position 450 --alternate T \
  --deletion-radius 300 --deletion-step 50 \
  -o design.json

editwitness expand-deletions design.json -o alternatives.json
editwitness analyze alternatives.json -o evidence.json --html report.html
```

The example position and alternate base are placeholders, not a proposed biological edit. Check them against your local reference.

Deletion expansion creates a bounded reference-haplotype deletion grid paired with one expected allele. It records the grid, source digest, event count and duplicate count. It does not assign event probabilities, choose heterozygous phase, or silently truncate alternatives to fit a limit. Candidate follow-up assays are supplied explicitly; this release does not synthesize new oligos or calculate primer thermodynamics.

## A usable interface for agents

```bash
editwitness capabilities
editwitness schema manifest
editwitness validate alternatives.json
editwitness analyze alternatives.json --compact
editwitness witness experiment.json --hypothesis hidden_primer_deletion --include-sequences
```

The final command uses a hypothesis from the bundled demo; select an ID from your own result for other manifests.

JSON goes to standard output; structured errors go to standard error. There are no prompts, progress decorations, telemetry, remote reference fetches, API keys or language-model calls. Full evidence and compact summaries are separate outputs. Files are not overwritten without `--force`, and input files are protected even then.

A successful process exit means the calculation completed. It does **not** mean that a clone passed validation. Use `--fail-on-ambiguity` to return exit code 4 after writing a result containing counterexamples.

See the [agent guide](docs/agent-guide.md), [machine-readable capabilities](docs/capabilities.json), [portable skill](skills/editwitness/SKILL.md), and [data contract](docs/data-contract.md).

## Python API

```python
from editwitness import analyze, compare_models, expand_deletions, load_manifest

manifest = load_manifest("experiment.json")
result = analyze(manifest)
for witness in result.witnesses:
    print(witness.hypothesis_id, witness.resolving_candidate_assays)

sensitivity = compare_models(manifest)
print(sensitivity["original_only"], sensitivity["exact_only"])
```

Public analysis functions revalidate their input, including model instances constructed using unchecked Pydantic helpers. Deterministic outputs include model/package/schema versions, input hashes, explicit assumptions, and replayable evidence.

## Scope and trust

The package currently supports a finite collection of local replacement alleles and diploid clonal hypotheses. It does not infer genotypes from BAM/FASTQ, interpret read fractions as genomic dosage, model mosaic samples, or estimate clinical risk. Whole-genome events and sites outside the supplied reference window are not assessed.

The deletion **geometry** scanner is deliberately separate from sequence-aware equivalence analysis. Its counters use the original-site model and are labeled accordingly, even when the manifest selects exact rematching for `analyze`.

Read the [audit](docs/audit-0.2.md), [validation plan](docs/validation.md), and [build evidence](BUILD_STATUS.md). No independent biological validation has yet been performed.

## Development and contribution

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
python scripts/generate_schemas.py --check
python scripts/check_style.py
python -m mypy src/editwitness
python -m build
```

Start with [AGENTS.md](AGENTS.md) and the [architecture](docs/architecture.md). The [roadmap](roadmap.json) gives bounded tasks and acceptance criteria. Experimental cases with known assay geometry and independently established outcomes are especially valuable; never post identifiable or restricted biological data in public issues.

Please cite the exact software version or commit using [CITATION.cff](CITATION.cff). This project does not yet have a paper or DOI.
