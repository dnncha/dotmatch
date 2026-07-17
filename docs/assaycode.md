# AssayCode Platform

AssayCode is the assay-level product identity built on the DotMatch deterministic
known-target assignment engine.

> **AssayCode expresses and operates the assay. DotMatch performs and validates
> the sequence assignment.**

This is an additive brand architecture. It does not rename or replace the
published `dotmatch` package, CLI, native library, output schemas, DOI, or
citation.

## Product Architecture

| Surface | Responsibility |
| --- | --- |
| AssayCode | Platform for specifying, validating, running, and diagnosing known-target sequencing assays |
| AssayScript | Human-reviewable declarative assay specification |
| DotMatch Engine | Native indexed assignment, ambiguity accounting, counting, demultiplexing, and validation |
| AssayCode Audit | Pre-run collision, input, and safety checks |
| AssayCode Autopsy | Evidence-backed diagnosis of offsets, unsafe correction, unmatched reads, and other failures |
| AssayCode Design | Barcode-panel design, optimization, simulation, layout, and export |
| AssayCode Watch | Experimental bounded-memory sequential QC over assignment event streams |
| AssayCode Pro | Commercial workbench, registry, signed-report, private-assay, audit-trail, and support boundary |

Names in this table are architecture and positioning. A surface is not a public
capability until its command, schema, tests, evidence, and release gate exist.

## AssayScript v2 Compilation

AssayScript v2 describes multiple target-bearing segments across R1, R2, I1,
and I2. Each segment declares a fixed position or anchor, target library,
orientation, jitter, metric, radius, ambiguity policy, and requirement status.
An optional constraint table declares allowed cross-segment combinations.

```bash
assaycode compile assay-v2.toml --out assay.plan.json
assaycode inspect assay.plan.json
```

The compiler validates and fingerprints every input, audits bounded target
sets, records review findings, selects a deterministic strategy per segment,
and emits a portable JSON plan. It does not yet claim a production universal
multi-segment execution runtime.

## Commands

Installing the Python distribution provides both command identities:

```bash
dotmatch --version
assaycode --version
```

The existing command remains authoritative and fully supported. AssayCode adds
shortcuts for assay-level work:

```bash
assaycode check assay.toml
assaycode plan assay.toml
assaycode run assay.toml
assaycode start assay.toml
```

These are exact convenience routes to the corresponding DotMatch AssaySpec
commands. They do not alter matching semantics or output contracts:

```text
assaycode check assay.toml
    == dotmatch assay check assay.toml
```

Specialized namespaces pass through unchanged:

```bash
assaycode crispr quickstart --library guides.csv --fastq '*.fastq.gz' --out run/
assaycode barcode autopsy --barcodes barcodes.tsv --reads reads.fastq.gz --out-dir autopsy/
assaycode panel check barcodes.tsv --k 1 --metric hamming --out-dir panel_check/
```

The explicit engine escape hatch is useful in scripts and documentation:

```bash
assaycode engine dist ACGT AGGT
assaycode engine validate --targets targets.tsv --reads reads.fastq.gz \
  --target-length 20
```

## Design-time simulation

Use the experimental digital twin to stress-test a fixed-length target panel
under a declared substitution-error rate before sequencing:

```bash
assaycode simulate \
  --targets targets.tsv \
  --reads-per-target 10000 \
  --error-rate 0.01 \
  -k 1 \
  --seed 42 \
  --out simulation.json
```

The JSON result reports correct and incorrect unique calls, ambiguous calls,
no-calls, usable yield, false-discovery rate, and a truth-by-call confusion
table. Results are deterministic for a seed. This substitution-only simulator
is an experimental design aid, not a replacement for held-out empirical
validation or a platform performance claim.

## Streaming QC

`assaycode watch` consumes assignment events as JSON Lines. It keeps bounded
state, emits periodic rate snapshots with a 95% Wilson interval, and returns
`insufficient_data`, `on_track`, or `review` decisions under explicit
thresholds.

```bash
assaycode watch assignments.jsonl --out watch.jsonl --every 100000
```

This is an experimental workflow primitive, not a sequencer-control or adaptive
sampling claim.

## Experimental Calibration

`dotmatch.calibration` contains a held-apart experimental decoder with
per-cycle empirical error fitting, Phred shrinkage, selective posterior and
likelihood-ratio thresholds, joint inference over allowed assay combinations,
Brier score, expected calibration error, and held-out FDR threshold selection.
Deterministic DotMatch assignment remains the production default.

A file-backed evaluation path makes this testable without integrating it into
the production FASTQ router:

```bash
assaycode calibrate trusted-observed-expected.tsv --out error-model.json
assaycode decode-quality \
  --reads held-out-windows.tsv \
  --targets targets.tsv \
  --model error-model.json \
  --out calls.tsv
```

Training rows must contain independently trusted `observed`, `expected`, and
`quality` columns. Calibration and false-discovery thresholds must be
evaluated on separate held-out truth data before any scientific claim.

## Python Identity

```python
import assaycode

assert assaycode.PLATFORM_NAME == "AssayCode"
assert assaycode.ENGINE_NAME == "DotMatch"
assert assaycode.__version__ == assaycode.engine.__version__

distance = assaycode.engine.distance("ACGT", "AGGT")
```

Scientific APIs continue to live under `dotmatch`. The AssayCode namespace
intentionally exposes the engine rather than silently copying its API and
creating two competing contracts.

## Compatibility Contract

The transition must preserve all of the following:

- `pip install dotmatch` and Bioconda package identity;
- `dotmatch` CLI commands and behavior;
- `import dotmatch` Python APIs;
- native C ABI and header/library artifacts;
- output schemas and provenance;
- DOI, CITATION.cff, and release citations;
- legacy `quickdna` compatibility where currently supported.

AssayCode may become the broader website or commercial identity, but scientific
methods should continue to cite the DotMatch release that performed assignment.

## Claim Boundary

AssayCode now includes an experimental multi-segment compiler, quality-aware
calibration module, joint tuple decoder, and streaming QC primitive. It does not
yet claim a production universal multi-segment runtime, calibrated public-data
superiority, sequencer control, or universal sequencing-platform coverage.

Those capabilities require implementation plus independent correctness,
calibration, performance, and public-data evidence before their names become
claims.
