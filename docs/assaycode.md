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
| AssayCode Watch | Reserved identity for future streaming run intelligence; not a current capability |
| AssayCode Pro | Commercial workbench, registry, signed-report, private-assay, audit-trail, and support boundary |

Names in this table are architecture and positioning. A surface is not a public
capability until its command, schema, tests, evidence, and release gate exist.

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

AssayCode currently organizes existing assay-level capabilities. It does not yet
claim a general assay compiler, calibrated probabilistic decoder, production
streaming monitor, or universal sequencing platform.

Those capabilities require implementation plus independent correctness,
calibration, performance, and public-data evidence before their names become
claims.
