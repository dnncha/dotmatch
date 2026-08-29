# Agent guide

This page routes a sequencing task to the narrowest supported DotMatch command.
It is written for coding agents, scientific agents, workflow authors, and people
who want a copy-paste starting point without guessing from the full CLI surface.

DotMatch performs deterministic fixed-window known-target short-DNA assignment
from FASTQ. It requires a finite target list and a known or reviewed read
window. It is not a genome aligner, basecaller, adapter trimmer, variant caller,
cell/UMI quantifier, or downstream CRISPR screen-statistics package.

## Machine-readable routes

Packages that include the agent interface can print the versioned capability
record:

```bash
dotmatch capabilities --json
```

The same record is published as
[`agent-capabilities.json`](https://dnncha.github.io/dotmatch/agent-capabilities.json)
and validated against
[`agent-capabilities.schema.json`](https://dnncha.github.io/dotmatch/agent-capabilities.schema.json).
Use an intent's exact `entrypoint`, then read its `inputs`, `outputs`, and
`limitations` before constructing a command. DotMatch 0.3.0 and later include
the same capability record in the installed package.

## Choose by task

| Intent or search phrase | Entry point | Required decision | Important limit |
| --- | --- | --- | --- |
| CRISPR guide counting; MAGeCK-compatible counts | `dotmatch crispr-count` | Guide start, length, and correction radius | Counting only; no downstream screen statistics |
| Inline barcode demultiplexing; split FASTQ by barcode | `dotmatch demux` | Barcode start, length, and correction radius | Starts from FASTQ; no basecalling |
| Feature-barcode assignment; TotalSeq feature reads | `dotmatch count` | Feature window and known feature list | Per-read assignment; no cell/UMI or Cell Ranger quantification |
| Perturb-seq guide capture; CRISPR guide-capture reads | `dotmatch count` | Guide window and known guide list | Public evidence is single-guide extraction; no guide-per-cell or perturbation effects |
| Barcode panel design or collision checking | `dotmatch panel design` or `dotmatch panel check` | Panel size, length, preset, and safety radius | Short barcode sets only; not probe, primer, or full assay design |
| Known-target FASTQ matching; whitelist counting | `dotmatch count` | Fixed window, metric, and target list | Finite known targets only; not general alignment |
| High unmatched or ambiguous barcode rate | `dotmatch barcode autopsy` | Plausible offset range and candidate `k` values | Suggestions require assay-context review |

## Install in a clean environment

```bash
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install dotmatch==0.3.0
dotmatch --version
```

Published Python wheels currently cover x86_64 and aarch64 Linux with glibc or
musl and macOS 11 or newer on Apple Silicon or Intel. Windows wheels are not
published.
Bioconda can lag PyPI, so inspect `dotmatch --version` when a workflow pins an
exact release.

## Small fixed-window workflow

From a clone of the repository, the committed fixture provides a complete
FASTQ count without downloading assay data:

```bash
mkdir -p smoke-output
dotmatch count \
  --targets demo-data/crispr_guides.tsv \
  --reads demo-data/reads.fastq \
  --sample-label smoke \
  --target-start 0 \
  --target-length 4 \
  --k 0 \
  --metric hamming \
  --out smoke-output/counts.tsv \
  --summary smoke-output/summary.json
```

The checked fixture has three reads. At `k=0`, one read is a unique exact match
and two are unmatched. The packaging gate builds the source distribution and
wheel, installs each into a fresh virtual environment, runs this equivalent
workflow, and checks the summary and counts.

## Safe defaults

- Coordinates are zero-based.
- Start with `k=0` unless correction is required.
- Use Hamming distance for equal-length windows and substitutions.
- Use Levenshtein only when short insertion or deletion rescue is intended.
- Before `k>=1`, run `dotmatch audit` on the same target list and radius.
- Count `unique` assignments only. Keep `ambiguous`, `none`, and `invalid`
  visible in assignments and QC.
- Request `--summary` and, when practical, `--assignments` for provenance and
  diagnosis.

## Error recovery

For a high unmatched rate, inspect recurring windows:

```bash
dotmatch inspect-unmatched \
  --targets targets.tsv \
  --reads sample.fastq.gz \
  --target-start 0 \
  --target-length 20 \
  --k 0 \
  --top 50 \
  --out top_unmatched.tsv
```

Check the read side, start, length, orientation, trimming, and target table.

For a high ambiguous rate or uncertain correction safety:

```bash
dotmatch audit \
  --targets targets.tsv \
  --k 1 \
  --audit-mode auto \
  --out-dir audit
```

Lower `k` or redesign colliding targets when the audit is unsafe. For many
`invalid` reads, confirm that start plus length fits the reads after trimming.

## Evidence and public boundaries

Each intent in the capability manifest names repository tests, checked example
outputs, or public benchmark records. Those files support only their recorded
conditions. In particular:

- the public feature-barcode lane supports fixed-window per-read assignment,
  not cell or UMI quantification;
- the public guide-capture lane contains one observed guide and supports
  extraction validation, not useful multi-guide or perturbation-effect claims;
- maintained workflow examples are not accepted external integrations unless
  [`workflow-adoption.json`](workflow-adoption.json) records one;
- package retrieval counts are download events, not unique users or adoption.

See [Scope and limitations](trust-and-scope.md),
[Scientific claims](scientific-claims.md), and [Output schemas](schemas.md)
before broadening a workflow or a claim.

## Why there is no MCP server

DotMatch already has a local, scriptable CLI, documented file contracts, and a
machine-readable capability command. The repository does not currently show a
consumer that needs a long-running tool server or remote protocol boundary.
Adding MCP would increase installation and trust surface without improving the
checked local workflow, so the supported agent interface remains the CLI and
ordinary files.
