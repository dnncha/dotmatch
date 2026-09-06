# DotMatch

**Count your guides. Account for every read.**

DotMatch turns FASTQ reads and a known guide library into count tables and
assignment QC. Use it for **CRISPR guide counting**, fixed-position **barcode
demultiplexing**, and other short-DNA assays with known targets. It runs locally
on Linux and macOS and writes MAGeCK-compatible counts.

[![CI](https://github.com/dnncha/dotmatch/actions/workflows/ci.yml/badge.svg)](https://github.com/dnncha/dotmatch/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/dotmatch?label=PyPI)](https://pypi.org/project/dotmatch/)
[![Documentation](https://readthedocs.org/projects/dotmatch/badge/?version=latest)](https://dotmatch.readthedocs.io/en/latest/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20541628.svg)](https://doi.org/10.5281/zenodo.20541628)

[Start counting guides](https://dnncha.github.io/dotmatch/crispr-guide-counting/) ·
[Check a library in your browser](https://dnncha.github.io/dotmatch/tools/library-safety/) ·
[Documentation](https://dotmatch.readthedocs.io/en/latest/) ·
[Methods and results](https://dotmatch.readthedocs.io/en/latest/benchmarks/crispr_comparison/README.html)

## Why use DotMatch?

An assigned-read percentage cannot tell you which targets gained counts, which
reads fit several targets, or whether a more permissive matching rule changed the
result. DotMatch keeps those decisions inspectable.

Every read has an explicit outcome: **unique**, **ambiguous**, **none** (unmatched),
or **invalid** (the requested window could not be extracted). Only unique calls
contribute to a target count. Choose the matching policy deliberately; a unique
call is not proof of biological origin.

Keep downstream screen statistics in the workflow you already use. DotMatch is not a genome aligner,
basecaller, cell/UMI pipeline or gene-level hit-calling package.

## Install

Release 0.5.0 includes the six `dotmatch agent` tools
described below:

```bash
python3 -m pip install dotmatch==0.5.0
dotmatch --version
```

Conda and container routes:

```bash
conda create -n dotmatch -c conda-forge -c bioconda dotmatch
conda activate dotmatch

# Or use the pinned release container:
docker run --rm ghcr.io/dnncha/dotmatch:v0.5.0 --version
```

Bioconda and its generated BioContainers images can lag PyPI/GHCR. When a
newly tagged version has not reached Bioconda yet, use PyPI or the source build.
Check the installed version. The Bioconda recipe includes `osx-arm64` for Apple Silicon.
Review the [packaging details](https://dotmatch.readthedocs.io/en/latest/packaging.html)
for platform and container verification. See the [installation guide](https://dotmatch.readthedocs.io/en/latest/getting-started.html)
for platform details, source builds and the third-party Homebrew tap. The optional
[desktop Workbench](https://github.com/dnncha/dotmatch-community) is maintained separately.

## Count a CRISPR screen

Prepare a guide CSV/TSV and FASTQ files. Start a new assay project:

```bash
dotmatch crispr quickstart \
  --library guides.csv \
  --fastq 'fastqs/*.fastq.gz' \
  --out crispr_screen/
```

This creates a draft project. Review `crispr_screen/inference_report.json` and
`assay.toml`: confirm the guide window, orientation, library and sample files.
After confirming the settings, change the top-level `status = "draft"` to
`status = "ready"` in `assay.toml`, then run and review:

```bash
dotmatch assay start crispr_screen/assay.toml

# After reviewing a completed run:
dotmatch assay handoff crispr_screen/assay.toml
```

The handoff carries configuration, QC, methods and checksums without copying raw
FASTQs. Follow the [complete CRISPR tutorial](https://dotmatch.readthedocs.io/en/latest/tutorials/crispr-count-first-run.html)
for inputs, direct CLI options and count-table outputs.

## Understand the effect of mismatch correction

The `dotmatch sensitivity` command, introduced in 0.5.0, compares exact, radius-one and
best-distance Hamming assignment using the same windows in **one FASTQ pass**.
It produces three count matrices, per-guide deltas, read-state transitions,
checksums and a self-contained HTML report. It never selects a policy for you.

Run the included synthetic example from a checkout of the v0.5.0 release:

```bash
python3 -m pip install dotmatch==0.5.0
dotmatch sensitivity \
  --targets examples/assignment_sensitivity/targets.tsv \
  --reads examples/assignment_sensitivity/reads.fastq \
  --target-start 0 --target-length 20 \
  --write-read-changes --out-dir sensitivity-example
```

The [nine-read synthetic example](https://github.com/dnncha/dotmatch/tree/main/examples/assignment_sensitivity)
shows why equal assigned totals can hide different per-guide counts.
[Read the output contract](https://dotmatch.readthedocs.io/en/latest/sensitivity.html).
This is sensitivity analysis, not an estimate of biological accuracy.

## Choose by task

| Task | Entry point | Workflow |
| --- | --- | --- |
| CRISPR guide counting | `dotmatch crispr-count` | [First run](https://dotmatch.readthedocs.io/en/latest/tutorials/crispr-count-first-run.html) |
| Inline barcode demultiplexing | `dotmatch demux` | [Getting started](https://dotmatch.readthedocs.io/en/latest/getting-started.html) |
| High unmatched or ambiguous barcode rate | `dotmatch barcode autopsy` | [Barcode diagnostics](https://dotmatch.readthedocs.io/en/latest/getting-started.html#diagnose-a-barcode-run) |
| Target-library collisions | `dotmatch audit` | [Browser checker](https://dnncha.github.io/dotmatch/tools/library-safety/) |
| Barcode panel design | `dotmatch panel design` | [Panel documentation](https://dotmatch.readthedocs.io/en/latest/barcode-panel-design.html) |
| Paired target counting | `dotmatch pair-count` | [Command reference](https://dotmatch.readthedocs.io/en/latest/command-reference.html) |
| Cell-by-feature matrix from extracted observations | `dotmatch feature matrix` | [scverse handoff](https://dotmatch.readthedocs.io/en/latest/tutorials/scverse-perturb-seq.html) |

Feature matrices require upstream cell identifiers and extracted feature windows.
They do not perform cell calling, UMI deduplication or perturbation-effect analysis.

## Reproduce the evidence

The [benchmark reports](https://dotmatch.readthedocs.io/en/latest/benchmarks/README.html)
include commands, hardware and assignment rules. Those reports cover the tested workloads;
they are not universal speed or biological-accuracy guarantees.

[Public CRISPR comparisons](https://dotmatch.readthedocs.io/en/latest/benchmarks/crispr_comparison/README.html)
record Yusa and Brunello inputs, methods, count differences, runtime and memory.
Exact, Hamming and Levenshtein results use different semantics and should be
compared separately. A comparison that completed successfully is not necessarily
an identical count matrix or biological validation.

The [GSE146194 direct-guide-capture case study](https://github.com/dnncha/dotmatch/tree/main/examples/perturb_seq_gse146194)
separates discovery and evaluation reads and checks per-read assignments against
independent reference implementations. It does not establish guide-per-cell or
perturbation-effect accuracy.

For an installation-free synthetic smoke demo, use [Binder](https://mybinder.org/v2/gh/dnncha/dotmatch/main?labpath=demo.ipynb)
or [Google Colab](https://colab.research.google.com/github/dnncha/dotmatch/blob/main/demo.ipynb).
For a shareable or de-identified workflow evaluation, see the [public validation invitation](https://github.com/dnncha/dotmatch/issues/82).
Do not post private reads or unpublished guide libraries.

## Pipelines, Python and local agents

DotMatch provides a [Python streaming API](https://dotmatch.readthedocs.io/en/latest/streaming-api.html),
[output schemas](https://dotmatch.readthedocs.io/en/latest/schemas.html), and
[workflow examples](https://github.com/dnncha/dotmatch/tree/main/examples/workflows).
The [ecosystem status ledger](https://dotmatch.readthedocs.io/en/latest/ecosystem-status.html)
distinguishes local examples from accepted upstream integrations.

The six structured agent tools are included in release 0.5.0:

```bash
dotmatch capabilities --json
dotmatch agent tools --json
dotmatch agent export-skill --target ./dotmatch-agent
```

They prepare, preflight, run, review and hand off local assays without accepting
free-form shell commands or uploading research data. Start with the
[Agent guide](https://dotmatch.readthedocs.io/en/latest/agent-guide.html),
[CRISPR agent route](https://dotmatch.readthedocs.io/en/latest/agent-crispr.html), or
[Perturb-seq agent route](https://dotmatch.readthedocs.io/en/latest/agent-perturb-seq.html).
Machine-readable discovery: [agent-capabilities.json](https://dnncha.github.io/dotmatch/agent-capabilities.json),
[agent-tools.json](https://dnncha.github.io/dotmatch/agent-tools.json), and the
[checked contract fixture](https://dnncha.github.io/dotmatch/agent-reference-crispr.json).

## Citation and contributing

DotMatch is [Apache-2.0 licensed](https://github.com/dnncha/dotmatch/blob/main/LICENSE).
Use `dotmatch citation` and [CITATION.cff](https://github.com/dnncha/dotmatch/blob/main/CITATION.cff)
to record the actual software version. Use the [methods and citation guide](https://dotmatch.readthedocs.io/en/latest/methods-and-citation.html)
to cite the actual release and configuration used. Improvements, discrepancy
fixtures and reproducible bug reports are welcome: [contributing guide](https://github.com/dnncha/dotmatch/blob/main/CONTRIBUTING.md).
