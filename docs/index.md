# DotMatch

**Known-target DNA read assignment and QC from FASTQ.**

DotMatch assigns a fixed window from each FASTQ read to a known table of short DNA sequences: CRISPR guides, inline barcodes, feature tags, primers, adapter prefixes, or panel targets. It runs locally on Linux and macOS and keeps ambiguous, unmatched, and invalid reads visible instead of collapsing them into target counts.

DotMatch is the guide-counting engine for **CRISPRWorks Count**; see the [CRISPRWorks family overview](crisprworks.md) for component status and experimental gene-level inference tools.

:::{admonition} Benchmark spotlight: faster one-mismatch counting with identical counts
:class: tip

In paired benchmarks on the source checkout, `dotmatch guide-counter count` ran **1.6–5.5× faster with 11.3–14.4× lower peak memory** than unmodified `guide-counter` 0.1.3 on tested one-mismatch workflows, producing identical full guide-by-sample count matrices.

Controlled 100k-read and 1M-read FASTQs, the public 87,437-guide Yusa library, one CPU thread, and five paired repeats define this result. `guide-counter` is faster in the tested exact-mode 1M-read cases; full experimental-screen performance and biological accuracy require separate validation.
[Explore the graphs, protocol, and raw results](benchmarks/guide_counter/README.md).
:::

For every read, DotMatch records one of four outcomes:

- `unique`: exactly one target fits the configured rule (added to target counts by default);
- `ambiguous`: two or more targets fit;
- `none`: no target lies within the selected distance;
- `invalid`: the configured read window could not be extracted.

## Choose the inputs and assignment rule

A reliable run starts with a known target library, a defined read window, and an explicit matching policy. Review those choices before interpreting a count matrix: a `unique` assignment reflects compatibility under the selected rules rather than proof of biological origin.

## Install and verify

Install the current PyPI package:

```bash
python3 -m pip install dotmatch
dotmatch --version
```

Release 0.8.0 is published on PyPI (GHCR has 0.7.0) and includes the six `dotmatch agent` tools and the local `dotmatch mcp` server described in the [Agent guide](agent-guide.md). Confirm that `dotmatch --version` reports `0.8.0` before relying on those tools; Bioconda and BioContainers can lag PyPI.

Before running on study data, verify your installation offline with `dotmatch demo --out-dir first-run/`, then follow [Getting started](getting-started.md) for your first count or demultiplexing run. Use the [Command reference](command-reference.md) to look up flags, or run `dotmatch agent tools --json` to inspect the structured tool contracts in your installed version.

## Choose a workflow

| I want to… | Start here |
| --- | --- |
| Count guides from a CRISPR screen | [CRISPR guide-counting first run](tutorials/crispr-count-first-run.md) |
| Fit gene effects from reviewed counts | [From counts to gene effects](tutorials/crispr-fit-first-run.md) |
| Reproduce a public CRISPR example | [Public CRISPR guide-counting example](../examples/crispr_guides/README.md) |
| Reproduce the multi-guide Perturb-seq case study | [GSE146194 direct-guide-capture report](benchmarks/perturb_seq_gse146194/README.md) |
| Try the workflow without a local install | [Binder](https://mybinder.org/v2/gh/dnncha/dotmatch/main?labpath=demo.ipynb) or [Google Colab](https://colab.research.google.com/github/dnncha/dotmatch/blob/main/demo.ipynb) |
| Compare existing count tables with DotMatch | [Count comparison](count-comparison.md) |
| Compare guide-counting workflows | [Guide-counting workflow comparison](usability-comparison.md) |
| Build a checked assay project | [AssaySpec workflows](assayspec.md) |
| Evaluate and hand off a run | [Lab evaluation and handoff](lab-evaluation.md) |
| Lock an independent comparison before results | [Independent evaluation packet](independent-evaluation-packet.md) |
| Split reads by inline barcode | [Getting started: demultiplexing](getting-started.md#demultiplex-inline-barcodes) |
| Assign feature-barcode reads | [Feature-barcode and guide-capture tutorial](tutorials/scverse-perturb-seq.md#count-guide-or-feature-barcode-reads) |
| Assign Perturb-seq guide-capture reads | [Feature-barcode and guide-capture tutorial](tutorials/scverse-perturb-seq.md) |
| Diagnose barcode failures | [Barcode run diagnosis](getting-started.md#diagnose-a-barcode-run) |
| Design or check a barcode panel | [Barcode panel design](barcode-panel-design.md) |
| Route an agent to a checked command | [Agent guide](agent-guide.md) |
| Run CRISPR counting with a local agent | [Agent CRISPR route](agent-crispr.md) |
| Run direct-guide Perturb-seq capture with a local agent | [Agent Perturb-seq route](agent-perturb-seq.md) |
| Build a cell-by-feature matrix from extracted observations | [scverse and feature barcodes](tutorials/scverse-perturb-seq.md) |
| Use DotMatch from Python | [Streaming Python API](streaming-api.md) |
| Add DotMatch to a pipeline | [Output schemas](schemas.md) |
| Evaluate DotMatch for a workflow | [Bioinformatics evaluation](bioinformatics-evaluation.md) |
| Check package and upstream integration status | [Ecosystem status](ecosystem-status.md) |
| Run an independent evaluation | [10x evaluation request](10x-evaluation-request.md) |
| Review the public Yusa evaluation draft | [Yusa independent evaluation request](yusa-independent-evaluation-request.md) |
| Review the 10x evidence gates | [10x goal baseline](10x-goal-baseline.md) |
| Record the software in a methods section | [Methods and citation](methods-and-citation.md) |

## Scope and limitations

DotMatch compares fixed read windows with a finite target list. It is not a genome aligner, basecaller, UMI pipeline, variant caller, or downstream screen-analysis package.

Performance and correctness results apply to the commands, datasets, and hardware recorded in the [benchmark reports](benchmarks/README.md). The strongest evidence covers native fixed-window assignment, public CRISPR guide-counting comparisons, checked inline-barcode lanes, and the bounded GSE146194 multi-guide direct-capture case study. Narrower experimental paths and their limits are documented in [Scope and limitations](trust-and-scope.md).

## Related research

[Cheerful Duck Research](https://cheerfulduck.com/research) collects our bioinformatics investigations, reproduction materials, and corrections. Cite the specific DotMatch release and assignment policy used in your analysis via the methods guide below.

## Help and citation

- Report a bug or request a feature on [GitHub Issues](https://github.com/dnncha/dotmatch/issues).
- Run `dotmatch citation` to print the installed release citation.
- See [Methods and citation](methods-and-citation.md) for `CITATION.cff`, DOI metadata, and methods wording.
- See [Packaging](packaging.md) for PyPI, Bioconda, containers, and source builds.

```{toctree}
:maxdepth: 2
:caption: Getting started

getting-started
agent-guide
agent-crispr
agent-perturb-seq
command-reference
sensitivity
tutorials/crispr-count-first-run
tutorials/crispr-fit-first-run
tutorials/scverse-perturb-seq
lab-evaluation
independent-evaluation-packet
count-comparison
```

```{toctree}
:maxdepth: 2
:caption: Workflows

assayspec
crispr-qc
barcode-panel-design
usability-comparison
```

```{toctree}
:maxdepth: 2
:caption: APIs and integration

streaming-api
python-dataframes
schemas
workbench
bioinformatics-evaluation
external-review-packet
agent-discovery-audit
ecosystem-status
10x-evaluation-request
yusa-independent-evaluation-request
10x-goal-baseline
10x-performance-control
```

```{toctree}
:maxdepth: 2
:caption: Reference

crisprworks
trust-and-scope
scientific-reliability
scientific-audit
benchmarks/README
benchmarks/hamming_seed_layout/README
benchmarks/guide_counter/README
methods-and-citation
packaging
```
