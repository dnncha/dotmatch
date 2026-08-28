# DotMatch

DotMatch assigns a short window from each FASTQ read to a known list of DNA
sequences. It is made for CRISPR guides, sample barcodes, feature tags, primers,
adapter checks, and other assays where the expected targets are already known.

For every read, DotMatch records one of four outcomes:

- `unique`: exactly one target is compatible;
- `ambiguous`: several targets are compatible;
- `none`: no target is within the selected distance;
- `invalid`: the requested window could not be extracted.

Ambiguous and invalid reads stay visible in the outputs. They are not quietly
added to a target count.

## Start with a real run

Install the current PyPI package:

```bash
python3 -m pip install dotmatch
dotmatch --version
```

Then follow [Getting started](getting-started.md) for a small count or
demultiplexing run. If you already know which command you need, go straight to
the [command reference](command-reference.md).

## Choose a workflow

| I want to… | Start here |
| --- | --- |
| Count guides from a CRISPR screen | [CRISPR guide-counting first run](tutorials/crispr-count-first-run.md) |
| Reproduce a public CRISPR example | [Public CRISPR guide-counting example](../examples/crispr_guides/README.md) |
| Try the workflow without a local install | [Binder](https://mybinder.org/v2/gh/dnncha/dotmatch/main?labpath=demo.ipynb) or [Google Colab](https://colab.research.google.com/github/dnncha/dotmatch/blob/main/demo.ipynb) |
| Compare guide-counting workflows | [Guide-counting workflow comparison](usability-comparison.md) |
| Build a checked assay project | [AssaySpec workflows](assayspec.md) |
| Split reads by inline barcode | [Getting started: demultiplexing](getting-started.md#demultiplex-inline-barcodes) |
| Diagnose barcode failures | [Barcode run diagnosis](getting-started.md#diagnose-a-barcode-run) |
| Design or check a barcode panel | [Barcode panel design](barcode-panel-design.md) |
| Use DotMatch from Python | [Streaming Python API](streaming-api.md) |
| Add DotMatch to a pipeline | [Output schemas](schemas.md) |
| Evaluate DotMatch for a workflow | [Bioinformatics evaluation](bioinformatics-evaluation.md) |
| Check package and upstream integration status | [Ecosystem status](ecosystem-status.md) |
| Record the software in a methods section | [Methods and citation](methods-and-citation.md) |

## Scope and limitations

DotMatch compares fixed read windows with a finite target list. It is not a
genome aligner, basecaller, UMI pipeline, variant caller, or downstream screen
analysis package.

Performance and correctness results are tied to the commands, datasets, and
hardware recorded in the [benchmark reports](benchmarks/README.md). The most
developed paths are native fixed-window assignment, public CRISPR guide-counting
comparisons, and checked inline-barcode examples. Other assay types and
experimental backends have narrower test coverage; those limits are described
in [Scope and limitations](trust-and-scope.md).

## Help and citation

- Report a bug or request a feature on
  [GitHub Issues](https://github.com/dnncha/dotmatch/issues).
- Use `dotmatch citation` for the installed release.
- See [Methods and citation](methods-and-citation.md) for `CITATION.cff`, DOI,
  and methods text.
- See [Packaging](packaging.md) for PyPI, Bioconda, containers, and source
  builds.

```{toctree}
:maxdepth: 2
:caption: Getting started

getting-started
command-reference
tutorials/crispr-count-first-run
tutorials/scverse-perturb-seq
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
schemas
workbench
bioinformatics-evaluation
external-review-packet
ecosystem-status
```

```{toctree}
:maxdepth: 2
:caption: Reference

trust-and-scope
benchmarks/README
methods-and-citation
packaging
```
