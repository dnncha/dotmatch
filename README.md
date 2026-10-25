# CRISPRWorks Count

Count CRISPR guides from FASTQ files and review how reads were assigned.

**CRISPRWorks Count, powered by DotMatch**, turns a known guide library and your
FASTQ files into MAGeCK-compatible count tables and assignment quality checks
(QC). It runs locally on Linux and macOS. Install the `dotmatch` package and use
the `dotmatch` command.

The DotMatch engine also supports fixed-position barcode demultiplexing and other
short-DNA assays with known targets.

[![PyPI](https://img.shields.io/pypi/v/dotmatch?label=PyPI)](https://pypi.org/project/dotmatch/)
[![Documentation](https://readthedocs.org/projects/dotmatch/badge/?version=latest)](https://dotmatch.readthedocs.io/en/latest/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20541628.svg)](https://doi.org/10.5281/zenodo.20541628)

[Count guides from your own reads](https://dnncha.github.io/dotmatch/crispr-guide-counting/) ·
[Check library collisions in your browser](https://dnncha.github.io/dotmatch/tools/library-safety/) ·
[Documentation](https://dotmatch.readthedocs.io/en/latest/) ·
[Methods and results](https://dotmatch.readthedocs.io/en/latest/benchmarks/crispr_comparison/README.html)

## Install and run the example

Install the published DotMatch 0.7.0 release from PyPI:

```bash
python3 -m pip install dotmatch==0.7.0
dotmatch --version
dotmatch demo --out-dir first-run/
```

Open `first-run/index.html` to review the counts, read outcomes, and reports.
The example includes synthetic reads and checks its results against expected
counts. It does not use study data or make network requests. Use a new output
directory each time.

To explore before installing, open the
[interactive synthetic example](https://dnncha.github.io/dotmatch/assignment-sensitivity/).
To compare your own raw-count tables, use
[`dotmatch compare-counts`](https://dotmatch.readthedocs.io/en/latest/count-comparison.html).

### Install with Conda or a container

```bash
conda create -n dotmatch -c conda-forge -c bioconda dotmatch
conda activate dotmatch

# Or use the pinned release container:
docker run --rm ghcr.io/dnncha/dotmatch:v0.7.0 --version
```

Release 0.7.0 is also published on GHCR. Bioconda and its generated BioContainers
images can lag PyPI and GHCR. When a
newly tagged version has not reached Bioconda yet, use PyPI or the source build.
Check the installed version. The Bioconda recipe includes `osx-arm64` for Apple Silicon.
Review the [packaging details](https://dotmatch.readthedocs.io/en/latest/packaging.html)
for platform and container verification. See the [installation guide](https://dotmatch.readthedocs.io/en/latest/getting-started.html)
for platform details, source builds and the third-party Homebrew tap.

## Understand the read outcomes

DotMatch reports what happened to each read under your chosen matching rule:

| Outcome | Meaning |
| --- | --- |
| **unique** | The rule assigns the read to one target. Adds one count to that target. |
| **ambiguous** | The rule cannot choose one target. Adds no count. |
| **none** (unmatched) | No target meets the rule. Adds no count. |
| **invalid** | The requested sequence window could not be extracted. Adds no count. |

Only unique assignments contribute to target counts. A unique assignment does
not prove biological origin. Review per-guide counts and ambiguous reads as well
as the overall assigned-read percentage.

Use the counts in your downstream analysis workflow. DotMatch does not perform
genome alignment, basecalling, cell calling, UMI deduplication, or gene-level
hit calling.

## Count a CRISPR screen

You need a guide CSV or TSV with target sequences and the FASTQ files from your
screen. The [first-run tutorial](https://dotmatch.readthedocs.io/en/latest/tutorials/crispr-count-first-run.html)
explains the library columns and provides a small example. Keep biological
sample names distinct from sequencing filenames when preparing the library and
sample configuration. Select the read that contains the guide sequence.
For paired-end data, do not include both R1 and R2 by default.

Start a new assay project:

```bash
dotmatch crispr quickstart \
  --library guides.csv \
  --fastq 'fastqs/*_R1*.fastq.gz' \
  --out crispr_screen/
```

Adjust the file pattern to match your guide reads. The command creates a draft
project. Review `crispr_screen/inference_report.json`, `samples.generated.tsv`,
and `assay.toml`. Confirm the guide window, orientation, matching rule, library,
and sample files.
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

## Check how matching rules change counts

`dotmatch sensitivity` compares exact matching with two rules that allow
substitutions: radius-one and best-distance Hamming assignment. It uses the same
read windows in **one FASTQ pass**.
It produces three count matrices, per-guide deltas, read-state transitions,
checksums and a self-contained HTML report. It never selects a policy for you.

The following command uses example files from a source checkout of tag `v0.7.0`.
Run it from the repository root; installing the package alone does not provide
these paths.

```bash
python3 -m pip install dotmatch==0.7.0
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

[Explore the interactive review example](https://dnncha.github.io/dotmatch/assignment-sensitivity/)
without installing anything. The example uses the public nine-read synthetic fixture.
The installed `dotmatch sensitivity` command continues to write its static
report. Both examples demonstrate software behavior, not biological accuracy.

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

## The CRISPRWorks family

CRISPRWorks is a family of open-source tools for pooled CRISPR screen analysis.
Count is available through DotMatch. The other components have separate release
status:

| Component | Purpose | Availability |
| --- | --- | --- |
| **Count** | Count guides and inspect read assignments | Available through the DotMatch command and Python package |
| **Fit** | Estimate gene effects, learn guide efficacy across screens, and calibrate results against control genes | [Experimental alpha available from source](https://github.com/dnncha/dotmatch/tree/main/packages/crisprworks-fit) |
| **Review** | Review QC, comparisons, and screen results | Planned; DotMatch already provides count and assignment-QC reports |

Read the [family overview](https://dotmatch.readthedocs.io/en/latest/crisprworks.html)
for component status and scientific scope. For gene-level analysis, follow
[From counts to gene effects](https://dotmatch.readthedocs.io/en/latest/tutorials/crispr-fit-first-run.html).
Fit requires a separate source install and has no published package release.
The desktop Workbench is developed separately as `dotmatch-community`.

## Research and reproducibility

[Cheerful Duck Research](https://cheerfulduck.com/research) publishes our
bioinformatics software investigations, reproducible examples, and follow-up
corrections. For DotMatch-specific methods and measurements, use the reports
below and record the release and assignment policy used in your own run.

### Reproduce a comparison

The benchmarked source checkout's `dotmatch guide-counter count` is
**1.6–5.5× faster with 11.3–14.4× lower peak memory** than unmodified
guide-counter 0.1.3 in the tested one-mismatch workflows. Every paired run
produced identical full guide-by-sample count matrices.

The comparison uses controlled 100k/1M-read FASTQs against the public
87,437-guide Yusa library, one CPU thread, and five paired repeats. Guide-counter
is faster in the tested exact-mode million-read cases. Full experimental-screen
performance and biological accuracy require separate validation.

The [benchmark reports](https://dotmatch.readthedocs.io/en/latest/benchmarks/README.html)
include commands, hardware and assignment rules. Those reports cover the tested workloads;
they are not universal speed or biological-accuracy guarantees.

The [guide-counter comparison](https://dotmatch.readthedocs.io/en/latest/benchmarks/guide_counter/README.html)
records complete-command throughput and peak memory with identical full counts.
The figures below show exact and one-mismatch modes, with median bars and the
observed range across five runs. They cover the benchmarked source checkout
and controlled 100k/1M-read FASTQs against the Yusa library.

![DotMatch versus guide-counter: complete-command throughput in exact and one-mismatch modes](https://raw.githubusercontent.com/dnncha/dotmatch/main/benchmarks/figures/guide_counter_throughput.svg)

![DotMatch versus guide-counter: peak memory for the same complete commands](https://raw.githubusercontent.com/dnncha/dotmatch/main/benchmarks/figures/guide_counter_memory.svg)

The current source engine's [Hamming improvement report](https://dotmatch.readthedocs.io/en/latest/benchmarks/hamming_seed_layout/README.html)
records 5.2–5.7× faster complete `k=3` counting against the previous engine on
30,000 simulated reads and the real Yusa guide library, with identical counts.
It also documents a corrected false-tie bug in Hamming queries containing
literal unknown bytes and 136,096 independent oracle checks.

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

Release 0.7.0 includes six structured agent tools. Inspect the available tools
and export instructions for your local agent:

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
