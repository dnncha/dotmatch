# CRISPRWorks

**Open-source software for pooled CRISPR screen analysis.**

CRISPRWorks brings related tools under one task-focused name. Its first
available component is **CRISPRWorks Count**, powered by the existing DotMatch
guide-counting engine. The family is intended to support counting, statistical
inference and review as separate, interoperable steps.

## Components and status

| Name | Task | Current status |
| --- | --- | --- |
| **CRISPRWorks Count** | Turn FASTQ reads and a known guide library into counts and assignment QC | Available through DotMatch |
| **CRISPRWorks Fit** | Estimate gene-level effects from guide counts and an experimental design | [Experimental alpha available from source](https://github.com/dnncha/dotmatch/tree/main/packages/crisprworks-fit) |
| **CRISPRWorks Review** | Inspect QC, comparisons and screen results | Planned as a broader tool; DotMatch already provides count and assignment-QC reports |

The Count documentation describes released DotMatch functionality. Fit is an
experimental source package with numerical regression checks, a tested
container build and a batch runner. Its [genome-wide HAP1 cohort measurement](https://github.com/dnncha/dotmatch/blob/main/packages/crisprworks-fit/benchmarks/HAP1.md)
covered 17,445 complete four-guide gene labels and measured 2.64× faster
end-to-end fitting with identical printed gene summaries. The cohort selection
and settings are recorded; this is not a universal speed or biological accuracy
claim. Fit has no PyPI release or independent lab evaluation yet.
The broader Review tool is planned. The desktop Workbench is developed separately in
`dotmatch-community`; its presence does not imply a released Review product.

## Use Count today

Install the published package and use its existing command:

```bash
python3 -m pip install dotmatch==0.6.4
dotmatch --version
dotmatch demo --out-dir first-run/
```

Continue with the [CRISPR guide-counting tutorial](tutorials/crispr-count-first-run.md).
Count writes MAGeCK-compatible counts for downstream analysis. It does not
perform gene-level inference.

The DotMatch engine also supports barcode demultiplexing and other known-target
short-DNA assays. Those capabilities remain documented under DotMatch.

## Try Fit from source

From a repository checkout:

```bash
python3 -m pip install ./packages/crisprworks-fit
crisprworks-fit demo --out-dir fit-demo/
```

Fit currently requires Python 3.10+, a C++ compiler and `make` for its MAGeCK2
dependency. The [Fit README](https://github.com/dnncha/dotmatch/tree/main/packages/crisprworks-fit)
also covers the tested container build, count-table analysis and batch workflow.
The synthetic demo checks installation; it does not validate biological hits.

## Names and compatibility

Use **CRISPRWorks** for the family, and **CRISPRWorks Count**, **CRISPRWorks Fit**
and **CRISPRWorks Review** for its components. Spell the family name as one word,
with uppercase `CRISPR` and an uppercase `W`.

For the existing engine, use **“CRISPRWorks Count, powered by DotMatch”** on
introductions and product pages. Use **DotMatch** for package names, commands,
Python imports, methods and citations.

The branding change preserves the existing identifiers:

| Identifier | Value |
| --- | --- |
| Python package and import | `dotmatch` |
| Command | `dotmatch` |
| GitHub repository | `dnncha/dotmatch` |
| Container | `ghcr.io/dnncha/dotmatch` |
| Documentation host | `dotmatch.readthedocs.io` |
| Release citation | Existing DotMatch release metadata and DOI |

Users can keep their current installations and workflows. Record the actual
DotMatch version and assignment policy in methods; the family name does not
replace the versioned software citation.

## Fit development criteria

The first Fit implementation targets the MAGeCK2 MLE inference step. Before a
broader production release or expanded performance claim:

1. Profile a pinned upstream version on reproducible workloads.
2. Compare inference using identical frozen count matrices, design matrices,
   settings and permutation inputs. Record numerical tolerances and differences
   in effect estimates, standard errors, p-values and adjusted p-values.
3. Measure runtime on a documented environment. Measure peak memory before
   making memory claims; report workload-specific limits.
4. Obtain independent lab workflow evaluation and investigate discrepancies
   before claiming production readiness.

Evaluate guide-counting changes separately from inference changes. DotMatch
count comparisons can help establish the former; they do not validate
gene-level statistics. Upstream collaboration and compatibility with existing
pipelines are development priorities.

## Rollout

Introduce the family through the README and documentation first. Retain the
published package names, links and citations while the additional components
are evaluated. Any future command aliases or distribution changes need their
own compatibility plan.

CRISPRWorks is the selected family name. Domain, package-registry and trademark
availability have not been established; this document does not claim ownership
of those names or namespaces.
