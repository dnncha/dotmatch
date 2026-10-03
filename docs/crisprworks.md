# CRISPRWorks

**Open-source software for pooled CRISPR screen analysis.**

CRISPRWorks connects guide counting, gene-effect inference and review through
ordinary count tables and explicit sample metadata. **CRISPRWorks Count** is
powered by DotMatch. **CRISPRWorks Fit** is a separate experimental source
package for gene-level inference.

## Components and status

| Name | Task | Current status |
| --- | --- | --- |
| **CRISPRWorks Count** | Turn FASTQ reads and a known guide library into counts and assignment QC | Available through DotMatch |
| **CRISPRWorks Fit** | Fit MAGeCK2 MLE effects, learn guide efficacy jointly across screens, and calibrate effects against control genes | [Experimental alpha available from source](https://github.com/dnncha/dotmatch/tree/main/packages/crisprworks-fit) |
| **CRISPRWorks Review** | Inspect QC, comparisons and screen results | Planned as a broader tool; DotMatch already provides count and assignment-QC reports |

The Count documentation describes released DotMatch functionality. Fit is an
experimental source package with a native MLE accelerator, a JACKS-derived
joint model, whole-gene control calibration and reproducible offline demos.
It has no PyPI release or independent lab evaluation yet. Read
[From counts to gene effects](tutorials/crispr-fit-first-run.md) for input
formats, model choice, outputs and interpretation.
The broader Review tool is planned. The desktop Workbench is developed separately in
`dotmatch-community`; its presence does not imply a released Review product.

## Use Count today

Install the published package and use its existing command:

```bash
python3 -m pip install dotmatch==0.7.0
dotmatch --version
dotmatch demo --out-dir first-run/
```

Continue with the [CRISPR guide-counting tutorial](tutorials/crispr-count-first-run.md).
Count writes MAGeCK-compatible counts for downstream analysis. It does not
perform gene-level inference.

The DotMatch engine also supports barcode demultiplexing and other known-target
short-DNA assays. Those capabilities remain documented under DotMatch.

:::{admonition} Count benchmark spotlight: faster one-mismatch workflows
:class: tip

The benchmarked source checkout's `dotmatch guide-counter count` is
**1.6–5.5× faster with 11.3–14.4× lower peak memory** than guide-counter 0.1.3,
with identical full count matrices. The comparison uses controlled 100k/1M-read
FASTQs, the public Yusa library, one CPU thread and five paired repeats.
Guide-counter is faster in the tested exact-mode million-read cases.
These measurements concern guide counting; gene-effect accuracy and full
experimental-screen performance require separate validation.
[See both modes, graphs and reproduction commands](benchmarks/guide_counter/README.md).
:::

## Try Fit from source

From a repository checkout:

```bash
python3 -m pip install ./packages/crisprworks-fit
crisprworks-fit demo --out-dir fit-demo/
crisprworks-fit demo --model joint --out-dir fit-joint-demo/
```

Fit currently requires Python 3.10+, a C++ compiler and `make` for its MAGeCK2
dependency. The [Fit README](https://github.com/dnncha/dotmatch/tree/main/packages/crisprworks-fit)
also covers the tested container build, count-table analysis and batch workflow.
The default demo runs MLE. The joint demo creates counts, a sample map and 32
declared synthetic control genes, then fits two conditions and calibrates them
with the default global BY adjustment. Both demos exercise software behavior;
they do not validate biological hits. Use a new directory for the joint demo.

## Choose a Fit workflow

| Task | Command | Input and interpretation |
| --- | --- | --- |
| Fit an explicit experimental design | `crisprworks-fit mle` | MAGeCK-compatible counts and a labeled design matrix; natural-log beta coefficients |
| Share guide-efficacy estimates across screens | `crisprworks-fit joint` | The same guide library, counts and a `Sample/Condition/Control` map; log2 effects and relative guide efficacies |
| Calibrate frozen fitted effects | `crisprworks-fit calibrate` | Full-precision MLE or joint JSON plus independently justified control genes; empirical p-values and BH/BY adjusted values |

Count's FASTQ sample sheet and Fit's sample metadata serve different purposes.
Count records which reads belong to a sample. Fit records which samples are
baselines and biological replicates. Preserve sample and guide IDs at the
handoff; do not treat sequencing lanes as independent biological replicates.
Joint inference requires at least two replicates per treatment condition and
64 guides for variance smoothing. The calibrator defaults to 20 control genes
per exact guide-count stratum, global adjustment and BY. It can return no hits
when control counts limit p-value resolution.

## What has been evaluated

The [native MLE HAP1 measurement](https://github.com/dnncha/dotmatch/blob/main/packages/crisprworks-fit/benchmarks/NATIVE.md)
covered 17,445 complete four-guide gene labels. Three paired runs per engine
measured 332.45 seconds for MAGeCK2 reference, 131.06 seconds for the NumPy
accelerator and 38.68 seconds for the native accelerator: 8.59× versus reference
on that cohort and runner. Printed summaries were identical. This is numerical
and runtime evidence, not an improvement in biological hit accuracy.

The [scientific comparison](https://github.com/dnncha/dotmatch/blob/main/packages/crisprworks-fit/benchmarks/SCIENTIFIC-VALIDATION.md)
includes BAGEL2 on two HAP1 screens and current Chronos/JACKS on three cell
lines in two libraries. Joint improved reference-essential ranking over Fit
MLE in all six library/cell comparisons and exceeded Chronos in four. The
panel is a curated, control-enriched subset; several paired intervals overlap
zero. Those results do not establish overall state of the art.

Whole-gene calibration retains correlated guide behavior that pooled-guide
resampling can lose. Its known-truth score simulation reports both false
discoveries and power. Finite tails and multiple-testing adjustment require a
representative null; they do not establish biological FDR control. Copy-number
correction, complete independent screens and biological perturbation evidence
remain outside the current validation.

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

The MLE accelerator preserves a pinned MAGeCK2 workflow; the joint engine and
gene-control calibrator have separate model and inference checks. Before a
broader production release or expanded scientific claim:

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
