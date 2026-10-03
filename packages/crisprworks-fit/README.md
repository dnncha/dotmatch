# CRISPRWorks Fit

**Experimental gene-level inference for pooled CRISPR screens.**

Fit runs MAGeCK2's count-table-to-gene-results workflow with faster numerical
kernels and worker reuse. A separate JACKS-derived engine learns guide efficacy
jointly across replicated screens, and a gene-control calibrator preserves
within-gene guide correlation. Alpha 2 adds a compiled C EM loop; the earlier
NumPy engine remains selectable for reproducible comparisons. It belongs to the CRISPRWorks family alongside
CRISPRWorks Count, powered by DotMatch. This is a separate installable package;
the existing DotMatch package remains independent of its dependencies.

The MLE adapter targets **MAGeCK2 0.3.0**. It checks the upstream functions against
the tested source before enabling acceleration. It is not a drop-in replacement
for every version of legacy MAGeCK.

## Install and try

From a checkout containing this directory:

```bash
python3 -m pip install ./packages/crisprworks-fit
crisprworks-fit --version
crisprworks-fit demo --out-dir fit-demo/
crisprworks-fit demo --model joint --out-dir fit-joint-demo/
```

MAGeCK2 is an installation dependency and currently builds its bundled C++
helpers, so its installation requires a C++ compiler and `make`. Python 3.10+
is required by Fit. Building the optional native loop requires a C compiler.
`--kernel auto` uses the native loop when available and otherwise NumPy;
`--kernel native` fails explicitly if the compiled engine is unavailable;
`--kernel numpy` selects the previous engine. The provenance manifest records
the engine actually used. The package has not been published to PyPI.

The default demo creates a small synthetic count table and design, then writes
MLE gene and guide summaries. The joint demo creates 48 genes, four guides per
gene, a sample map and 32 independently declared synthetic control genes. It
fits two conditions and writes whole-gene calibration with default global BY
adjustment. A `demo.json` completion marker records input/output hashes. The
joint demo requires a new output directory and uses fixed seed 1729. Both
demos demonstrate software behavior, not biological accuracy.

The [Count-to-Fit tutorial](https://dotmatch.readthedocs.io/en/latest/tutorials/crispr-fit-first-run.html)
covers model choice, file formats, named conditions and result interpretation.
`--backend`, `--kernel` and `--threads` select MLE demo behavior; joint inference
uses its own NumPy model.

For an isolated installation with pinned numerical dependencies:

```bash
docker build -t crisprworks-fit:alpha packages/crisprworks-fit
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD:/data" crisprworks-fit:alpha demo --out-dir fit-demo --threads 2
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD:/data" crisprworks-fit:alpha mle \
  -k counts.tsv -d design.tsv -n results/screen --seed 42
```

This builds locally; no published container image is implied. The container
uses Python 3.12, MAGeCK2 0.3.0, NumPy 2.3.5 and SciPy 1.17.0. CI builds it
and checks an installed two-worker demo. The base image is a version tag rather
than an immutable digest, so rebuilds are not byte-for-byte reproducible.

## Analyze counts

Use a MAGeCK-format count table and a labeled design matrix:

```bash
crisprworks-fit mle \
  -k counts.tsv -d design.tsv -n results/screen \
  --threads 1 --seed 42 --write-fit-details
```

Count tables have `sgRNA`, `Gene`, then sample columns. The first design column
is all ones; the first row must be a baseline sample with zeros in every
condition column. Sample labels select the corresponding count-table columns.
DotMatch's MAGeCK-compatible counts can be used as input.

For several screens, the [batch runner](examples/README.md) accepts a JSON
screen list, checks inputs and output collisions before running, and retains a
separate provenance manifest for each screen.

Outputs:

| File | Contents |
| --- | --- |
| `screen.gene_summary.txt` | Upstream-compatible beta, z-score, permutation p-value/FDR and Wald p-value/FDR columns |
| `screen.sgrna_summary.txt` | Upstream-compatible guide results |
| `screen.crisprworks.json` | Backend, seed, options, input/output hashes, versions, BLAS configuration and elapsed time |
| `screen.fit-details.json` | Optional full-precision estimates, efficiencies, p-values/FDR, condition labels and effect units for comparison or calibration |
| `screen.log` | MAGeCK2 workflow log |

Use `--backend reference` to run the same workflow with the original MAGeCK2
functions. Use identical input files, options and seed when comparing backends.
Both modes default to one BLAS thread per process. `--threads` controls worker
processes; `--blas-threads` controls BLAS threads within each worker.

### Finite permutation inference with negative controls

For an analysis with independently justified nonessential control genes:

```bash
crisprworks-fit mle -k counts.tsv -d design.tsv -n results/screen \
  --control-gene nonessential-controls.txt --norm-method control \
  --permutation-pvalues finite --update-efficiency --seed 42 \
  --permutation-round 10 --write-fit-details
```

The control file contains one gene name per line, without a header. Choose
controls independently of the effects you are testing. This uses their guides
for normalization and the grouped permutation background. A matched null can
matter much more for hit sensitivity than a high ranking score; the
[held-out HAP1 benchmark](benchmarks/SCIENTIFIC-VALIDATION.md) measures both.
Controls are assay-specific: an annotation of nonessentiality elsewhere does
not prove that a gene is neutral in your treatment or cell line.

`--permutation-pvalues finite` counts ties and uses `(extreme + 1)/(draws + 1)`
in each directional tail. Two-sided p-values double the smaller tail, capped
at one. An observation outside a finite simulated null never receives `p=0`;
an all-tied null receives `p=1`. Nonfinite observations or null columns and
genes skipped by the upstream guide-count threshold receive `p=1`. The
manifest records the null source, draw counts, p-value resolution and failed
fits. More draws improve numerical resolution; they do not repair an invalid
null model. Finite control inference rejects `--no-permutation-by-group`
because that upstream branch ignores the control-guide background.

The default `--permutation-pvalues legacy` retains MAGeCK compatibility,
including strict tails and zero p-values. Finite mode changes permutation
p-values and their FDR, while preserving beta estimates, guide efficiencies
and Wald statistics. It also works with `--backend reference`, where only the
fitting is the original implementation. The add-one safeguard follows
[Phipson and Smyth (2010)](https://doi.org/10.2202/1544-6115.1585);
valid randomization inference still requires exchangeability. Pooled-guide
pseudo-genes do not establish that assumption or biological FDR control.

### Learn guide efficacy jointly across screens

For screens using the same guide library, provide a tab-separated sample map:

```text
Sample	Condition	Control
baseline	BASE	BASE
cell_a_1	CELL_A	BASE
cell_a_2	CELL_A	BASE
cell_b_1	CELL_B	BASE
cell_b_2	CELL_B	BASE
```

Every count-table sample must appear once. Baseline conditions map to
themselves; treatment conditions name their matched baseline and require at
least two replicates. Baselines may have one replicate. Different treatment
conditions can use different baseline groups. The table needs at least 64
guides for variance smoothing; counts must be finite and nonnegative.

```bash
crisprworks-fit joint -k counts.tsv --sample-map samples.tsv \
  --control-gene nonessential-controls.txt -n results/joint
```

This writes `joint.joint.tsv` and full-precision `joint.joint-details.json`.
Effects use log2 relative-abundance units. Guide efficacies are continuous
relative weights shared across conditions, rather than editing probabilities.
Posterior standard deviations and z scores describe the model; they are not
calibrated biological hit probabilities. The model adapts
[JACKS](https://doi.org/10.1101/gr.238923.118), with pseudocount 32 and its
default variational updates. No pretrained efficacy reference, hierarchical
effect prior or copy-number correction is enabled. Omitting the control file
uses all-guide median normalization. MLE remains separately selectable.

### Calibrate frozen effects with whole control genes

```bash
crisprworks-fit calibrate --fit-details results/joint.joint-details.json \
  --control-gene nonessential-controls.txt -n results/calibrated
```

The same command accepts MLE's `screen.fit-details.json`. It compares each
effect with entire control-gene effects having the same guide count, retaining
within-gene correlation. Ties use inclusive add-one empirical tails. Training
controls cannot become hits and are excluded from the testing family. Missing
support, nonfinite fits and excluded guide-count strata receive p=1 with an
explicit status. The default requires 20 controls per stratum and excludes
genes with 40 or more guides; those policy limits are configurable.

New MLE details preserve design-coefficient labels, and joint details preserve
sample-map condition labels. Calibration carries these names into its TSV as
`condition_label`, alongside a one-based `condition` index. Older details without
labels retain numeric indices; keep the corresponding design or sample map.

Outputs are `calibrated.calibration.tsv` and a JSON manifest with input/output
hashes, strata, p-value resolution and testing-family size. `--score effect`
uses the fitted effect; `--score z` uses its descriptive posterior/Wald z
score. The default `--family global` adjusts every tested gene/condition
hypothesis together. `--family condition` declares a separate family per
condition. Each alternative has its own family; use two-sided results when
choosing effect direction after seeing the data.

The default `--adjust by` uses Benjamini–Yekutieli adjustment for arbitrary
dependence of valid marginal p-values. It can lose substantial power: it made
no calls on the bounded public panel because control counts limit resolution.
`--adjust bh` offers more power when independence or suitable positive
dependence is justified. Neither adjustment repairs an invalid control null.
Controls must remain exchangeable with tested null genes after the complete
normalization and fitting procedure; selecting controls from observed effects
or using unrepresentative annotations can violate that requirement.

The [correlated-guide simulation](benchmarks/NULL-CALIBRATION.md) demonstrates
why preserving gene correlation matters and quantifies the power tradeoff.

## What is faster

- Execute the EM loop in compiled C, with working storage reused across
  iterations. For finite inputs, process exact nonzero design entries while
  preserving row summation order; nonfinite inputs use the dense calculation.
- Apply IRLS weights by vector multiplication and solve the smaller ridge
  system directly. No dense diagonal weight matrix is constructed.
- Calculate the sandwich covariance after the last iteration, and compute the
  hat-matrix trace terms using the smaller Gram matrix. No observation-by-
  observation hat matrix is constructed.
- Skip discarded diagnostic likelihood evaluations when efficiency updates and
  debug output do not require them.
- Evaluate required negative-binomial log likelihoods directly using SciPy
  special functions, preserving reference arithmetic and distribution support.
- Run single-worker fitting directly; reuse a spawned worker pool across
  fitting stages and permutation rounds for multiple workers.
- Sort permutation-null columns once and use binary searches. Legacy mode's
  strict greater-than/less-than comparisons preserve the reference treatment
  of ties; finite mode uses inclusive comparisons and the add-one correction.

In MLE acceleration, the dispersion model, normalization, convergence threshold,
efficiency rules,
permutation grouping, gene skip threshold and FDR calculation retain upstream
conventions. The adapter restores upstream functions and the design cache when
it exits. Its scoped patches are process-global; concurrent unwrapped MAGeCK2
calls in the same interpreter are unsupported.

## Validation and measurements

The [scientific inference evaluation](benchmarks/SCIENTIFIC-VALIDATION.md)
includes current Chronos and JACKS on a frozen Broad/Sanger panel, plus
BAGEL2 on two public HAP1 complete-four-guide cohorts. The production joint
engine improves AP over Fit MLE in all six panel comparisons (three cell
lines, two libraries) and exceeds Chronos in four. This is a curated vignette
subset enriched for controls; several paired intervals overlap zero. Exact
same-dependency JACKS parity, versions and cross-version differences are
recorded in the [frontier report](benchmarks/FRONTIER.md).

For the earlier HAP1 analysis, five external gene folds exclude evaluation genes from both BAGEL2 training
and Fit's nonessential controls. With finite tails, using control genes for
normalization and the permutation null raised reference-essential recall at
nominal directional FDR 0.05 from **1.6% to 94.5%** and **57.2% to 98.7%**.
Held-out nonessential call rates were **0.7% and 1.2%**. These class-specific
false-positive frequencies do not estimate genome-wide FDR. Control-based
ranking AP was **0.9983 and 0.9984**, comparable to BAGEL2's **0.9980 and
0.9985** on these annotations. The screens share a cell line and data source;
this evidence does not establish superiority across assays or overall SOTA.

The regression suite checks 54 unrounded fits across three designs, three guide
counts, two random seeds, and disabled/fixed/updated efficiency modes. It also
checks dense-versus-compact covariance algebra, restarts, strict permutation
tails including ties and NaNs, FDR, source compatibility, failure restoration,
invalid designs, the public upstream count-table fixture and spawned workers.

The [native-engine HAP1 measurement](benchmarks/NATIVE.md) covered 17,445
complete four-guide gene labels and 69,780 guides. Three paired runs per engine
measured **332.45 s MAGeCK2 reference, 131.06 s NumPy accelerator, and 38.68 s
native accelerator**: **8.59× versus reference and 3.39× versus NumPy**.
All nine printed gene summaries were byte-identical; permutation p-values and
FDR matched exactly. Maximum full-precision beta differences were `2.49e-14`.
These timings include startup and diagnostic output on one shared CI runner.

This cohort selects complete four-guide labels before normalization and
excludes larger control bins. The unfiltered 71,090-guide, 18,056-label table
also matched the reference in a separate cross-environment full-precision
comparison. That comparison establishes numerical agreement, with no
unfiltered-table timing ratio claimed. These checks preserve the reference
results; they do not establish superior biological hit accuracy or a universal
speedup. The native timing experiment compares MLE engines; the separate
frontier report evaluates Chronos/JACKS ranking and calling choices.

The [earlier NumPy measurement](benchmarks/HAP1.md) remains archived with its
original implementation commit and runner timings. Compare engines within each
paired experiment rather than combining times from different environments.

Initial measurements are documented in [benchmarks/RESULTS.md](benchmarks/RESULTS.md).
The initial records cover two bounded workloads at commit `e900ecdb`; they do
not establish a universal speedup or industry readiness. Full-library public
parity and performance are also checked in CI using a pinned HAP1 screen:

```bash
python3 packages/crisprworks-fit/benchmarks/public_hap1.py \
  --out-dir hap1-benchmark --repeats 3
```

The script downloads the Hart Lab BAGEL HAP1 count table at an immutable commit,
checks its SHA-256, and analyzes 71,090 guides across 18,056 gene labels with
efficiency updates and mean-variance modeling. Labels with 40 or more guides
retain upstream's default permutation skip behavior. See the generated
`provenance.json` and `paired/benchmark.json`. This checks numerical agreement
on real counts; it does not establish biological hit accuracy.

Normal push CI uses `--cohort four-guide`: 17,445 complete four-guide gene labels
and 69,780 guides. This excludes incomplete labels and larger control bins
before normalization, and compares the same resulting counts in both backends.
The full table retains all seven guide-count groups and is a separate, more
expensive stress test. The manual CI workflow defaults to that full table.

From the repository root:

```bash
python3 -m unittest discover -s packages/crisprworks-fit/tests -v
python3 packages/crisprworks-fit/benchmarks/run.py \
  --out-dir benchmark/ --genes 1000 --repeats 3
```

The benchmark alternates backend order, checks full-precision results at
`rtol=1e-7, atol=1e-8`, requires identical printed gene summaries, and records all
commands and input hashes. The underlying covariance algebra is tested at
`rtol=1e-10, atol=1e-12`.

## Prototype scope

CNV correction and `--debug-gene` are rejected until separately evaluated.
Experimental Bayes options retain upstream's explicit rejection. The upstream
active fitting path does not apply `--remove-outliers`; Fit preserves that
behavior. Fit offers a compiled C EM engine and a NumPy fallback, with BLAS
used for covariance calculations. The native loop omits exact zero design
entries for finite inputs and retains a dense path for nonfinite inputs. It does not implement counting or change
DotMatch's read-assignment rules.

Keep counting comparisons separate from inference comparisons. A faster fit
with frozen counts does not establish that a different counting policy is
biologically accurate.

## Attribution

The EM control flow and statistical conventions derive from
[MAGeCK2](https://github.com/davidliwei/mageck2), copyright Wei Li, under
BSD-3-Clause. The notice is retained in [LICENSE](LICENSE). The reference release
commit is `630aea0b6fc152a81006435f21d629297275c911`. Cite the upstream methods
and actual software versions used in analyses. CRISPRWorks Fit does not imply
upstream endorsement.

The joint variational updates and posterior variance smoothing are adapted
from [JACKS](https://github.com/felicityallen/JACKS), copyright Felicity Allen
and Leopold Parts, commit `dd5c4be5e83baa5ee79a589b0a8c3a2ac3f7d6ad`.
The source package's MIT notice and repository-root Apache-2.0 license are
bundled under [licenses](licenses/README.txt); adapted source files record their
changes. Cite [JACKS's method](https://doi.org/10.1101/gr.238923.118) for this
model as well as the software versions used.

Joint details include per-gene iteration counts, termination reasons and final
stopping-statistic changes. The CLI exposes `--max-iterations` (default 50) and
`--tolerance` (default 0.1). Meeting this inherited stopping rule is not an
independent posterior-accuracy guarantee. Calibration rejects records explicitly
marked nonconverged and verifies any recorded companion summary hash before
writing outputs. Older standalone files without these records remain supported
with the corresponding verification limitations.

Calibration's `resolution` manifest block reports optimistic attainable q-value
bounds for its actual declared family and eligible strata. `--fdr-alpha` selects
the diagnostic threshold (default 0.05), without changing p/q values. A warning
identifies families where no discovery is attainable; zero calls in that case
must not be interpreted as evidence of no effects. These bounds do not predict
power or validate the negative-control assumptions.

Joint and calibration guard publication with a nonblocking OS lock. Overlapping
publishers at one prefix fail before replacement; distinct prefixes retain both
runs. Hidden `.publish.lock` files remain in place, and the OS releases the lock
on process exit. Input provenance hashes the exact bytes parsed; input changes
or late output/input aliases abort publication before replacing prior results.
Completion manifests and hashes must still be checked. Filesystem crash durability remains separate work.

MLE now runs upstream fitting from private input copies and stages both summary
tables plus optional full-precision details. Failed fitting/staging preserves
an earlier completed bundle. Publication validates input hashes, locks the
prefix, invalidates the old completion manifest, replaces the outputs and
publishes a new complete manifest last. A failure during replacement leaves no
completion manifest; individual table existence does not establish completion.
A successful rerun without `--write-fit-details` removes stale details. Failed
attempts report errors without replacing an earlier completion manifest with a
failure record. Check the command exit status as well as the manifest/hashes.

MLE validates its private count snapshot before fitting: columns must be distinct
and nonempty, rows must have matching widths, guide IDs must be globally unique,
and counts must be finite and nonnegative. It preserves upstream literal TSV/CSV
parsing and filename-suffix delimiter selection. These checks validate input
structure, not the biological assumptions of the fitted model.
