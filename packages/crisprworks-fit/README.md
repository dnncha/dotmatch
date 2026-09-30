# CRISPRWorks Fit

**Experimental acceleration of MAGeCK2 MLE for pooled CRISPR screens.**

Fit runs MAGeCK2's count-table-to-gene-results workflow with faster numerical
kernels and worker reuse. It belongs to the CRISPRWorks family alongside
CRISPRWorks Count, powered by DotMatch. This is a separate installable package;
the existing DotMatch package remains independent of its dependencies.

The prototype targets **MAGeCK2 0.3.0**. It checks the upstream functions against
the tested source before enabling acceleration. It is not a drop-in replacement
for every version of legacy MAGeCK.

## Install and try

From a checkout containing this directory:

```bash
python3 -m pip install ./packages/crisprworks-fit
crisprworks-fit --version
crisprworks-fit demo --out-dir fit-demo/
```

MAGeCK2 is an installation dependency and currently builds its bundled C++
helpers, so its installation requires a C++ compiler and `make`. Python 3.10+
is required by Fit. The package has not been published to PyPI.

The demo creates a small synthetic count table and design, then writes gene
and guide summaries. It demonstrates software behavior, not biological accuracy.

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
| `screen.fit-details.json` | Optional full-precision estimates, efficiencies, p-values and FDR for numerical comparisons |
| `screen.log` | MAGeCK2 workflow log |

Use `--backend reference` to run the same workflow with the original MAGeCK2
functions. Use identical input files, options and seed when comparing backends.
Both modes default to one BLAS thread per process. `--threads` controls worker
processes; `--blas-threads` controls BLAS threads within each worker.

## What is faster

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
- Sort permutation-null columns once and use binary searches. Strict greater-
  than/less-than comparisons preserve the reference treatment of ties.

The dispersion model, normalization, convergence threshold, efficiency rules,
permutation grouping, gene skip threshold and FDR calculation retain upstream
conventions. The adapter restores upstream functions and the design cache when
it exits. Its scoped patches are process-global; concurrent unwrapped MAGeCK2
calls in the same interpreter are unsupported.

## Validation and measurements

The regression suite checks 54 unrounded fits across three designs, three guide
counts, two random seeds, and disabled/fixed/updated efficiency modes. It also
checks dense-versus-compact covariance algebra, restarts, strict permutation
tails including ties and NaNs, FDR, source compatibility, failure restoration,
invalid designs, the public upstream count-table fixture and spawned workers.

The [genome-wide HAP1 cohort measurement](benchmarks/HAP1.md) covered 17,445
complete four-guide gene labels: median wall time was **199.74 s reference vs
75.63 s accelerated (2.64×)** across three runs per backend. Printed gene
summaries were byte-identical, permutation p-values/FDR matched exactly, and
full-precision beta differences were below `1.25e-14`. This selects complete
four-guide labels before normalization; it is not a timing for the unfiltered
table or evidence of biological hit accuracy.

The unfiltered 71,090-guide table also passed a separate cross-environment
full-precision comparison with identical printed gene results; see the HAP1
record for how the CI reference and local accelerated outputs were compared.
No unfiltered-table timing claim is made.

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
behavior. Fit is currently an optional Python/NumPy accelerator, with native
linear algebra provided by BLAS. It does not implement counting or change
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
