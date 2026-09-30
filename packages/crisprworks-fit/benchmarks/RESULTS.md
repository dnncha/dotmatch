# Initial measurements, 2026-09-30

For the newer measured genome-wide cohort, see [HAP1.md](HAP1.md): 17,445
complete four-guide gene labels, three paired repetitions, 2.64× end-to-end
speedup, identical printed gene summaries and full-precision parity.

These measurements compare CRISPRWorks Fit's accelerated and reference backends
at commit `e900ecdb903fe2eb25087ca7544a93412c62f49e`, using MAGeCK2 0.3.0.
Input counts, design, settings and permutation seed are held
fixed. They establish a working prototype on two bounded workloads, not
biological accuracy or a universal speed guarantee.

| Workload | Reference median | Accelerated median | Speedup | Printed gene results |
| --- | ---: | ---: | ---: | --- |
| Synthetic: 1,000 genes, 6 guides/gene, 6 samples | 7.972 s | 3.532 s | 2.26× | Identical bytes |
| Upstream fixture: 100 genes, 999 guides, 4 samples | 3.619 s | 1.816 s | 1.99× | Identical bytes |

Each median uses three runs per backend, alternating execution order. Times
include Python startup and full-precision diagnostic output. Each run uses one
worker, one BLAS thread, two permutation rounds and seed 42. The synthetic run
uses the upstream default constant dispersion and fixed efficiencies; the
fixture run enables efficiency updates and mean-variance modeling on 20 genes.

Environment: shared Linux x86_64 virtual runner, AMD EPYC 9V74, Python 3.12.14,
NumPy 2.3.5, SciPy 1.17.0 and OpenBLAS 0.3.30. CPU time was not reserved; no
peak-memory measurement or industry-scale throughput claim is made.

Full-precision comparisons passed `rtol=1e-7, atol=1e-8` across all runs:

| Maximum absolute difference | Synthetic | Upstream fixture |
| --- | ---: | ---: |
| Beta estimates, including guide baselines | 1.16e-14 | 9.33e-15 |
| Wald z-scores | 8.42e-14 | 3.29e-14 |
| Guide efficiency estimates | 0 | 5.56e-15 |
| Wald p-values | 4.91e-14 | 2.34e-14 |
| Wald FDR | 5.42e-14 | 2.69e-14 |
| Permutation p-values and FDR, both tails | 0 | 0 |

The separate regression suite checks standard-error parity over 54 fitting
cases. Compact-versus-dense covariance calculations use a tighter algebraic
tolerance of `rtol=1e-10, atol=1e-12`.

## Reproduce

Install the package from the source checkout, then run from the repository root:

```bash
python packages/crisprworks-fit/benchmarks/run.py \
  --out-dir synthetic-benchmark --genes 1000 --repeats 3

python packages/crisprworks-fit/benchmarks/run.py \
  --out-dir fixture-benchmark \
  --count-table packages/crisprworks-fit/tests/data/mageck2-counts.tsv \
  --design-matrix packages/crisprworks-fit/tests/data/design.tsv \
  --update-efficiency --genes-varmodeling 20 --repeats 3
```

The generated `benchmark.json` contains commands, input hashes, individual run
times, software/BLAS configuration and numerical differences. The committed
[synthetic record](synthetic-1000.json) and [fixture record](upstream-fixture.json)
retain the measurements, settings and hashes from the runs summarized here.
Fixture provenance is recorded in [tests/data/README.md](../tests/data/README.md).

Genome-wide public screens, larger designs, numerical boundary cases and
independent workflow evaluations remain release gates.

The newer `public_hap1.py` harness performs a pinned full-library HAP1
comparison in CI. Results from that run must be assessed separately from the
initial measurements above; the required-likelihood optimization was added
after these initial records.
