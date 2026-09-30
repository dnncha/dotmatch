# Genome-wide HAP1 cohort, 2026-09-30

CRISPRWorks Fit at commit `88dacfa271fe2be63f72d68db3b393f638abd86a`
analyzed a public genome-wide HAP1 cohort of **17,445 gene labels and 69,780
guides**. Three paired runs per backend used identical input counts, design,
settings and random seed, alternating execution order.

| Backend | Median wall time | Individual wall times |
| --- | ---: | --- |
| MAGeCK2 0.3.0 reference | 199.740 s | 199.740, 198.434, 200.622 s |
| CRISPRWorks Fit accelerated | 75.633 s | 76.567, 74.727, 75.633 s |

**Measured speedup: 2.64×.** Wall time includes Python startup, the count-table
workflow and full-precision diagnostic output. Printed gene summaries were
byte-identical across all six runs. Permutation p-values and FDR, including
both tails, matched exactly.

| Full-precision field | Maximum absolute difference |
| --- | ---: |
| Beta estimates, including guide baselines | 1.25e-14 |
| Wald z-scores | 3.41e-14 |
| Guide efficiency estimates | 1.45e-15 |
| Wald p-values | 1.97e-14 |
| Wald FDR | 2.11e-14 |
| Permutation p-values and FDR | 0 |

## Data and selection

The source is Hart Lab's BAGEL `reads_hap1.txt` at commit
`53388adbb4fb0931e5c9dda135502be19e4555f0`. Its complete table contains 71,090
guides and 18,056 gene labels, including controls. Source SHA-256:
`7638bc6237cf8a3e2302fcd969d014b041819669b37f2f87bc1f7cfbc45ca8a7`.

This measurement selects labels with exactly four guides, excluding incomplete
labels and larger control bins **before normalization**. Cohort SHA-256:
`bf24bf255b75cde8cb01279290bdbbbbbba5773596245b8c0300b07b879c36bf`.
It is a genome-wide cohort measurement, not a timing for the entire unfiltered
table. The full table requires separate permutation nulls for seven guide-count
groups; it remains an additional stress test.

The design uses HAP1 T0 as one baseline and three T18 replicate samples as one
treatment effect. Runs enable efficiency updates, model mean-variance on 1,000
genes, use two permutation rounds and seed 42, and use one worker and one BLAS
thread per process. Two rounds are the upstream default; the upstream suggested
ten rounds cost more and are not measured here.

Environment: shared GitHub Actions `ubuntu-latest` Linux x86_64 runner, Python
3.12.14, NumPy 2.3.5, SciPy 1.17.0 and OpenBLAS 0.3.30. CPU resources were not
reserved. Peak memory was not measured. These results establish numerical
agreement and a runtime improvement for this workload; they do not establish
biological hit accuracy, independent lab validation or a universal speedup.

## Reproduce and inspect

```bash
python packages/crisprworks-fit/benchmarks/public_hap1.py \
  --out-dir hap1-benchmark --cohort four-guide --repeats 3
```

The harness verifies the source hash, records the selection and analyzes the
same selected matrix with both backends. Full-precision comparisons use
`rtol=1e-7, atol=1e-8`. The [committed record](hap1-four-guide.json) retains
commands, run times, manifests, hashes, provenance and differences.
The [successful CI run](https://github.com/dnncha/dotmatch/actions/runs/36713005248)
also retains the six complete result sets in its public HAP1 artifact.

## Additional full-table numerical comparison

The **unfiltered 71,090-guide, 18,056-label table** was also compared at full
precision. The complete MAGeCK2 reference output came from the CI stress run;
the complete accelerated output came from a local Linux run using the same
input hashes, options and seed. Printed gene summaries were byte-identical,
permutation p-values/FDR matched exactly, and maximum absolute differences were
`7.50e-15` for beta and `6.17e-14` for Wald z-scores.

This is a cross-environment numerical comparison, **not a timing comparison**.
The original 30-minute CI job finished its reference run but timed out near the
end of its first accelerated run. It did not finish three paired repetitions.
The [full-table parity record](hap1-full-parity.json) includes both completed
run manifests, settings, hashes and numerical differences. No full-table
speedup is claimed. One 96-guide label retains upstream's fitting skip and
permutation fallback behavior; all output labels are included in the comparison.

For a paired rerun, use `public_hap1.py --cohort full --repeats 1` with the
output directory option. The manual CI workflow supports the full-table case
and has a larger time budget. Three full-table timing repetitions and
independent lab evaluation remain further checks.
