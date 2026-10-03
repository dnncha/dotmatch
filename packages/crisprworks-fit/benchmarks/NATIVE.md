# Native-engine HAP1 benchmark, 2026-09-30

CRISPRWorks Fit alpha 2 at implementation commit
`2938b8e6e2219e80404280facf14f323e55ee0c3` analyzed **17,445 gene labels and
69,780 guides** from the pinned public Hart Lab BAGEL HAP1 table. All three
engines used identical counts, design, options and random seed. Execution order
rotated across three repetitions per engine on the same shared CI runner.

| Engine | Median wall time | Individual wall times |
| --- | ---: | --- |
| Unmodified MAGeCK2 0.3.0 fitting functions | 332.452 s | 336.551, 330.837, 332.452 s |
| Previous Fit NumPy engine | 131.060 s | 131.060, 129.386, 132.549 s |
| Fit native engine | 38.681 s | 38.838, 38.590, 38.681 s |

**8.59× faster than reference; 3.39× faster than the NumPy accelerator.**
Wall time includes interpreter startup, the count-table workflow and
full-precision diagnostic output. All nine printed gene summaries were
byte-identical. All permutation p-values and FDR values, including both tails,
matched exactly. Maximum absolute differences were `2.49e-14` for beta,
`3.27e-14` for Wald z-scores and `2.01e-14` for Wald FDR.

The design has one T0 baseline and three T18 samples with one treatment effect.
Settings: updated guide efficiencies, mean-variance modeling using 1,000 genes,
two permutation rounds, seed 42, one worker and one BLAS thread. This cohort
selects labels with exactly four guides before normalization, excluding
incomplete labels and larger control bins. It is not the unfiltered table.

The [machine-readable record](hap1-native-four-guide.json) includes dataset
provenance, hashes, all commands and run manifests, versions, BLAS configuration
and numerical differences. The [successful CI run](https://github.com/dnncha/dotmatch/actions/runs/36724219194)
retains all nine complete result sets. Artifact ID: `11102737825`; archive
SHA-256: `cf30d6e98a92efb8e2252e66592e9c59fd6afcc3bf8e4a2c8645e9c954b097ae`.

To reproduce with a compiled install:

```bash
python packages/crisprworks-fit/benchmarks/public_hap1.py \
  --out-dir hap1-native --cohort four-guide --repeats 3 \
  --kernel native --compare-numpy
```

The native implementation uses a checked SciPy special-function C API,
partial-pivot LU without extra regularization, and compiler settings that
exclude fast-math and floating-point contraction. It omits exact zero design
entries for finite inputs while retaining a dense path for nonfinite inputs.
Tests cover multiple designs, efficiency modes, low counts, matrix adapters,
nonfinite initialization, singular solves, invalid buffers, engine fallback and
state restoration. Linux, macOS, installed wheels and the container passed CI.
Local fitting tests also passed AddressSanitizer and undefined-behavior checks.

The [unfiltered-table parity record](hap1-native-full-parity.json) separately
covers all 71,090 guides and 18,056 labels. It compares a completed CI reference
with local native outputs using identical input hashes, settings and seed.
Printed gene summaries matched byte for byte, all permutation statistics
matched exactly, and maximum absolute differences were below `5.53e-14`.
It is a cross-environment numerical comparison; no timing ratio is claimed.
Upstream's default skip and permutation fallback for the 96-guide label remain.

These results establish a measured performance improvement with numerical
compatibility on this dataset. They do not establish universal performance,
biological superiority or overall state of the art. Chronos, JACKS, legacy
MAGeCK releases and independent multi-screen hit-quality benchmarks remain
outside this comparison. Earlier NumPy timing records use different runners
and must not be combined with these times to calculate a speedup.
