# Reproducible 10x performance control

The `bench-10x-control` target measures the reusable indexed Hamming API on a
deterministic synthetic workload. It is a local hot-path control, not an
end-to-end scientific workflow claim and not a replacement for the public
Sanson/Brunello comparison.

```bash
make bench-10x-control
```

The default control is 1,000,000 reads, 4,096 targets, 20-base windows,
Hamming `k=1`, and five timed repeats. Override the workload when diagnosing a
different lane:

```bash
DOTMATCH_10X_READS=1000000 \
DOTMATCH_10X_TARGETS=4096 \
DOTMATCH_10X_LENGTH=20 \
DOTMATCH_10X_K=1 \
DOTMATCH_10X_REPEATS=5 \
make bench-10x-control
```

Before timing, the harness compares a 256-read prefix against the exhaustive
scan API. Each repeat reports throughput, candidate rates, and a checksum;
different checksums or a failed preflight invalidate the run. The target
library is built once outside the timed region.

For an auditable comparison, record the exact commit, compiler and flags,
machine model, operating-system version, workload arguments, repeat rows,
peak-memory method, and whether other CPU-heavy work was running. A noisy or
partially captured run is a control failure, not a speed result.

The project speed gate remains the repeated full Sanson/Brunello Hamming
`k=1` lane in `docs/10x-goal-baseline.md`. A local synthetic control can
identify a hot path and reject regressions; it cannot establish a 10x
end-to-end improvement by itself.
