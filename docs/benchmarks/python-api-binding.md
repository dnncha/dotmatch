# Python API Binding Throughput

This is a local developer measurement of the Python binding path: argument
marshalling, native kernel, result-object construction, and (for the stream
case) FASTQ parsing. It is not a public claim, a cross-machine guarantee, or a
replacement for the workflow benchmarks in the evidence notes.

Harness:

```bash
PYTHONPATH=python DOTMATCH_LIB=$PWD/libdotmatch.dylib \
  python scripts/bench_python_api.py --workload hamming --repeats 11 --json
```

Workload: 40,000 reads, 4,096 targets, length 20, 70% one-substitution reads,
median of 31. The stream case feeds the same reads through `stream_assign` from
a temporary FASTQ. Host: Apple M5, macOS, CPython 3.12.12. Raw data with every
sample, min, and max: `benchmarks/raw/python_api_binding_20261027.json`.

| workload | baseline reads/s | current reads/s | change |
|---|---:|---:|---:|
| Hamming, k=1 | 1,273,487 | 1,683,587 | +32.2% |
| Levenshtein, k=1 | 596,400 | 699,997 | +17.4% |
| Exact, k=0 | 1,297,718 | 1,968,851 | +51.7% |
| Stream Hamming | 296,564 | 395,915 | +33.5% |

Changes measured: result records became `NamedTuple`s with a status lookup
table; read marshalling got a fast `str`/`bytes` path and result records use
`tuple.__new__`; `stream_assign` batches fixed-length reads at 256 instead of
4096.
