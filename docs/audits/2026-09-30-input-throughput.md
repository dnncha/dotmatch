# Python input throughput audit

Baseline: `89cc76efecfb803a42e3ce034e855e4b85019d41` (DotMatch 0.6.4).

Two Python input paths did unnecessary work for every read or count cell:

- FASTQ sequence and quality validation ran a Python generator and `ord()` for
  each character. A compiled regex now scans for anything outside ASCII 33–126.
  The parser retains literal symbol handling, sequence uppercasing, structural
  checks, error messages, and hashing before newline/case normalization.
- Ordinary unsigned integer counts went through regex and Decimal parsing.
  ASCII integers of at most 64 characters now use exact `int()` conversion.
  Longer inputs, signed values, decimals, and scientific notation retain the
  existing Decimal path and limits. Unicode digits remain invalid; unselected
  samples are still validated. Counts above 2**53 remain exact.

## Measured evidence

Linux x86_64, CPython 3.12.14; five alternating baseline/candidate runs after
untimed warm-ups. Input generation and imports are outside the timings; parsing
and output hashing are inside. Both versions run in the same process.

| Synthetic workload | Baseline median | Candidate median | Ratio |
| --- | ---: | ---: | ---: |
| 100,000 FASTQ records, 150 bases, lowercase sequence and CRLF | 2.058 s | 0.437 s | 4.71x |
| 10,000 guides x 32 count columns, including integers above 2**53 | 0.468 s | 0.176 s | 2.66x |

Every timed run verifies that the complete parsed output hash matches the
baseline. Raw timings, source hashes, input hashes, and output hashes are in
[2026-09-30-input-throughput.json](2026-09-30-input-throughput.json).

These are synthetic Python parsing measurements. They do not measure whole
assay throughput, native CLI performance, gzip decompression, disk latency,
or peak memory. Deployment workloads may see smaller overall gains.

Reproduce from a source checkout:

```sh
make shared
python3 scripts/bench_python_input.py \
  --baseline-ref 89cc76efecfb803a42e3ce034e855e4b85019d41
```

The benchmark uses the standard library and the source package. Existing
optional package dependencies are imported outside the measured interval.

## Validation

- `make all shared test`: native build and both C test executables pass.
- `make cli-test`: all CLI fixtures and native/Python FASTQ validation parity pass.
- `make python-test`: **1,538 passed**, including optional pandas, Polars,
  AnnData and MultiQC integrations in the cloud environment.
- All **537 new boundary tests** also pass when run against the original
  baseline functions. They cover all byte-range characters except CR/LF
  (existing tests cover line endings), representative Unicode, literal symbols,
  long/large/invalid counts, conversion limits, and unselected columns.
- `git diff --check`: clean.
- `CC=/usr/bin/cc make python-package-test`: reproducible source distribution,
  wheel, clean installs, and installed CLI/evaluation lifecycle checks pass.
  The explicit compiler selects the available GCC rather than the cloud Python
  runtime's absent default `clang` executable.

Validation uses Linux and Python 3.12. The repository's complete GitHub Actions
matrix, other Python versions, macOS, containers and web checks remain required
before merge. The change does not modify native matching or scientific policies.
