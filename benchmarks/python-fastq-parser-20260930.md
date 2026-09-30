# Python FASTQ parser component benchmark — 30 September 2026

The shared Python parser checks sequence and Phred+33 quality strings with a
compiled invalid-character expression instead of Python character loops.
Both fields still accept exactly ASCII 33–126. Record structure, identifiers,
normalization, decompressed input hashing, and error messages are unchanged.
This does not change or measure the native counting engine.

Baseline: `89cc76efecfb803a42e3ce034e855e4b85019d41`.
Raw evidence: `raw/python-fastq-parser-20260930.csv`.

```sh
python scripts/bench_fastq_parser.py \
  --baseline-ref 89cc76efecfb803a42e3ce034e855e4b85019d41 \
  --output benchmarks/raw/python-fastq-parser-20260930.csv
```

Five separate-process repetitions per variant and length, alternating order,
20,000 synthetic records each. Input has lowercase literal DNA, CRLF newlines,
repeated identifiers, comments, and variable printable quality characters.
Every parsed output is checked against the expected record. Python 3.12.14,
Linux x86_64; source/input hashes and process peak RSS are in the CSV. No
resource limits were set. RSS includes interpreter, fixture, and StringIO setup.

| Read length | Baseline median seconds | Candidate median seconds | Ratio |
|---|---:|---:|---:|
| 50 | 0.167502 | 0.062167 | 2.69× |
| 150 | 0.510131 | 0.093594 | 5.45× |
| 300 | 0.889991 | 0.100162 | 8.89× |

These are in-memory Python parsing measurements. They do not establish FASTQ
disk/gzip throughput, native assignment speed, or production pipeline gains.
Longer reads benefit more because the removed interpreter work scales with
sequence and quality length.

The regression suite exercises every byte value and selected Unicode values
in both fields, including embedded newline/control characters and surrogates.
Existing structure, hashing, gzip, and workflow regressions also pass.
