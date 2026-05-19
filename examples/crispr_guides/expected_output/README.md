# Expected Output

The small fixture is deterministic. Running `../run.sh` should produce files that match the tracked TSV snapshots here:

- `counts.tsv`
- `counts.mageck.tsv`
- `assignments.tsv`
- `ambiguous.tsv`
- `unmatched.tsv`
- `mageck_skipped.txt`

`summary.json` includes timing fields, so `summary.stable.json` records only the stable assignment fields that should match across machines.

The exact full-data counts are not checked into the repo because the public FASTQ files are large.
