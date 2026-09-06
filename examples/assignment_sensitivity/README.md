# Why equal totals can hide different guide counts

This **synthetic software fixture** contains nine reads and five target IDs,
including two IDs for the same sequence. It exercises a close pair, an isolated
target, substitutions, a literal N, an unmatched read and an invalid short window.
There is no experimental biological claim attached to these data.

From a checkout of the v0.5.0 release:

```bash
python3 -m pip install dotmatch==0.5.0
dotmatch sensitivity \
  --targets examples/assignment_sensitivity/targets.tsv \
  --reads examples/assignment_sensitivity/reads.fastq \
  --target-start 0 --target-length 20 \
  --write-read-changes --out-dir sensitivity-example
```

Open `sensitivity-example/report.html`. Use a new output directory
for each run; existing directories are deliberately refused.

| Policy | Unique | Ambiguous | Unmatched | Invalid |
| --- | ---: | ---: | ---: | ---: |
| Exact | 3 | 1 | 4 | 1 |
| Radius k=1 | 3 | 4 | 1 | 1 |
| Best distance, k=1 | 5 | 2 | 1 | 1 |

Exact and radius-one have the same unique total but different per-guide counts.
Five reads change state across the three policies. The higher best-distance count
is not evidence of greater biological accuracy. Inspect `read_changes.tsv` and
`guide_deltas.tsv` to understand each transition.

The landing page consumes `public/assignment-demo.json`, generated from these
inputs by `scripts/generate_assignment_demo.py`. CI checks it against the actual
native matcher and independently enumerates candidate targets. The example is
not an unexecuted mock report.
