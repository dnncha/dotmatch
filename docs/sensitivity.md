# Assignment-policy sensitivity

The `dotmatch sensitivity` command, introduced in 0.5.0, compares exact, Hamming radius-one,
and best-distance Hamming assignments on the same fixed windows of one FASTQ.
Install the matching release with `python3 -m pip install dotmatch==0.5.0`.

```bash
dotmatch sensitivity \
  --targets guides.tsv --reads sample.fastq.gz \
  --target-start 23 --target-length 20 \
  --sample-label sample_1 --out-dir sensitivity/
```

Open `sensitivity/report.html`. Start with read outcomes and per-guide deltas,
not the assigned fraction alone. No policy is automatically selected.

## Semantics

All three policies use the same library, orientation and fixed extraction window.
The comparison uses the native Hamming index, with one query per distinct valid
window in each batch. The exact and radius results are derived from the full
radius-one result, not from a heuristic candidate sample. Batches bound read
memory; target/index memory scales with the supplied library.

`exact` requires one target at distance zero. `radius_k1` requires exactly one
target anywhere within distance one. `best_k1` requires one nearest target within
distance one. Ties and duplicate sequences remain ambiguous; an exact hit with a
nearby alternative can be unique under best-distance but ambiguous under radius.

Only unique calls add counts. Short windows are invalid, not silently dropped.
N and IUPAC symbols are literal symbols, not wildcards. Quality characters are
validated but do not weight assignments. There is no indel model, reverse-
complement search, offset inference, cell/UMI processing or downstream statistics.

## Files

| File | Contents |
| --- | --- |
| `report.html` | Self-contained review, no external scripts or tracking |
| `{exact,radius_k1,best_k1}.counts.tsv` | Raw guide counts with `sgRNA`, `Gene`, and sample columns |
| `guide_deltas.tsv` | Every supplied target, including unchanged and zero-count targets |
| `transitions.tsv` | Complete four-state transitions for each policy pair |
| `sample_qc.tsv` | All four outcomes under each policy |
| `summary.json` | Versioned completion manifest, settings, counts, hashes and comparison scope |
| `read_changes.tsv` | Optional changed record ordinals/IDs and policy calls; no raw sequence or quality |

Enable the last file with `--write-read-changes`. Read IDs can repeat; the 1-based
record ordinal disambiguates occurrences. Report tables show at most 50 changed
guides; TSVs retain all guides. The source files are never altered.

The summary uses schema `dotmatch.sensitivity.v1`. `changed_reads` counts records
whose status or unique target differs across policies. `counts_identical` compares
all per-guide counts, not just the grand total. Equal count matrices can still
have different read membership; these are separate observations.

The FASTQ digest is calculated over the original file bytes during the input
pass (compressed bytes for gzip). The report records the analysis implementation
and native binary hashes as well as the package version. Those distinguish a
source build from a published package with the same version string. Artifact
hashes exclude the summary itself; do not invent a self-referential checksum.

## Failure handling

The output directory must not already exist. Unsupported configuration, duplicate
or empty target IDs, truncated FASTQs and invalid qualities stop the run. An empty
FASTQ is rejected. Duplicate target *sequences* with distinct IDs are retained.
All target lengths must equal `--target-length`.

Work is staged inside an exclusively created directory. Final files are linked
without replacing existing paths, with the completion manifest published last.
An interrupted process can leave a pending directory; do not treat it as a
completed report without `summary.json` containing `completion: complete`.
Normal input errors clean up the run's temporary files. Default batch size is
4,096, configurable from 1 to 65,536. No data is uploaded.

## Reproduce the example

```bash
python3 -m pip install dotmatch==0.5.0
dotmatch sensitivity \
  --targets examples/assignment_sensitivity/targets.tsv \
  --reads examples/assignment_sensitivity/reads.fastq \
  --target-start 0 --target-length 20 \
  --write-read-changes --out-dir sensitivity-example
```

The nine-read synthetic example has exact outcomes 3 unique, 1 ambiguous,
4 unmatched, 1 invalid; radius-one outcomes 3, 4, 1, 1; and best-distance outcomes
5, 2, 1, 1. Five reads change outcomes across policies. The website's example is
regenerated and checked with `python scripts/generate_assignment_demo.py --check`.
These are algorithmic checks, not biological validation or a competitor benchmark.
