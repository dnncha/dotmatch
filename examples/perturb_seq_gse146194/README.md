# Replogle Direct-Guide-Capture Perturb-seq Case Study

This is DotMatch's reproducible multi-guide public-data case study. It uses the
32-guide UPR GBC experiment from GSE146194 and a bounded prefix of the dedicated
guide-enrichment run SRR11214031.

The frozen protocol is in `protocol.json`. Machine-readable results and
provenance are in `results.json`, `expected-results.json`, and
`provenance.json`. The human report is
`../../docs/benchmarks/perturb_seq_gse146194/README.md`.

## Reproduce the public run

The workflow downloads the 1.7 MiB primary-paper supplement and streams only
the first 50,000 complete FASTQ records from the 423 MiB ENA archive. It keeps
the 2,000 discovery records out of the 48,000-read evaluation. Allow about
25 MiB of temporary working space; the checked run used 21 MiB.

```bash
make dotmatch
python3 scripts/run_perturb_seq_gse146194.py public --dotmatch ./dotmatch
python3 scripts/check_perturb_seq_gse146194.py --require-public --require-work
```

The full archive MD5 is recorded from ENA but is not locally reverified in
bounded mode. Derived target and read inputs are verified with SHA-256. Raw
reads and the publisher workbook remain in the ignored `work/` directory and
are not redistributed by this repository.

## Run the no-network fixture

```bash
make perturb-seq-case-study-fixture-gate
```

The synthetic fixture checks exact, corrected, ambiguous, unmatched, and
invalid calls. It validates the harness and CLI contract only; it is not
biological evidence.

## Evidence boundary

The public result supports per-read fixed-window assignment against the checked
32-guide barcode list, with matched exact and exhaustive Hamming oracles. It
does not support cell-barcode correction, UMI deduplication, guide-per-cell
calls, expression quantification, perturbation effects, Cell Ranger parity,
full-run prevalence, biological validity, or speed superiority.

The source records are publicly accessible without credentials. This workflow
does not relicense them. [NCBI's molecular-data
policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) places no NCBI
restrictions on use or distribution while noting that NCBI cannot transfer any
rights asserted by submitters. The workflow retrieves the publisher supplement
and archive mirror at run time, cites the paper and accessions, commits only
aggregate results and hashes, and does not redistribute the raw reads or
workbook. Source terms apply to the supplement; this workflow asserts no
redistribution license.
