# GSE146194 Direct-Guide-Capture Perturb-seq Case Study

This bounded public-data case study tests DotMatch at its evidence boundary: per-read fixed-window assignment to a known multi-guide barcode library. It uses the UPR GBC sample from Replogle et al. and holds window-discovery reads out of the reported evaluation.

## Result

The publisher supplement yielded `32` guide barcodes. A frozen discovery rule selected `forward` orientation at zero-based start `41` for an `18`-base window, then `48,000` held-out reads were evaluated.

| Rule | Unique | Exact | Corrected | Ambiguous | Unmatched | Invalid | Distinct guides | Oracle mismatches | Agreement |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Hamming k=0 radius | 46189 | 46189 | 0 | 0 | 1811 | 0 | 32 | 0 | 100.00% |
| Hamming k=1 radius | 46910 | 46189 | 721 | 0 | 1090 | 0 | 32 | 0 | 100.00% |

The unmatched and ambiguous columns are part of the result, not discarded failures. The small deterministic fixture separately requires exact, corrected, ambiguous, unmatched, and invalid outcomes. If this public slice contains zero ambiguous reads, that is reported as zero rather than manufactured.

## Scientific question and criterion

**Question.** On an accessioned multi-guide direct-capture Perturb-seq guide-enrichment run, can DotMatch reproduce independent fixed-window per-read guide-barcode assignments while keeping ambiguous and unmatched reads explicit?

**Ground truth.** For k=0, a transparent exact-slice hash classifies target multiplicity. For k=1, an independent exhaustive Hamming-radius implementation checks every guide. The frozen pass criterion is zero held-out per-read differences in status, target identifier, and distance.

The authors' published guide caller is not used as a per-read comparator: it combines guide-aligned reads with Cell Ranger-corrected cell barcodes and UMIs, then makes cell-level threshold calls. Those units and rules are not matched to this bounded per-read assignment question.

## Dataset and provenance

- GEO series/sample: [GSE146194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146194) / [GSM4367979](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4367979)
- SRA study/experiment/run: [SRP251252](https://www.ncbi.nlm.nih.gov/sra/?term=SRP251252) / [SRX7826824](https://www.ncbi.nlm.nih.gov/sra/SRX7826824) / [SRR11214031](https://www.ncbi.nlm.nih.gov/sra/SRR11214031)
- BioProject/BioSample: [PRJNA609688](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA609688) / [SAMN14258014](https://www.ncbi.nlm.nih.gov/biosample/SAMN14258014)
- Primary paper: Replogle JM et al., *Nature Biotechnology* 38, 954-961 (2020), [doi:10.1038/s41587-020-0470-y](https://doi.org/10.1038/s41587-020-0470-y), [PMID 32231336](https://pubmed.ncbi.nlm.nih.gov/32231336/)
- Access and reuse: [NCBI's molecular-data policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) places no NCBI restrictions on use or distribution but does not transfer any rights asserted by submitters. The [publisher supplement](https://support.springernature.com/en/support/solutions/articles/6000210902-supplementary-information) is retrieved from its article link; this workflow asserts no redistribution license and does not commit the workbook.

Guide workbook SHA-256:

```text
3dd53733987890bb7577c3dac77215c9874267d4d69f732ddf7525f356914085
```

Target table SHA-256:

```text
cc17aef8352611de2fa3d51007792b2721d159205ab393b1268db2486d499629
```

Evaluation FASTQ SHA-256:

```text
93ed9633f44e2ba0c60b1312492e90b101aa6249ea45766ce753c2908b410fe5
```

Selected-prefix uncompressed SHA-256:

```text
8346a01785a1cf65e38789b3ae703353734e24df9bf14775ea31e7825dab7b01
```

The raw read archive is not committed. The workflow records ENA's full-file MD5 and byte count, streams only the frozen prefix, and verifies the derived prefix and evaluation files with SHA-256. The full archive MD5 is not claimed as locally reverified in bounded mode.

## Reproduce

```bash
make dotmatch
python3 scripts/run_perturb_seq_gse146194.py public --dotmatch ./dotmatch
python3 scripts/check_perturb_seq_gse146194.py --require-public
```

For the no-network deterministic harness:

```bash
make perturb-seq-case-study-fixture-gate
```

Expected public outputs are `examples/perturb_seq_gse146194/results.json`, `examples/perturb_seq_gse146194/provenance.json`, `benchmarks/raw/perturb_seq_gse146194.csv`, and this report. Large reads and per-read work products stay under the ignored `examples/perturb_seq_gse146194/work/` directory.

## Methods

The workflow extracts `32` 18-base GBCs from Supplementary Table 2, scans every valid start and both target orientations on the first `2000` reads, and applies deterministic score and tie-break rules from `protocol.json`. Those discovery reads are excluded. DotMatch uses Hamming distance with conservative radius ambiguity at k=0 and k=1. Independent matched-rule oracles produce per-read status, target, and distance for the held-out prefix. Commands, versions, hashes, and resources are recorded in machine-readable artifacts.

## QC interpretation

- k=0: assignment `96.23%`, ambiguous `0.00%`, unmatched `3.77%`, invalid `0.00%`; `32` distinct guides received unique reads.
- k=1: assignment `97.73%`, ambiguous `0.00%`, unmatched `2.27%`, invalid `0.00%`; `32` distinct guides received unique reads.
- Guide-library audit: minimum pairwise Hamming distance `8`; `0` pairs within distance 1 and `0` pairs within distance 2.
- Resource record for this bounded run: `12.708` wall seconds, `133.7` MiB peak RSS, and `1094540` compressed bytes read from the FASTQ stream.

## What this proves

This proves reproducible, matched-rule per-read assignment and explicit QC behavior for the checked held-out SRR11214031 prefix and the published 32-guide GBC list. It closes the earlier single-guide extraction gap with a real multi-guide direct-capture dataset.

## What this does not prove

It does not prove guide-per-cell accuracy, UMI deduplication, Cell Ranger parity, expression quantification, perturbation effects, biological validity, full-run prevalence, or a speed advantage. The bounded prefix and fixed-position assumption can fail to represent other runs or protocols.

## Next step

A core facility or Perturb-seq workflow maintainer can reuse the manifest and held-out protocol on a complete guide-enrichment run, then add matched cell-barcode/UMI aggregation and compare guide-per-cell calls with the authors' published `cell_identities.csv` under explicitly matched thresholds.
