# Amplicon/Panel Assignment Evidence

This report covers panel-style target assignment evidence for DotMatch's known-target counting layer.

The synthetic lane checks exact, ambiguous, and unmatched fixed-window target assignment. The public lane uses nf-core viralrecon Illumina ARTIC V3 amplicon sample R1 and validates DotMatch k=0 against a transparent exact-prefix hash baseline over the selected full-primer length group.

Current status: public primer-start assignment evidence only. This is not amplicon consensus, variant calling, primer trimming, or clinical validation evidence.

## Synthetic Command

```bash
dotmatch count --targets benchmarks/work/amplicon_panel/panel_targets.tsv --reads benchmarks/work/amplicon_panel/panel_reads.fastq --sample-label amplicon_panel_fixture --target-start 0 --target-length 12 --k 1 --metric hamming --format dotmatch --out benchmarks/work/amplicon_panel/panel_counts.tsv --summary benchmarks/work/amplicon_panel/panel_summary.json --assignments benchmarks/work/amplicon_panel/panel_assignments.tsv --ambiguous report --sample-qc benchmarks/work/amplicon_panel/panel_sample_qc.tsv
```

## Raw Rows

| tool | workflow | status | targets | reads | start | length | k | metric | assigned | exact | corrected | ambiguous | unmatched | validation mismatches |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| dotmatch_count | synthetic_amplicon_panel_fixture | smoke | 4 | 6 | 0 | 12 | 1 | hamming | 4 | 4 | 0 | 1 | 1 | 0 |
| dotmatch_count | public_nfcore_artic_v3_amplicon_primer | supported | 80 | 20000 | 0 | 22 | 0 | hamming | 6454 | 6454 | 0 | 0 | 13546 | 0 |
| dotmatch_count | public_nfcore_artic_v3_amplicon_primer | supported | 80 | 20000 | 0 | 22 | 1 | hamming | 6545 | 6454 | 91 | 0 | 13455 | 0 |
| exact_prefix_hash | public_nfcore_artic_v3_amplicon_primer | supported | 80 | 20000 | 0 | 22 | 0 | exact | 6454 | 6454 | 0 | 0 | 13546 | 0 |

## Public Amplicon/Panel Lane

- Dataset: nf-core/test-datasets viralrecon Illumina amplicon `sample1_R1.fastq.gz`.
- Primer panel: ARTIC V3 SARS-CoV-2 primer FASTA from the same nf-core/test-datasets branch.
- Fixture semantics: the fetcher selects the full-primer length group with the most exact R1 prefix assignments, then counts fixed-position primer starts.
- Comparator semantics: the exact-prefix baseline counts reads whose R1 prefix exactly matches the selected primer sequence. It validates per-read known-primer assignment semantics, not consensus generation, primer trimming, variant calling, or clinical interpretation.

## Public Commands

```bash
dotmatch count --targets examples/amplicon_panel/data/artic_v3_primers_len22.tsv --reads examples/amplicon_panel/data/nfcore_viralrecon_sample1_R1.subsample20000.fastq.gz --sample-label nfcore_viralrecon_artic_v3_sample1_R1 --target-start 0 --target-length 22 --k 0 --metric hamming --format dotmatch --out benchmarks/work/amplicon_panel/public_nfcore_amplicon_k0_counts.tsv --summary benchmarks/work/amplicon_panel/public_nfcore_amplicon_k0_summary.json --assignments benchmarks/work/amplicon_panel/public_nfcore_amplicon_k0_assignments.tsv --ambiguous report --sample-qc benchmarks/work/amplicon_panel/public_nfcore_amplicon_k0_sample_qc.tsv
```

```bash
dotmatch count --targets examples/amplicon_panel/data/artic_v3_primers_len22.tsv --reads examples/amplicon_panel/data/nfcore_viralrecon_sample1_R1.subsample20000.fastq.gz --sample-label nfcore_viralrecon_artic_v3_sample1_R1 --target-start 0 --target-length 22 --k 1 --metric hamming --format dotmatch --out benchmarks/work/amplicon_panel/public_nfcore_amplicon_k1_counts.tsv --summary benchmarks/work/amplicon_panel/public_nfcore_amplicon_k1_summary.json --assignments benchmarks/work/amplicon_panel/public_nfcore_amplicon_k1_assignments.tsv --ambiguous report --sample-qc benchmarks/work/amplicon_panel/public_nfcore_amplicon_k1_sample_qc.tsv
```

```bash
python3 scripts/bench_amplicon_panel.py --include-public --metadata examples/amplicon_panel/data/metadata.json
```


## Evidence Boundary

Use these lanes to verify fixed-window known-target panel assignment plumbing, explicit ambiguity handling, and narrow public ARTIC primer-start per-read assignment. Broader amplicon/panel benchmark wording requires public full-assay comparator semantics, consensus or variant-call validation where relevant, exact commands, validation artifacts, and a passing gate.
