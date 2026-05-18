# DotMatch Barcode Panel Design

DotMatch does not merely design barcodes. It designs barcode panels that come with proof of how they will behave during assignment.

`dotmatch panel` is for barcode design for known-target assignment. It is not general genome alignment, not UMI entropy generation, not basecalling, and not a replacement for downstream biological interpretation.

## Commands

```bash
dotmatch panel design --n 96 --length 16 --preset illumina-inline-strict --seed 42 --out-dir panel/
dotmatch panel check panel/barcodes.tsv --k 1 --metric hamming --out-dir panel_check/
dotmatch panel optimize vendor_barcodes.tsv --n 24 --out-dir optimized/
dotmatch panel simulate panel/barcodes.tsv --reads 1000000 --out-dir simulation/
dotmatch panel layout panel/barcodes.tsv --plate 96 --out plate_layout.tsv
dotmatch panel export panel/barcodes.tsv --format illumina-samplesheet --out-dir sample_sheet_templates/
```

Every designed panel includes a machine-checkable safety certificate:

- `panel_summary.json`
- `target_safety.tsv`
- `collision_pairs.tsv`
- `ambiguous_error_spheres.tsv`
- `reverse_complement_warnings.tsv`
- `prefix_collisions.tsv`
- `suffix_collisions.tsv`
- `cycle_balance.tsv`
- `context_risk.tsv`
- `flanked_sequences.tsv`
- `panel_report.html`

## Safety Semantics

Panel safety is checked under DotMatch assignment outcomes: `unique`, `ambiguous`, `none`, and `invalid`.

DotMatch never marks a panel safe for a configured correction radius when any query inside the checked error sphere maps ambiguously. It also reports silent assignment risk, where a possible error from one barcode would be assigned uniquely to another barcode.

The v1 certificate exhaustively enumerates configured error spheres up to `k=2`. Requests above `k=2` are refused rather than downgraded to a partial proof.

For Hamming distance, a minimum pairwise distance `d` gives substitution correction up to `floor((d - 1) / 2)` and detection up to `d - 1`. The certificate still matters because DotMatch checks the exact assignment rules that will later be used by demux/counting.

## Outputs From `panel design`

`panel design` writes:

- `barcodes.tsv`
- `design_report.json`
- `design_trace.tsv`
- `panel_check/`
- `assignment_safety.tsv`
- `collision_pairs.tsv`
- `ambiguous_error_spheres.tsv`
- `flanked_sequences.tsv`
- `plate_layout.tsv`
- `plate_layout.svg`
- `neighbor_distance.tsv`
- `lab_picklist.csv`
- `sample_sheet_templates/`
- `report.html`
- `README_FOR_LAB.md`

The barcode table is intentionally auditable rather than a two-column list. It includes GC, homopolymer length, nearest-neighbor distances, reverse-complement collision flags, self-complement score, status, warnings, and the certified DotMatch command.

## Presets

| Preset | Intended use | Default length | Core guarantee |
| --- | --- | ---: | --- |
| `strict-24x12` | small panels | 12 | high spacing |
| `strict-96x16` | general inline barcodes | 16 | safe one-edit correction target |
| `strict-384x20` | large panels | 20 | high capacity |
| `illumina-inline-96` | fixed-window inline demux | 10 | Hamming-first |
| `illumina-dual-384` | paired i7/i5 indexing | 10+10 | pair-aware safety |
| `nanopore-indel-robust-24` | noisy long-read contexts | 18 | Levenshtein/seqlev spacing |
| `ont-rna004-signal-12` | experimental raw-signal candidate narrowing | 16 | symbolic safety first |

Signal-aware output is experimental and cannot override symbolic safety.
