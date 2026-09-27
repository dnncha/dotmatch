# MAGeCK/Yusa Public Reference Checks

DotMatch 0.6.1 exactly reproduced the six sample counts displayed in MAGeCK's
public Yusa tutorial: three guides across the `plasmid` (ERR376998) and `ESC1`
(ERR376999) samples. This is a small, independently published-output check, not
a whole-table comparison.

| guide | gene | plasmid | ESC1 | result |
|---|---|---:|---:|---|
| chr19:5884430-5884453 | SLC25A45 | 13 | 32 | match |
| chr11:58831475-58831498 | OLFR312 | 94 | 108 | match |
| chr4:49282352-49282375 | E130309F12RIK | 85 | 128 | match |

## Protocol

- Baseline: the immutable DotMatch 0.6.1 release commit
  `5c82bd405fc49d52e04f429fab429d044ec86997`.
- External reference: the three count rows displayed by the
  [MAGeCK public-data tutorial](https://sourceforge.net/p/mageck/wiki/demo/).
- Inputs: the tutorial's ERR376998/ERR376999 FASTQs, with ENA MD5s verified,
  and MAGeCK's Yusa library.
- Semantics: exact 19-base guide counting at fixed offset 23, with the tutorial
  sample labels `plasmid,ESC1`.
- Failure criteria: any accession/hash mismatch, missing guide, gene mismatch,
  or count mismatch.

The run processed 20,394,663 reads and 87,437 library guides. All three rows
and all six displayed sample counts matched. The complete DotMatch count table
has SHA-256
`88ad3ef533f6c285d36e06fd55a5b562a738a8413512d51e03fc967dc9cbb4c3`;
it is not committed because it is a derived 87,437-row artifact.

The checked details and provenance are recorded in
`benchmarks/raw/mageck_yusa_public_reference_spotcheck.csv` and
`benchmarks/raw/mageck_yusa_public_reference_spotcheck.json`. Recheck any
MAGeCK-format output with:

```bash
python3 scripts/check_mageck_yusa_reference_rows.py counts.tsv
```

## Complete nf-core subset comparison

A separate whole-table check uses nf-core's committed 10,000-read subsets of
the same two accessions and its independently recorded MAGeCK 0.5.9.5 output
MD5. A local MAGeCK 0.5.9.5 rerun reproduced that MD5 exactly. After joining by
`sgRNA` to make row order irrelevant, DotMatch 0.6.1 matched all 87,437 guide
IDs, every gene value, and all 174,874 count cells; there were no missing,
extra, or mismatching values. MAGeCK and DotMatch both assigned 8,499 reads per
input.

The raw files have different compressed hashes but decompress to identical
10,000-record FASTQs. The two columns are therefore duplicate software checks,
not independent biological samples. This result strengthens implementation
agreement on a complete external fixture while leaving the full-run external
table gate open.

Pinned inputs, tool artifacts, commands, hashes, results, and this negative
finding are recorded in
`benchmarks/raw/mageck_yusa_nfcore_subset_full_table.json`. Recheck complete
tables without relying on their row order with:

```bash
python3 scripts/compare_mageck_count_tables.py \
  mageck.count.txt dotmatch.count.txt --samples test,test2
```

## Complete full-run implementation comparison

The full public workflow provides a stronger implementation check. A pinned
MAGeCK 0.5.9.5 package and the immutable DotMatch 0.6.1 release independently
counted the complete verified ERR376998 and ERR376999 FASTQs against the same
Yusa library. The run processed 20,394,663 reads. Both tools assigned
8,615,587 plasmid reads and 8,475,790 ESC1 reads.

After joining by `sgRNA`, all 87,437 guide IDs, every gene value, and all
174,874 integer count cells matched exactly. There were no missing, extra, or
mismatching values. As a separate check, preserving the header and sorting
both tables by guide produced byte-identical canonical files with SHA-256
`dabac8493ef031802522f18194c43617d0b35193748aa1cedd938c1369cfc500`.

The complete commands, ENA MD5s, input and tool-package hashes, output hashes,
mapped totals, runtime, failure criteria, and instrumentation limitation are
recorded in
`benchmarks/raw/mageck_yusa_full_implementation_comparison.json`. The derived
87,437-row tables are not committed. Recheck locally produced full tables with:

```bash
python3 scripts/compare_mageck_count_tables.py \
  escneg.count.txt dotmatch-0.6.1.counts.tsv --samples plasmid,ESC1
```

## Claim boundary

The tutorial result supplies independently published values for six displayed
sample counts across three guides. The nf-core result supplies a committed
external whole-table artifact for a 10,000-read subset. The complete full-run
comparison establishes cross-implementation agreement, but its MAGeCK output
was produced during this validation rather than published independently.
Together these checks do not determine which implementation is biologically
correct, validate downstream result equivalence, provide a controlled
performance benchmark, or demonstrate release adoption or independent use of
DotMatch.
