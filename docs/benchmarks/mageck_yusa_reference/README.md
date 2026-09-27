# MAGeCK/Yusa Public Reference Spot Check

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

## Claim boundary

This result establishes exact agreement only for six externally displayed
sample counts across three guides. The tutorial does not publish its complete
count table, so this cannot establish whole-table agreement. It also does not
establish which implementation is biologically correct, downstream result
equivalence, performance, release adoption, or independent use of DotMatch.
