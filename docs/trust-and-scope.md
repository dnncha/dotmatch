# Trust, Scope, and Evidence

DotMatch is built around a small promise: for a configured short read window and
a known target table, it reports how each read relates to the targets under
explicit edit-distance rules. The software is useful when that promise matches
the assay. It should not be stretched into claims about tasks it does not
perform.

## What DotMatch Does

DotMatch supports deterministic known-target assignment for workflows such as:

- CRISPR guide counting before downstream screen statistics;
- inline barcode demultiplexing from FASTQ or FASTQ.gz;
- feature-barcode and guide-capture fixed-window assignment;
- primer-start, amplicon-panel, adapter-prefix, and whitelist-style assays;
- barcode panel design, checking, simulation, and lab export;
- target-library audits before one-edit or multi-edit correction.

Each read is assigned one of four states:

- `unique`: exactly one target is compatible and can be counted or written to
  the matching output.
- `ambiguous`: more than one target is compatible, so the read should not be
  forced into a single target count without an explicit compatibility policy.
- `none`: no target is close enough.
- `invalid`: the requested read window could not be extracted.

## What DotMatch Does Not Do

DotMatch is not a genome aligner, basecaller, adapter trimmer, cell/UMI
quantifier, variant caller, consensus builder, or screen-level statistical
analysis tool. It does not produce SAM/BAM, CIGAR strings, expression matrices,
variant calls, guide-level effect sizes, or hit-calling statistics.

DotMatch also does not currently make posterior-probability or calibrated
likelihood calls. `--max-correction-qual` is a deterministic Phred-quality
filter for selected one-edit correction paths; it can reject high-quality
observed mismatches or insertions, but it is not a probabilistic confidence
score.

For CRISPR screens, DotMatch stops at guide assignment, count matrices, and QC.
Use tools such as MAGeCK, BAGEL, drugZ, CERES, or other appropriate downstream
methods for biological effect analysis.

## Why Ambiguity Is Not Hidden

Many short-target workflows fail quietly when a read is compatible with more
than one target. DotMatch treats that as a first-class output state. Under the
default `radius` policy, a read is counted only when exactly one target lies
inside the configured radius.

This makes correction auditable. If one-mismatch rescue would cause collisions,
the target audit and run outputs show that risk rather than converting it into
apparently confident counts.

## Evidence Boundaries

The repository keeps benchmark reports, public-data examples, comparator notes,
and gate scripts alongside the code. Public claims should stay within that
checked evidence.

Important evidence pages:

- [Evidence Notes](scientific-claims.md): current supported, experimental, and
  intentionally unsupported claims.
- [Benchmark Overview](benchmarks/README.md): reports with data sources,
  commands, comparator settings, and checked outputs.
- [Barcode Validation Notes](barcode-science-readiness.md): public fixed-window
  barcode examples and failure-mode checks.
- [Native Comparator Scope](native-comparator-scope.md): exactly which native
  comparator evidence is currently available.
- [Methods and Citation](methods-and-citation.md): methods text that matches
  the current assignment semantics.

## Practical Review Checklist

Before trusting a production run, check:

1. The read window start and length match the assay design.
2. The target table has no duplicate or unsafe near-neighbor records for the
   selected correction radius.
3. The edit metric matches the biology: Hamming for fixed-length substitution
   correction, Levenshtein when short indels are intended.
4. The ambiguity policy is recorded and appropriate.
5. `sample_qc.tsv` does not show unexpected ambiguous, unmatched, invalid, or
   low-assignment fractions.
6. Any public or comparative statement is backed by a report and gate listed in
   the evidence pages.
