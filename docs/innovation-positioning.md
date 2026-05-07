# DotMatch Innovation Positioning

DotMatch should not claim to invent edit distance, approximate DNA matching, or CRISPR guide counting.

The precise innovation target is:

> DotMatch is a general indexed exact-assignment engine for known short-DNA target sets. It supports `k=1` Levenshtein assignment, including substitutions and single-base insertions/deletions, while preserving exhaustive-oracle `unique`, `ambiguous`, and `no-match` semantics and emitting workflow-ready FASTQ count outputs.

## What Is Not New

- Edit distance and Myers-style bit-vector alignment are established; Edlib is the exact pairwise oracle.
- Error-tolerant adapter/short-sequence matching is established; Cutadapt is a relevant workflow comparator.
- Trie/Levenshtein search and barcode clustering are established; Starcode-style methods are relevant algorithmic context.
- CRISPR guide counting with one mismatch is established; guide-counter is a serious comparator.

## Strong DotMatch Wedge

DotMatch should compete on the combination of:

- exact `k=1` Levenshtein semantics, not only Hamming/one-mismatch semantics;
- deterministic assignment states: `unique`, `ambiguous`, `no-match`;
- candidate-verification collapse: Edlib scan verifies all targets per read; DotMatch should verify a small candidate set;
- target-set audit for safe one-edit correction;
- native FASTQ/FASTQ.gz streaming and count-table outputs;
- explicit `hamming` and `levenshtein` count modes so guide-counter-style comparisons do not blur semantics;
- optional guide-offset detection for CRISPR FASTQ workflows;
- reusable native CLI, C ABI, and Python API across CRISPR, barcode, amplicon, and whitelist-style workflows.

Current algorithmic implementation:

- `k=0` uses exact encoded lookup for A/C/G/T targets up to 32 bp.
- `k=1` uses exact neighbor-key lookup for substitutions and one-base insertions/deletions, then verifies every candidate with exact Levenshtein distance before assignment.
- Non-ACGT and longer targets fall back to exact verification paths.

This is a strong systems implementation, but it should be described as indexed candidate generation unless later work proves a genuinely novel algorithmic contribution.

## Competitor Framing

| Tool class | Role in paper | Claim boundary |
| --- | --- | --- |
| Edlib scan | exact semantic oracle | pairwise exhaustive baseline, not a workflow tool |
| guide-counter | CRISPR one-mismatch workflow comparator | mismatch-only/no-indel comparator per current docs |
| Cutadapt | adapter/search workflow comparator | useful where adapter semantics match, not assignment oracle |
| Bowtie2 | aligner workflow comparator | over-general mapping workflow, not known-target assignment oracle |
| hash lookup | exact `k=0` baseline | exact-only baseline, not a `k=1` competitor |
| BK-tree / neighbor lookup | approximate lookup baselines | algorithmic baselines for candidate pruning |

## Paper Claim

Preferred claim:

> DotMatch turns exact one-edit known-target assignment into a lookup-like operation while preserving exhaustive Edlib assignment semantics and producing auditable count outputs for real FASTQ workflows.

Avoid:

- "novel edit-distance algorithm" unless a genuinely new algorithm is added;
- "replacement for guide-counter" unless the benchmark directly proves it for guide-counter's supported semantics;
- "state of the art aligner" or broad alignment claims.

## What Would Make The Claim Hard To Dismiss

- Real-data benchmark against native Edlib scan with zero mismatches.
- Workflow comparison against guide-counter on a CRISPR dataset, clearly noting guide-counter's mismatch-only/no-indel semantics.
- Synthetic and real indel-containing tests showing DotMatch handles single insertions/deletions where mismatch-only counters do not.
- Target audit report showing unsafe libraries and ambiguity risk.
- Reproducible figures showing verified candidates/read remains near flat as target count grows.
