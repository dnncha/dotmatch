# ML Backend Optimizer Design

## Goal

Add a benchmark-informed optimizer around DotMatch AssaySpec workflows so the system can recommend the fastest safe execution backend without allowing ML to change scientific assignments.

## Scope

The optimizer is advisory for this release. CPU deterministic assignment remains the authority. The optimizer can recommend an experimental Apple Metal GPU lane only when the assay is compute-compatible, public GPU evidence is validated for that assay class, and the plan records that CPU validation remains required.

## Architecture

- `python/dotmatch/assayspec.py` owns the first implementation because AssaySpec already validates backend mode, builds reliability summaries, and records GPU eligibility.
- The optimizer extracts workload features from the spec: mode, assay type, metric, edit radius, target count, target length, A/C/G/T packability, GPU evidence status, and optional read estimates.
- A small local benchmark-informed scorer uses checked evidence from current GPU benchmark rows and deterministic rules. This is intentionally not a cloud model and not a learned assignment classifier.
- `dotmatch assay optimize assay.toml` writes `backend_optimization.json` and prints a concise human-readable plan.
- Reliability summaries include optimizer metadata so assay reports explain why CPU or GPU was selected or deferred.

## Boundaries

- No read assignment is performed by ML.
- No automatic production GPU switch happens in this release.
- GPU recommendations must include an accuracy gate: CPU checksum/count agreement is required before GPU speed is considered useful.
- Non-Hamming, `k != 1`, non-A/C/G/T, variable-length, pair-count, and unsupported assay evidence cases remain CPU-required.

## Success Criteria

- CRISPR Hamming `k=1` fixed-window A/C/G/T count specs report a GPU candidate with an expected speedup band.
- Inline barcode demux specs that are compute-compatible but lack public GPU gates report GPU as gated, not selected.
- Levenshtein and pair-count specs report CPU-required with clear reason codes.
- `dotmatch assay optimize` writes JSON and human text.
- Existing assay check/run behavior remains deterministic and CPU-authoritative.
