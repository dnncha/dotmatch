# Closed-Loop Assay Reliability Design

Date: 2026-05-21

## Purpose

DotMatch should deliver the most value by becoming the evidence-bounded control
plane for known-target sequencing assays. The product goal is to prevent failed
runs, run assignments quickly, explain suspicious outputs, and hand scientists a
report they can trust.

The first shippable slice is closed-loop assay reliability for fixed-window
known-target assays. It extends the existing AssaySpec, panel certificate,
autopsy, public evidence, and GPU benchmark systems instead of creating a
parallel workflow.

## Product Promise

For a CRISPR, feature-barcode, inline-barcode, guide-capture, amplicon-panel, or
adapter-prefix run, a user should be able to:

1. describe or infer the assay;
2. certify whether the target/barcode set is safe for the intended correction
   radius;
3. run assignment with deterministic semantics;
4. detect likely configuration or assay failures early;
5. receive an evidence-bounded report with supported claims and explicit
   boundaries.

This is not a promise to replace MAGeCK, Cell Ranger, bcl-convert, a genome
aligner, a basecaller, or downstream biological interpretation tools.

## Approaches Considered

### Recommended: Reliability Layer On Existing AssaySpec

AssaySpec remains the user-facing contract. Reliability analysis is added as a
checked preflight, runtime QC, autopsy, and report layer around existing native
commands.

This is the best path because it reuses the current tested surfaces, keeps CLI
and Workbench aligned, and makes the new value available to workflow systems.

### Alternative: GPU-First Product Track

A GPU-first track would prioritize Metal/CUDA kernels and streaming throughput
before improving assay correctness or reporting.

This is useful later, but it is too narrow as the main product bet. Faster wrong
assumptions do not save failed experiments. GPU acceleration should be one
execution backend inside reliability-gated workflows.

### Alternative: ML-First Assay Intelligence

An ML-first track would infer assay type, offsets, failure modes, and biological
signals using learned models.

This can add value around the deterministic core, but it should not replace
assignment semantics. The first ML uses should be anomaly detection, offset
recommendations, and report explanation where model output is advisory and
auditable.

## First Shippable Slice

The first slice adds an assay reliability report and machine-readable reliability
summary to `dotmatch assay check`, `dotmatch assay plan`, and `dotmatch assay
run`.

New artifacts:

- `reliability_summary.json`
- `reliability_findings.tsv`
- `reliability_report.html`
- `reliability_manifest.summary.tsv`

The report aggregates existing evidence where possible:

- AssaySpec validation status;
- target audit and panel safety certificate status;
- correction-radius ambiguity risk;
- inferred or configured fixed-window confidence;
- sample-level assignment, ambiguous, unmatched, invalid, and correction rates;
- automatic autopsy findings;
- public-evidence claim boundary for the assay type;
- backend selection summary, including whether GPU was eligible, used, skipped,
  or unavailable.

## Architecture

### AssaySpec Reliability Contract

AssaySpec gains an optional `[reliability]` section. Defaults are conservative
and match current automatic autopsy behavior.

Example:

```toml
[reliability]
profile = "production"
fail_on_unsafe_targets = true
fail_on_draft_inference = true
min_assignment_rate = 0.80
max_ambiguous_rate = 0.05
max_unmatched_rate = 0.15
max_invalid_rate = 0.02
require_public_evidence_boundary = true

[backend]
mode = "auto"
allow_gpu = true
```

`profile = "exploratory"` records findings without failing the run.
`profile = "production"` fails on unsafe preflight conditions and marks runtime
threshold violations as failed reliability checks after outputs are written.

### Reliability Engine

Add a small reliability engine in the Python workflow layer. It should consume
existing artifacts instead of reimplementing assignment:

- parsed AssaySpec;
- target audit outputs;
- panel certificate outputs when present;
- `sample_qc.tsv`;
- `summary.json`;
- CRISPR QC artifacts where applicable;
- autopsy outputs;
- GPU evidence and backend-selection metadata when present;
- `docs/assay-evidence.json` claim boundaries.

The engine emits normalized findings:

```text
finding_id
severity: info | warning | error | blocked
stage: preflight | run | postrun | evidence
sample_id
metric
observed
threshold
message
recommended_action
source_artifact
```

### Backend Selector

Backend selection should be explicit and recorded. Initial modes:

- `cpu`: current production path;
- `gpu-metal-experimental`: eligible only for fixed-length A/C/G/T Hamming
  `k=1` lanes with a validated CPU checksum gate;
- `auto`: starts with CPU production path and records GPU eligibility without
  using GPU unless the workflow has an enabled and passing real-workload gate.

The first production slice should not silently switch assignment to GPU. It
should report GPU eligibility and keep CPU as the authority until the GPU path
has fallback handling, reproducibility checks, and real workload gates for the
specific assay type.

### Workbench Integration

Workbench should read the same reliability artifacts. It does not need a new
analysis engine.

Expected views:

- preflight readiness;
- run status and failed thresholds;
- sample-level findings;
- target or barcode safety;
- backend eligibility;
- evidence boundary and supported claims.

Workbench remains local-only and continues to run explicit DotMatch commands.

## Data Flow

1. `dotmatch assay check assay.toml` validates AssaySpec and runs reliability
   preflight checks that do not require reads.
2. `dotmatch assay plan assay.toml` prints native commands and includes a
   reliability plan section.
3. `dotmatch assay run assay.toml` runs target audit, assignment, QC, optional
   autopsy, and reliability aggregation.
4. The reliability engine reads generated artifacts and writes JSON, TSV, HTML,
   and manifest summary outputs.
5. Evidence gates validate that supported public statements remain aligned with
   repository artifacts.

## Error Handling

Preflight blocked conditions:

- draft inferred AssaySpec under production profile;
- target audit unsafe at configured correction radius when
  `fail_on_unsafe_targets = true`;
- panel certificate missing or invalid when the assay declares a designed panel;
- unsupported backend forced by user configuration;
- missing target/read paths.

Runtime error conditions:

- native command failure;
- unreadable generated artifacts;
- malformed QC tables;
- GPU row mismatch or count delta when GPU validation is enabled.

Postrun reliability failures:

- assignment, ambiguous, unmatched, invalid, or correction rates crossing
  configured thresholds;
- autopsy identifying likely wrong offset or unsafe rescue;
- claim boundary missing for assay type when required.

Errors should not hide partial outputs. When assignment ran and produced
artifacts, reliability failure should mark the run as scientifically unsafe, not
delete the evidence needed to diagnose it.

## Evidence And Claims

The reliability report must preserve the existing evidence discipline:

- every public claim links to a gate, raw artifact, or generated report;
- unsupported scope is stated in the report;
- assay-specific comparator semantics are shown where available;
- ML and GPU outputs are advisory unless backed by deterministic validation.

New public claims require:

- raw evidence under `benchmarks/raw/`;
- generated report under `docs/benchmarks/`;
- gate script;
- test coverage for both passing and failing rows;
- update to `docs/scientific-claims.md`.

## ML And Statistics Scope

ML should initially support the human operator, not replace assignment.

Allowed first uses:

- anomaly scoring from reliability metrics;
- likely offset recommendation from unmatched-read structure;
- assay type suggestion during `assay infer`;
- natural-language explanation generated from deterministic findings.

Required guardrails:

- deterministic findings remain the source of truth;
- model output is labeled advisory;
- no learned assignment replaces exact target matching without a held-out public
  evidence gate and deterministic oracle comparison.

## Testing Strategy

Unit tests:

- AssaySpec `[reliability]` parsing and defaults;
- threshold evaluation;
- finding severity mapping;
- malformed artifact handling;
- backend eligibility decisions.

CLI tests:

- `assay check` writes preflight reliability artifacts;
- `assay plan` includes reliability plan output;
- `assay run` writes reliability JSON, TSV, HTML, and manifest summary;
- production profile fails unsafe targets;
- exploratory profile records unsafe targets without failing.

Evidence gate tests:

- reliability report contains supported and unsupported claim boundaries;
- public evidence rows are required before README-level claims are allowed;
- GPU backend cannot be promoted to production on synthetic evidence alone.

Regression fixtures:

- wrong offset;
- duplicate target;
- unsafe one-edit correction;
- high unmatched rate;
- high ambiguous rate;
- invalid extraction window;
- GPU unavailable;
- GPU mismatch.

## Rollout Plan

Phase 1: Reliability artifact MVP.

- Add `[reliability]` parsing.
- Aggregate existing audit, QC, and autopsy outputs.
- Generate JSON, TSV, and HTML report.
- Add tests and a release gate.

Phase 2: Pre-run assay certification.

- Integrate panel certificates and target safety more directly into
  `assay check`.
- Add production versus exploratory profiles.
- Improve recommendations for unsafe correction and wrong offset risk.

Phase 3: Backend intelligence.

- Record CPU/GPU eligibility in every run.
- Keep CPU as production authority.
- Add real-workload GPU gates for feature-barcode and BCL before any production
  promotion.

Phase 4: Advisory intelligence.

- Add deterministic anomaly scoring first.
- Add ML-assisted offset and failure explanations only after enough checked
  fixtures exist.

## Non-Goals

- replacing downstream screen statistics;
- Cell Ranger-compatible cell/UMI quantification;
- broad BCL/CBCL replacement claims;
- clinical variant interpretation;
- cloud-hosted data upload;
- black-box learned assignment;
- silent GPU backend switching.

## Design Decisions

1. `assay check` writes preflight reliability artifacts even when read-level QC
   is unavailable. Read-dependent fields are recorded as unavailable with a
   source-stage explanation.
2. `profile = "production"` fails fast before assignment on blocked preflight
   conditions. Runtime threshold violations are evaluated after assignment and
   leave generated outputs in place for diagnosis.
3. CLI reliability artifacts ship before Workbench reliability views. Workbench
   consumes the finalized JSON and TSV artifacts rather than creating a separate
   reliability model.
