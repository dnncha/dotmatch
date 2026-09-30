# Technical specification and trust boundaries

## Data flow and failure states

The operator selects a local project root and a manifest containing a biological
question, sample-to-file mapping, versioned reference, analysis method, study
design, reporting limitations, and expected SHA-256 digests. Sequence Doctor
hashes complete file bytes and checks a bounded FASTQ prefix. A mismatch,
malformed prefix, missing/unsafe file or exceeded budget produces refusal.
Matching bytes produce `needs_human_review`, never automatic approval.

The Analysis Contract has six decisions: biological question, sample identity,
reference choice, study design, analysis method and reporting limits. Its approval
binds the exact inspection and all decision values and rationales. A fresh file
inspection and reread of the saved contract precede authorization. Changed files,
manifest, decisions or approvals invalidate readiness. A new inspection never
erases an older refusal or approval from the session history.

The current gate validates declarations and their review, not the biological
truth of the declarations. For example, it detects a reference different from an
independently expected hash; it does not infer that a matching reference is the
right genome or guide library. Similarly, it requires sample/replication metadata
but does not independently identify donors or adjudicate biological replication.

## Passport and verification

The versioned JSON envelope contains `payload` and `integrity`. The payload binds
project snapshots, declared and observed file hashes, reference versions, bounded
QC, all six decisions, approval records, tool source fingerprint, self-reported
client identity, compilation inputs and a hash-linked event history. It contains
no biological results because this is a pre-execution pilot. Events retain the
complete approved contract/compilation snapshot, so later replacement of the
current contract does not destroy the older decision's reconstructibility.

Each event contains a contiguous index, UTC timestamp, kind, subject hash,
details, previous-event hash and its own hash. The verifier checks these links,
inspection order, report/contract/plan bindings, required prior approval, current
readiness, version identifiers and explicit limitations. The JSON Schema checks
structure; the Python semantic verifier is required for cross-field semantics,
strict parsing, hashes and state transitions. Timestamps are observations from
the local clock, not independently certified time. A producer source fingerprint
is a byte identity, not proof that trusted software produced the document.

`dotmatch-json-v1` is UTF-8 JSON with lexicographically sorted Python string keys,
compact separators, preserved Unicode, no floats, no non-finite values and
integers restricted to the portable exact range. Duplicate keys are rejected on
input. This is a specifically named restricted format, **not RFC 8785/JCS**.
Non-Python implementations must reproduce its sorting and escaping rather than
silently assuming another canonicalization scheme.

An external expected payload digest detects modification even after an attacker
recalculates every internal hash. Without that independent value, an attacker
can rewrite the entire envelope. The returned assurance is correspondingly
`internally_consistent_only` or `verified_against_independent_digest`. The
verifier explicitly reports that identities are unauthenticated, biological
validity is not established and external file bytes were not rechecked. A
separate consumer must rehash mounted inputs against the passport at actual
execution time; this pilot cannot prevent changes after its compilation gate.

## Operational limits

Python 3.10+ and POSIX no-follow descriptor operations are required. There are
at most 128 declared files and samples, a 1 GiB per-file limit, a 2 GiB accepted
project-size limit, 1 MiB FASTQ lines, an 8 MiB decompressed inspection-prefix
budget and a 1,000-record inspection cap. FASTQ syntax beyond that prefix and
some compressed-stream tail errors are not checked. Full raw/compressed file
bytes are hashed, but that is not complete FASTQ syntax validation.

Passport/JSON input is capped at 2 MiB, nesting at 64 levels, session history at
256 reports and 4,096 events. These are bounded pilot limits, not a suitable
production-scale FASTQ envelope. The accepted aggregate size does not constitute
an exact global I/O-time budget. Reinspection deliberately rereads bytes rather
than relying on path/mtime caches. No production throughput claim is made.

## Threat model

The server rejects root traversal and symlinked path components and reads only
regular files using directory descriptors. It does not accept shell commands,
network destinations or file writes via MCP. This is not a sandbox against the
same OS user, hard-linked external data, independently granted shell/file tools,
compromised Python/runtime, falsified expected hashes or a malicious reviewer.
The operator is responsible for a dedicated project root, least-privilege OS
access, trust in the original manifest and an appropriate approval channel.

Named approval is a self-asserted local record. The separate CLI is a workflow
boundary, not authentication: an agent with another write path can forge such a
record. Production organizational enforcement requires an authenticated signer
and a policy authority outside the agent's permissions. Neither a name nor a
hash chain establishes nonrepudiation. The demo's reviewer is explicitly fictional.

Only events observed in the current server session are included; this is not a
complete record of all project history, external commands, edits or identities.
Raw sequences and absolute paths are not copied by the exporter, but user-entered
question/rationale/sample text can itself contain identifying information. Users
must review those fields before sharing. This is not an automatic anonymizer.

## Interoperability and research limits

The stdio MCP implementation is exercised through an actual scripted subprocess,
with initialize, initialized notification, discovery and calls. It emits only
newline JSON-RPC on stdout and refuses tool-side arguments. No live LLM or named
external client has been validated. `compile_workflow` produces a typed handoff
with `external_adapter_required` and `not_executed`, not a Nextflow/Snakemake run.

This schema is not a claimed Workflow Run RO-Crate, BioCompute, Sigstore or MCP
conformance certification. Post-award mapping/adapters should reuse established
standards. It is research-use-only and makes no clinical/diagnostic, regulatory
compliance, biological accuracy, paid-pilot or user-adoption claim.
