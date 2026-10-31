# Sequence Doctor: Research Evidence Passport pilot

A source-only, local-first pre-execution trust gate for known-target FASTQ projects.
This is a review candidate, not a released DotMatch feature or a biological validator.
The counting engine and published CLI are unchanged.

## Run the complete demonstration

From a source checkout on Linux or macOS, with Python 3.10 or later:

```sh
python -m unittest discover -s python/tests -p test_evidence_passport.py -v
python scripts/evidence_passport_demo.py --out /tmp/sequence-doctor-demo
```

Choose an output directory that does not already exist. The example uses one
synthetic FASTQ and a deliberately incorrect reference. It starts a real stdio
MCP subprocess; no LLM account, native DotMatch build, network or paid compute is
required. The scripted client discovers the tools, observes a refused reference,
tries and fails to compile, observes the restored expected reference, requests
human review, and still cannot compile before the separate approval ceremony.
It then obtains an authorized compilation input and exports the passport.
The demonstration uses a clearly fictional reviewer, not a real human attestation.

The output includes `research-evidence-passport.json`, its payload digest,
`refused-evidence-passport.json`, `workflow-compilation-input.json`,
`mcp-transcript.json`, a deliberately rehashed alteration, and `demo-result.json`.
The final passport retains the earlier refusal and both inspections.

## Inspect a passport independently

```sh
python python/assaycode/evidence_passport.py verify \
  --passport /tmp/sequence-doctor-demo/research-evidence-passport.json \
  --expected-sha256 YOUR_INDEPENDENTLY_RETAINED_PAYLOAD_DIGEST
```

The payload digest is SHA-256 of the restricted canonical JSON payload, NOT a
checksum of the pretty-printed file. Without an independently retained digest,
verification establishes internal consistency only. Keeping the digest beside
the passport is convenient but is not independent anchoring. Hashing is not a
signature, human authentication, biological validation, or proof of completeness.

## Connect a local MCP client

Set the client's stdio command to `python`, with arguments:

```json
["/absolute/path/to/evidence_passport.py", "serve", "--root", "/absolute/path/to/project"]
```

The process root and `project.json` are selected by the operator, not by tool
arguments. The five read-only tools are `inspect_project`,
`draft_analysis_contract`, `check_analysis_contract`, `compile_workflow`, and
`export_evidence_passport`. The protocol version is 2025-11-25. A real scripted
stdio exchange is tested; interoperability with a named third-party client has
not been certified. There is no approval, filesystem-write, shell or execution
tool. An agent with separate OS/shell access is outside this server's boundary.

`compile_workflow` returns an engine-neutral **compilation input**; an external
adapter still has to compile and execute a real workflow. It does not currently
run Nextflow, Snakemake, or a biological analysis.

## Actual human review, outside the agent tool surface

Use the generated synthetic `project/project.json` as a format example. Supply
independently established expected file hashes; the tool cannot decide whether a
biologically wrong but self-consistent reference declaration is correct.

```sh
python python/assaycode/evidence_passport.py draft --root PROJECT > contract-draft.json
# Review the project and all six decision values and copy the draft into PROJECT:
cp contract-draft.json PROJECT/analysis-contract.json
python python/assaycode/evidence_passport.py approve --root PROJECT \
  --reviewer 'YOUR REAL REVIEWER NAME' \
  --rationale 'YOUR ACTUAL REVIEW RATIONALE' \
  --confirm-all-reviewed > contract-reviewed.json
mv contract-reviewed.json PROJECT/analysis-contract.json
python python/assaycode/evidence_passport.py compile --root PROJECT
```

Never redirect approval output onto its own input file: the shell would truncate
it before the program can read it. Approval is currently a named, self-asserted
local-operator record, not an authenticated identity. For finer-grained rationale,
edit the individual decisions and rebuild the approval through a reviewed human
procedure; changing any decision after approval invalidates its binding. A
blocked inspection cannot be overridden by naming a reviewer.

See [TECHNICAL.md](TECHNICAL.md) for the threat model and validation boundaries.
The schema is at `docs/schemas/research-evidence-passport-v1.schema.json`.

## Scope discipline

The immediate delivery is this narrow source demo, export, offline verifier and
tests. Authenticated organizational approvals, an independently anchored signature
service, workflow-engine adapters, field pilots, and corroboration policy are
not silently claimed as shipped. HMMForge and EditWitness are not integrated.
