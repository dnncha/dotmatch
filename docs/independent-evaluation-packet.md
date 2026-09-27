# Independent evaluation packet

Use this packet when an independent laboratory, core facility, workflow
maintainer, or bioinformatics team compares DotMatch with an established
workflow. It separates decisions made **before** results from observations made
afterward.

The packet does not send data anywhere. Keep it in the evaluator-controlled
workspace and replace private paths or names with approved references before
publishing any part of it.

The helper currently runs from a DotMatch source checkout and is not included in
the published 0.6.1 wheel. Record the source commit used for the packet as well
as the released DotMatch version used for the scientific comparison.

## 1. Create the protocol

```bash
python3 scripts/independent_evaluation_packet.py template evaluation.protocol.json
```

Replace every `REPLACE_ME`. The protocol requires:

- target-library revision and SHA-256;
- sample roles, read mates, and FASTQ SHA-256 values;
- immutable DotMatch and comparator versions, install routes, commands, and
  settings;
- primary endpoints with expected results and failure criteria;
- an explicit downstream endpoint, or an explicit statement that none is
  planned;
- interpretation limits.

For MAGeCK-style downstream comparisons, compare keyed scores, p-values, FDRs,
and declared hit sets. Raw ordinal rank or a top-N list is not a safe primary
endpoint when identical statistics can tie. Use
`compare_exact_unless_identical_statistical_tie` for ordinal ranks and either
`not_primary` or `report_boundary_ties` for top-N lists.

## 2. Lock it before viewing outcomes

```bash
python3 scripts/independent_evaluation_packet.py lock \
  evaluation.protocol.json evaluation.locked.json
python3 scripts/independent_evaluation_packet.py verify evaluation.locked.json
```

`lock` rejects placeholders, missing input identities, duplicate samples or
endpoints, undeclared downstream settings, and unsupported tie policies. It
writes a canonical SHA-256 over the evaluation ID and full protocol. Commit or
timestamp `evaluation.locked.json` in the evaluator's normal record system
before running either workflow. The hash detects later edits; the external
record establishes when the protocol existed.

## 3. Run both workflows unchanged

Use the locked commands and inputs. Preserve raw outputs privately. Do not tune
thresholds after seeing which result looks more favorable. If the protocol is
wrong, stop, explain why, lock a new evaluation ID, and retain the superseded
record.

## 4. Attach results

Create `evaluation.results.json` with this shape:

```json
{
  "schema_version": 1,
  "evaluation_id": "the-locked-evaluation-id",
  "protocol_sha256": "the-hash-from-evaluation.locked.json",
  "completed_at": "2026-09-27T12:00:00Z",
  "protocol_record_reference": "commit, DOI, timestamp receipt, or private record ID",
  "observations": [
    {
      "endpoint_id": "count_agreement",
      "observed": "state the measured result",
      "status": "pass"
    }
  ],
  "operational_metrics": [],
  "conclusion": "bounded conclusion owned by the evaluator",
  "limitations": ["software-output agreement is not biological correctness"],
  "public_use_permission": "state exactly what may be named or published"
}
```

Every locked primary endpoint must appear once. Allowed statuses are `pass`,
`fail`, `inconclusive`, and `not_run`; negative and incomplete results remain
visible.

```bash
python3 scripts/independent_evaluation_packet.py complete \
  evaluation.locked.json evaluation.results.json evaluation.completed.json
python3 scripts/independent_evaluation_packet.py verify evaluation.completed.json
```

Completion never changes the locked protocol or its hash. It rejects a result
with the wrong evaluation ID, hash, endpoint set, or unsupported status.

## What a completed packet establishes

A valid packet establishes that the recorded observations were attached to an
unchanged, internally complete protocol. It does not verify that the evaluator
actually ran the commands, that either workflow is biologically correct, that
the data are representative, or that public-use permission exists. Those
claims require independent evidence and approval.
