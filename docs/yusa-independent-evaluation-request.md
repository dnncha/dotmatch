# Draft independent Yusa evaluation request

This is an **unsent draft**, not an adoption record. Choose and verify a
recipient separately before sending it. The recipient must control protocol
locking, execution, interpretation, and permission to publish.

## Draft message

**Subject:** Independent check of DotMatch 0.6.1 on the public Yusa workflow

Hello,

Would you be willing to independently compare released DotMatch 0.6.1 with
MAGeCK 0.5.9.5 on the public Yusa ERR376998/ERR376999 workflow?

This is a request for a disconfirming reproduction, not an endorsement. The
inputs are public and already have fixed SHA-256 identities. Before viewing
outcomes, you would replace the two evaluator fields in the linked draft,
review or change every command and endpoint, and lock it in your own record
system. Please stop if an input hash differs.

The primary checks are exact identifier-keyed guide/gene/count agreement and,
if you run `mageck test`, keyed non-rank statistics plus FDR≤0.05 hit sets.
Ordinal ranks may differ only inside blocks with identical score, p-value, and
FDR; an arbitrary top-N boundary is not a primary result.

You own the conclusion, including a negative or inconclusive result. Raw
outputs can remain private. Nothing will identify you or your organization
without approval of the exact wording.

The scientific run remains frozen to released DotMatch 0.6.1. The packet helper
is included in DotMatch 0.6.2; if it is used only to lock the record, record
that package version and artifact identity separately.

Would this bounded public reproduction fit your review process?

## Ready-to-review materials

- [Draft protocol JSON](../examples/evaluation_packets/yusa-independent-protocol.draft.json)
- [Independent packet instructions](independent-evaluation-packet.md)
- [Existing Yusa evidence and limitations](benchmarks/mageck_yusa_reference/README.md)

The draft deliberately contains `REPLACE_ME` in `evaluation_id` and
`protocol_owner`, so `dotmatch evaluation-packet lock` refuses it. The
evaluator must replace those fields, review every prefilled identity, command,
endpoint, expectation, and failure rule, and then create the prospective lock.

The prefilled public identities are:

| Input | SHA-256 |
| --- | --- |
| ERR376998 FASTQ | `7f79b76cec12b70319744417282f963c00818a5f0ae61497bd7b64790ac55f2f` |
| ERR376999 FASTQ | `cf2bc10938e178d16dfb81ca2f9fda805cae892290fac3cfc243bb637f8cda17` |
| Yusa library | `d41a4122a46fdedb47deca61081dde9037d461a1fce1e6f9744522e21005ceba` |

These values are input identities, not expected-output evidence. DotMatch's
prior project-run outputs must not serve as the evaluator's oracle.

## Permission choices

The evaluator should choose one result-record boundary:

1. private result only;
2. anonymized workflow and aggregate outcome;
3. named evaluator or organization with approved wording and evidence link.

Silence, a download, repository traffic, a project-authored rerun, or an unsent
draft does not count as independent use.
