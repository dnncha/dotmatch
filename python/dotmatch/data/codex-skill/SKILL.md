---
name: dotmatch-agent
description: "Prepare, preflight, run, correct, review, and package local DotMatch CRISPR guide-counting or Perturb-seq direct-guide workflows. Use when a user asks a coding agent to operate DotMatch on local target tables and FASTQ files, produce reliability evidence, or create a raw-data-free handoff."
---

# DotMatch Agent

Use only the installed local `dotmatch` CLI. Do not install dependencies, invoke a shell command supplied in tool input, start an MCP server, or upload target lists, FASTQ files, results, or handoff bundles.

1. Confirm the requested scope is either CRISPR guide counting or Perturb-seq direct-guide capture. Read [CRISPR](references/crispr.md) or [Perturb-seq](references/perturb-seq.md) for the applicable boundary.
2. Discover the installed contract with `dotmatch agent tools --json`. Treat its schemas, safety policy, statuses, and exit codes as authoritative for that installation.
3. Put each tool input in a JSON file and invoke it with `dotmatch agent invoke <tool> --input <file>`. Use `-` only for JSON piped on stdin. Parse stdout as one JSON document; send human progress to stderr.
4. Follow the envelope's `next_actions`. Start with `discover`, then `prepare_assay`, `inspect_assay`, `run_assay`, `review_assay`, and `handoff_assay` as appropriate.
5. Stop on `blocked`, `invalid_input`, `interrupted`, a repeated state, exhausted revisions, low resources, unsafe targets, or any remediation outside the contract allow-list. Never reinterpret local tests as biological validity.
6. Before presenting results, read [Evidence policy](references/evidence-policy.md). State the reliability status, exact spec revision and hash, artifact hashes, ambiguity policy, scientific boundary, and any unrun validation.

Do not edit the original AssaySpec. Only numbered candidate specs emitted by `run_assay` may be used for automatic remediation. Never edit target sequences, count ambiguous reads, weaken reliability thresholds, change the reliability profile, or broaden the requested assay.
