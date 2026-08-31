# Evidence policy

Treat every output as a bounded local evidence record.

- `passed` means the configured deterministic and reliability gates passed for the declared inputs. It is not proof of biological validity or scientific adoption.
- `needs_review` means counts or preflight evidence exist but require explicit review before interpretation.
- `blocked`, `failed`, `invalid_input`, and `interrupted` are not production-ready outcomes.
- Record the exact DotMatch version, tool-contract version, candidate spec revision and SHA-256, artifact SHA-256 values, reliability status, and findings.
- Preserve all ambiguous assignments as ambiguous. Never turn them into target counts.
- Handoff bundles exclude raw FASTQ. They contain hashes of declared inputs so a reviewer in the originating controlled workspace can verify lineage.
- CPU assignment is authoritative unless the workload-specific GPU validation described by the reliability record passes.
- A fixture, local run, package check, merged upstream contribution, released upstream integration, public deployment, download, and approved scientific use are different proof classes. Describe only the class actually evidenced.
