# Commercial Boundary

DotMatch uses an open-core strategy. The public repository should stay useful,
auditable, and scientifically credible on its own. Commercial work should add
team, provenance, support, and evidence-management capabilities without
crippling the open-source engine or changing its deterministic assignment
semantics.

The open-source DotMatch project includes:

- the deterministic assignment engine and CLI;
- local CRISPR guide counting, barcode demultiplexing, panel checking, and
  fixed-window workflows;
- public report generation and machine-readable outputs;
- schemas, examples, tutorials, benchmark reports, and evidence boundaries;
- local workflow integration examples for reproducible pipelines.

DotMatch Pro is the commercial assay reliability workbench for teams that need:

- hosted or team workspaces;
- run registries and searchable run history;
- signed reports and reviewer-ready evidence packets;
- private assay registries and private assay specifications;
- enterprise connectors for storage, identity, workflow systems, and LIMS-like
  handoffs;
- commercial support, onboarding, validation review, and audit packs.

The open-source engine remains available under Apache-2.0. DotMatch Pro may
package or call the open engine, but it must not silently force assignments,
hide `ambiguous` or `none` outcomes, or weaken the public contract that read
states are reported as `unique`, `ambiguous`, `none`, or `invalid`.

Commercial features should strengthen reliability, provenance, auditability, or
evidence-backed claims. They should not add clinical, diagnostic, FDA, CE-IVD,
patient-care, or treatment-decision claims unless the project maintainer
explicitly changes the documented scope.

## Repository Readiness Checklist

Use this checklist before publishing public-repo changes that affect release,
governance, documentation, or commercial boundary wording:

- license present and still Apache-2.0;
- security policy present;
- no raw customer assay data, private FASTQ/BAM/BCL data, or unminimized
  biological inputs committed;
- docs build passes with `make docs-ready`;
- tests pass with the relevant native, CLI, and Python commands for the change.
