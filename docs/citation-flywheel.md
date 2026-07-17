# Citation Flywheel

This page tracks the highest-return work for making DotMatch easy to cite,
reuse, and recognize in external workflows. It is intentionally about adoption
mechanics, not scientific claims. Public claims still need the evidence gates
listed in `docs/scientific-claims.md`.

## Priority Actions

| Priority | Action | Done When |
| --- | --- | --- |
| P0 | Upstream official nf-core modules | A DotMatch module PR is merged into `nf-core/modules`, with tests, metadata, containers, and a workflow adoption record. |
| P0 | Publish a MultiQC module | MultiQC automatically discovers DotMatch outputs without a custom config, renders core assignment metrics, and has a released package or upstream integration. |
| P1 | Publish Galaxy and Snakemake wrappers | Reviewed wrappers are available outside this repository, with tiny fixtures and documented outputs. |
| P1 | Generate run citation artifacts | Done in AssaySpec runs: each production run writes `methods.md`, `CITATION.bib`, and `software_versions.yml` beside summary and report outputs. |
| P1 | Add copyable methods and citation UI | Workbench and generated HTML reports expose a single "Copy methods and citation" action backed by the generated methods artifact. |
| P1 | Document external evaluations | Three to five independent labs, cores, or workflow maintainers approve public use records with scope notes. |
| P1 | Maintain the workflow integration kit | Homepage audience routes, workflow handoffs, citation links, and public-use record rules stay aligned in `docs/workflow-integration-kit.md`. |
| P1 | Maintain the workflow integration roadmap | `docs/workflow-integration-roadmap.md`, `docs/workflow-integration-plan.json`, and the homepage evaluation surface stay aligned. |
| P2 | Publish a methods and benchmark preprint | The manuscript demonstrates speed and the consequences of ambiguity, unsafe rescue, and incorrect barcode windows. |
| P2 | Maintain public use records | Approved independent examples are listed in `docs/adopters/` and mirrored in the public docs. |

## Operating Rules

- Record external integrations in `docs/workflow-adoption.json` only after they
  are accepted, released, or published outside this repository.
- Add public use records only with approved wording and a public URL.
- Keep citation text generated from checked project metadata, not manually
  copied into each report surface.
- Treat ambiguity, unsafe rescue, shifted windows, and barcode-window mistakes
  as first-class methods details whenever they affect interpretation.

## Existing Owners

- External workflow submission checklist:
  [Workflow Submission Pack](workflow-submissions.md).
- Methods and citation language:
  [Methods and Citation Template](methods-and-citation.md).
- Public use records:
  [DotMatch Public Use Records](adopters/README.md).
- Evidence limits for public claims:
  [DotMatch Evidence Notes](scientific-claims.md).

## Implementation Notes

The next implementation slice should avoid one-off text in reports. Prefer a
small shared artifact writer that can be called by CLI commands, Workbench
exports, workflow wrappers, and report generators. A complete run should emit:

- `methods.md`: copyable methods prose with DotMatch version, command, input
  semantics, edit-distance metric, radius, ambiguity policy, alphabet policy,
  extraction window, and relevant warnings.
- `CITATION.bib`: BibTeX generated from `CITATION.cff`, release DOI, and package
  version metadata.
- `software_versions.yml`: versions for DotMatch, Python package, native
  library, workflow wrapper, and any comparator or report tool used in the run.

The Workbench and HTML reports should copy the contents of `methods.md` plus the
canonical citation entry instead of maintaining their own text templates.
