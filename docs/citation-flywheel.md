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
| P1 | Recruit documented external pilots | Three to five independent labs, cores, or workflow maintainers approve public pilot records with scope notes. |
| P1 | Maintain the outreach and integration kit | Homepage audience routes, outreach copy, workflow handoffs, citation links, and adopter-record rules stay aligned in `docs/industry-exposure.md`. |
| P1 | Execute the next 10 distribution actions | `docs/industry-next-wins.md`, `docs/industry-exposure-plan.json`, and the homepage action section stay aligned. |
| P2 | Publish a methods and benchmark preprint | The manuscript demonstrates speed and the consequences of ambiguity, unsafe rescue, and incorrect barcode windows. |
| P2 | Maintain a public used-by page | Approved independent examples are listed in `docs/adopters/` and mirrored in the public docs. |

## Operating Rules

- Do not count private conversations, internal examples, or unmerged PRs as an
  accepted external-use record.
- Record external integrations in `docs/workflow-adoption.json` only after they
  are accepted, released, or published outside this repository.
- Keep citation text generated from checked project metadata, not manually
  copied into each report surface.
- Treat ambiguity, unsafe rescue, shifted windows, and barcode-window mistakes
  as first-class methods details whenever they affect interpretation.

## Existing Owners

- External workflow submission checklist:
  [Workflow Submission Pack](workflow-submissions.md).
- Methods and citation language:
  [Methods and Citation Template](methods-and-citation.md).
- Approved external-use records:
  [DotMatch Adopter Notes](adopters/README.md).
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
