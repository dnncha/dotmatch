# Bioinformatics Evaluation Packet

This packet is for bioinformatics teams, core facilities, workflow maintainers,
and assay-development groups evaluating whether DotMatch is appropriate for a
known-target sequencing workflow. It records what is available now, what is
validated by repository checks, and what should not be inferred from the current
release.

DotMatch should be evaluated as a local command-line and Python package for
deterministic assignment of fixed read windows to known short DNA targets. It is
not a genome aligner, basecaller, variant caller, adapter trimmer, cell/UMI
quantifier, or screen-level statistical analysis package.

## Current Package Surface

For release `0.3.0`, the package surface is:

| Surface | Current state | Evidence or source |
| --- | --- | --- |
| PyPI | Version 0.3.0 is published; release CI smoke-tested the sdist and platform wheels | `docs/distribution-release.json` |
| Bioconda | Version 0.2.2 remains published for linux-64, linux-aarch64, osx-64, and osx-arm64; 0.3.0 is not published there yet | `docs/distribution-release.json` |
| GHCR container | Public v0.3.0 multi-architecture manifest is verified; release CI smoke-tested the amd64 and arm64 native CLI paths | `docs/distribution-release.json` |
| BioContainers | The 0.2.2 images remain public on Quay; 0.3.0 depends on Bioconda propagation | `docs/distribution-release.json` |
| Documentation | Sphinx docs and public schemas in repository | `docs/index.md`, `docs/schemas.md` |
| Citation | `CITATION.cff`, Zenodo concept DOI, and generated run artifacts | `docs/methods-and-citation.md` |
| Python API | `dotmatch` package, streaming helpers, pandas/polars/AnnData interop, and a MultiQC parser entry point; 0.3.0 includes automatic MultiQC discovery | `pyproject.toml`, `docs/streaming-api.md`, `docs/ecosystem-status.md` |
| R interface | Reticulate-backed package skeleton and vignette | `R/`, `vignettes/dotmatch.Rmd` |

Do not describe a channel as verified until `make distribution-channels` passes
for that channel and the release record is updated.

## Validated Assay Scope

Use `docs/assay-evidence.json` and `docs/scientific-claims.md` as the source of
truth. A short evaluator summary:

| Assay or workflow area | Current public status | Boundary |
| --- | --- | --- |
| CRISPR guide counting | Supported for checked public lanes and comparator semantics | Not screen-level effect analysis |
| Fixed-position inline barcode demultiplexing | Supported for checked SRP009896 exact-prefix and fixed-length Hamming lanes | Not arbitrary adapter trimming or raw BCL replacement |
| Feature-barcode assignment | Supported for checked per-read 10x TotalSeq-B fixed-window assignment | Not Cell Ranger feature-matrix parity |
| Perturb-seq guide capture | Gated fixed-window per-read evidence | Not guide-per-cell quantification or perturbation-effect inference |
| Amplicon or panel primer-start assignment | Supported for checked fixed-window ARTIC primer-start lane | Not consensus generation, variant calling, or clinical interpretation |
| Oligo or adapter-prefix assignment | Supported for checked fixed-window adapter-prefix assignment | Not adapter trimming, read merging, or UMI grouping |
| Classic BCL milestone | Narrow checked tiny-BCL milestone | Not production Illumina demultiplexing or CBCL/NovaSeq support |

If a statement would affect a procurement decision, manuscript claim, pipeline
replacement, or regulated workflow, link it to the relevant raw artifact,
benchmark report, and gate.

## Minimum Local Evaluation

Start from the released package, not an unpublished checkout, unless the
evaluation is explicitly for development work:

```bash
python3 -m pip install dotmatch==0.3.0
dotmatch --version
dotmatch dist ACGT AGGT
```

For Conda-based environments:

```bash
conda create -n dotmatch -c conda-forge -c bioconda dotmatch=0.2.2
conda activate dotmatch
dotmatch --version
```

For a first source-level review:

```bash
make test
make cli-test
make python-test
make workflow-examples-ready
make repository-ready
```

For scope-sensitive review, run the specific evidence gates named in
`docs/methods-and-citation.md` and `docs/scientific-claims.md`; do not infer
broader assay coverage from a generic test pass.

## Outputs To Inspect

A useful evaluation should inspect the files a lab or workflow system would
actually receive:

- `sample_qc.tsv`: assignment rate, ambiguous reads, unmatched reads, invalid
  windows, rescue counts, and representation fields.
- `summary.json`: command, version, policy, target-window configuration, and
  provenance.
- `assignments.tsv`: per-read states when requested.
- `top_unmatched.tsv`: recurring unmatched sequences for assay review.
- `assay_report.html` or workflow report HTML: human-readable reliability
  review.
- `methods.md`, `CITATION.bib`, and `software_versions.yml`: methods and
  software citation artifacts produced by AssaySpec run paths.
- `assay_manifest.summary.tsv`: workflow-facing contract for run artifacts.

The central contract is that `ambiguous`, `none`, and `invalid` outcomes remain
visible instead of being silently collapsed into target counts.

## Workflow Integration Status

`docs/workflow-adoption.json` is currently `not_ready` because no external
workflow-manager integration has been accepted yet. This is separate from
Spack package-recipe acceptance and from PyPI or Bioconda distribution.

Public submissions are open in nf-core/modules, Galaxy IUC, Snakemake wrappers,
and MultiQC. Their exact heads, check state, and acceptance gates are recorded
in `docs/ecosystem-status.md` and `docs/workflow-submissions.md`. None of those
open pull requests is described as accepted or released. The bio.tools metadata
remains an unsubmitted local draft, and no exact DotMatch workflow was found in
WorkflowHub on 2026-08-28.

Record accepted integrations only in `docs/workflow-adoption.json`. Public use
records require approved wording and a public URL.

## Language Rules For Public Surfaces

Use plain technical language:

- say "known-target assignment", not broad "sequencing analysis";
- say "checked public lane", "gated", or "experimental" when evidence is
  narrow;
- say "workflow artifact", "schema", "gate", "report", and "command" when
  there is a concrete file or check;
- avoid broad replacement claims, superlatives, and launch copy;
- do not imply external workflow integration before
  `docs/workflow-adoption.json` records an accepted public integration.

The project should feel like a bioinformatics package first: installable,
auditable, scoped, reproducible, and cautious about unsupported claims.
