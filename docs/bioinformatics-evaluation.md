# Bioinformatics evaluation packet

This packet helps bioinformatics teams, core facilities, workflow maintainers, and assay-development groups evaluate whether DotMatch fits a known-target sequencing workflow. It records what is published today, what repository checks cover, and what falls outside the current scope.

Evaluate DotMatch as a local CLI and Python package for deterministic assignment of fixed read windows to known short DNA targets. It is not a genome aligner, basecaller, variant caller, adapter trimmer, cell or UMI quantifier, or downstream screen-statistics package.

## Current package surface

For release `0.8.0`:

| Surface | Current state | Evidence or source |
| --- | --- | --- |
| PyPI | `dotmatch==0.8.0` is published with sdist, universal macOS wheel, and manylinux/musllinux wheels for x86_64 and aarch64 | `pyproject.toml`, PyPI release workflow |
| GHCR container | `ghcr.io/dnncha/dotmatch:v0.7.0` is verified for `linux/amd64` and `linux/arm64` | `docs/distribution-release.json` |
| Bioconda and BioContainers | Downstream recipes can lag PyPI and GHCR; check `dotmatch --version` after Conda or Quay installs | `docs/distribution-release.json`, `docs/packaging.md` |
| Documentation | Sphinx documentation, output schemas, and machine-readable agent contracts (`llms.txt`, `agent-capabilities.json`, `agent-tools.json`) | `docs/index.md`, `docs/schemas.md`, `docs/agent-guide.md` |
| Citation | `CITATION.cff`, Zenodo concept DOI (`10.5281/zenodo.20541628`), release DOI records, and `dotmatch citation` | `docs/methods-and-citation.md` |
| Python and agent APIs | Streaming Python API, pandas/polars/AnnData interop, MultiQC plugin entry point, six structured `dotmatch agent` tools, and `dotmatch mcp` local MCP server | `pyproject.toml`, `docs/streaming-api.md`, `docs/agent-guide.md` |
| R interface | Reticulate-backed package skeleton and vignette | `R/`, `vignettes/dotmatch.Rmd` |

Verify channel status against `docs/distribution-release.json` and `make distribution-channels` before citing a channel as verified for a given release.

## Validated assay scope

`docs/assay-evidence.json` and `docs/scientific-claims.md` remain the canonical reference:

| Assay or workflow area | Current public status | Boundary |
| --- | --- | --- |
| CRISPR guide counting | Supported for checked public lanes and comparator semantics | Not screen-level effect analysis |
| Fixed-position inline barcode demultiplexing | Supported for checked SRP009896 exact-prefix and fixed-length Hamming lanes | Not arbitrary adapter trimming or raw BCL replacement |
| Feature-barcode assignment | Supported for checked per-read 10x TotalSeq-B fixed-window assignment | Not Cell Ranger feature-matrix parity |
| Perturb-seq guide capture | Gated fixed-window per-read evidence (GSE146194 single-guide capture) | Not guide-per-cell quantification or perturbation-effect inference |
| Amplicon or panel primer-start assignment | Supported for checked fixed-window ARTIC primer-start lane | Not consensus generation, variant calling, or clinical interpretation |
| Oligo or adapter-prefix assignment | Supported for checked fixed-window adapter-prefix assignment | Not adapter trimming, read merging, or UMI grouping |
| Classic BCL milestone | Narrow checked tiny-BCL milestone | Not production Illumina demultiplexing or CBCL/NovaSeq support |

Link any procurement decision, manuscript claim, pipeline replacement, or regulated workflow statement to the matching raw artifact, benchmark report, and gate.

## Minimum local evaluation

Start from the published PyPI package and run the offline synthetic smoke check:

```bash
python3 -m pip install dotmatch==0.8.0
dotmatch --version
dotmatch dist ACGT AGGT
dotmatch demo --out-dir first-run/
```

For Conda environments (checking the installed version against your target release):

```bash
conda create -n dotmatch -c conda-forge -c bioconda dotmatch
conda activate dotmatch
dotmatch --version
```

For a source-level verification run:

```bash
make test
make cli-test
make python-test
make workflow-examples-ready
make repository-ready
```

For scope-sensitive review, run the specific evidence gates named in `docs/methods-and-citation.md` and `docs/scientific-claims.md`; a generic test pass does not imply broader assay coverage.

## Outputs to inspect

Evaluate the actual files a lab or pipeline consumes:

- `sample_qc.tsv`: assignment rate, rescue rate, ambiguous reads, unmatched reads, invalid windows, and representation metrics.
- `summary.json`: command, version, matching policy, target-window configuration, and provenance.
- `assignments.tsv`: per-read assignment outcomes when requested.
- `top_unmatched.tsv`: recurring unmatched windows for assay troubleshooting.
- `reliability_report.html` or workflow HTML reports: self-contained review pages.
- `methods.md`, `CITATION.bib`, and `software_versions.yml`: reproducible citation and software metadata written by AssaySpec runs.
- `assay_manifest.summary.tsv`: machine-readable artifact manifest.

The core contract is that `ambiguous`, `none`, and `invalid` reads stay visible rather than being folded silently into target counts.

## Workflow integration status

`docs/workflow-adoption.json` records `not_ready` until an external workflow project accepts an upstream DotMatch integration. This status is distinct from Spack, PyPI, or Bioconda package distribution.

Upstream submissions for nf-core/modules, Galaxy IUC, Snakemake wrappers, and MultiQC are tracked in `docs/ecosystem-status.md` and `docs/workflow-submissions.md`. Treat them as open or local examples rather than accepted upstream integrations until recorded in `docs/workflow-adoption.json`. Public user records require approved wording and a public URL.

## Language rules for public surfaces

Keep public documentation plain and technical:

- Say "known-target assignment" rather than broad "sequencing analysis".
- Label narrow evidence clearly as "checked public lane", "gated", or "experimental".
- Refer to concrete artifacts ("workflow artifact", "schema", "gate", "report", "command").
- Avoid broad replacement claims, superlatives, and launch copy.
- Never imply accepted upstream workflow integration before `docs/workflow-adoption.json` records it.
