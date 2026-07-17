# Reviewer Evidence Matrix

This page maps common JOSS-style reviewer concerns to concrete DotMatch
artifacts. It is intended to be pasted into a resubmission PR or release notes
and kept current as evidence changes.

| Reviewer concern | Evidence artifact | How to verify |
| --- | --- | --- |
| Runnable example | Tiny deterministic quickstart and reproducibility harness | Run `make repro`; inspect `repro/small/repro_manifest.json` |
| Tests and CI | Cross-platform CI, Python compatibility checks, native C tests, package smoke tests, evidence gates | See `.github/workflows/ci.yml`; run `make python-test` and `make repository-ready` locally |
| Packaging | PyPI metadata, Bioconda recipe, wheel build configuration, versioned CLI | See `pyproject.toml`, `packaging/bioconda/meta.yaml`, `.github/workflows/release.yml`, and `dotmatch --version` |
| Reproducibility | Deterministic tiny fixture with stable counts, assignments, summary, and checksums | Run `python3 scripts/repro_small.py --out-dir repro/small` |
| Scientific claims | Scoped evidence notes and machine-checked benchmark/readiness gates | See `docs/scientific-claims.md`, `docs/benchmarks/README.md`, and `docs/evidence-gallery/README.md` |
| Public CRISPR example | Public Yusa/CRISPR guide-counting fixture and comparison reports | See `examples/crispr_guides/`, `examples/crispr_sanson_brunello/`, and `docs/benchmarks/public_crispr/README.md` |
| Workflow reuse | nf-core-style, Snakemake, MultiQC, and Galaxy integration examples | See `examples/workflows/` and `docs/workflow-adoption.json` |
| Citation and archival record | CITATION metadata and Zenodo DOI badge | See `CITATION.cff` and `docs/methods-and-citation.md` |
| Public use records | External confirmations | Track approved public confirmations in `docs/adopters/README.md` |

## Minimal Reproduction Command

```bash
make repro
```

The command writes all generated files under `repro/small/`, including
`README.md`, `repro_manifest.json`, `reviewer-repro-packet.json`,
`resubmission-matrix.csv`, and one log tail per focused command. The generated
manifest records the git commit, dirty-worktree status, command outcomes, and
SHA-256 checksums for the evidence files reviewers should inspect first.

CI also runs:

```bash
make repro-small
```

The `ci` workflow uploads `repro/small/` as the `reviewer-repro-packet`
artifact from a dedicated Linux job.

## Resubmission Matrix

`repro/small/resubmission-matrix.csv` is generated from
`docs/assay-evidence.json`. Each row ties an assay lane to its current support
status, gate command, validated scope, and next public-evidence item.
The compact target runs native C tests, CLI fixture tests, scientific-readiness
guardrails, evidence-gallery validation, and assay-evidence validation. Full
public-data regeneration and comparator-backed gates remain in the
benchmark-specific `make bench-*` and `make *-gate` targets because those runs
can require larger downloads or external comparators.

Workflow adoption is tracked separately in `docs/workflow-adoption.json`.
Current in-repository Galaxy, Nextflow, nf-core-style, Snakemake, and MultiQC
examples are workflow-example evidence, not accepted external integrations. The
workflow manifest remains `not_ready` until an external workflow record is
accepted, released, or published.

## Resubmission Checklist

- Refresh CI and release badges in `README.md`.
- Run `make repro`, `make repository-ready`, and the relevant evidence gates for
  any claim being made in the submission text.
- Link to release artifacts: PyPI, Bioconda, GitHub Release, and Zenodo DOI.
- Add only approved public use records to `docs/adopters/README.md`.
- Keep claims aligned with `docs/scientific-claims.md`; move aspirational or
  unverified statements to roadmap text.
