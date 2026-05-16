# Release Process

DotMatch releases should be specific, reproducible, and evidence-bounded.

## Pre-Tag Checks

Run the consolidated local pre-tag gate:

```bash
make pretag-ready
```

It is a local readiness gate and intentionally does not include
`make distribution-channels`, `make workflow-adoption-status`, or
`make bcl-comparison-gate`. Those gates require public/external evidence and
are listed below.

The target runs:

```bash
make test
make cli-test
make python-test
make python-package-test
make repository-ready
make release-ready
make assay-evidence-ready
make alphabet-policy-ready
make citation-metadata-ready
make native-comparator-scope-ready
make distribution-record-ready
make bioconda-recipe-ready
make coverage
make public-crispr-evidence-gate
make crispr-comparison-gate
make barcode-comparison-gate
make feature-barcode-public-gate
make perturb-seq-public-gate
make amplicon-panel-public-gate
make bcl-tiny-public-gate
make oligo-adapter-public-gate
make workflow-examples-ready
npm run lint
npm audit --audit-level=moderate
npm run build
NEXT_OUTPUT=export NEXT_PUBLIC_BASE_PATH=/dotmatch NEXT_PUBLIC_SITE_URL=https://dnncha.github.io/dotmatch npm run build
```

`make bcl-comparison-gate` requires additional real-data and comparator evidence. Keep release notes within the evidence that is checked into the repository.

## Tagging

Use annotated tags:

```bash
git tag -a v0.1.0 -m "DotMatch v0.1.0"
git push origin v0.1.0
```

Pushing `v*` tags runs `.github/workflows/release.yml`. The workflow starts
with a preflight job that runs `make test`, `make cli-test`,
`make python-test`, `make repository-ready`, `make release-ready`, and
`make python-package-test`; artifact publication jobs depend on that preflight.
The workflow builds:

- raw Linux wheel release artifact;
- macOS wheel;
- source distribution;
- repaired manylinux/musllinux Linux wheels for PyPI;
- GHCR container image;
- `SHA256SUMS.txt`;
- PyPI publication through trusted publishing for the sdist, macOS wheel, and repaired Linux wheels;
- a draft GitHub release with generated notes.

Keep the GitHub release as a draft until the release notes, artifacts, checksums, `CITATION.cff`, and `codemeta.json` have been checked.

## Release Notes

Lead with:

- exact known-target short-DNA assignment;
- deterministic `unique`, `ambiguous`, `none`, and `invalid` semantics;
- CRISPR guide-counting, exact-prefix inline-barcode, narrow feature-barcode assignment, narrow CRISPR guide-capture assignment, and narrow ARTIC amplicon primer-start assignment evidence only where gates pass;
- package/install improvements;
- clear scope boundaries.

Avoid:

- genome-aligner language;
- universal guide-counter replacement language;
- broad barcode, feature quantification, amplicon consensus/variant-calling, or BCL comparative wording without matching evidence.

## Distribution Follow-Up

- Create a Zenodo archive and add the DOI to `CITATION.cff` when available.
- Publish the PyPI source distribution, native macOS wheel, and repaired manylinux/musllinux wheels through trusted publishing; do not upload raw `linux_x86_64` wheels.
- Track the Bioconda recipe PR after `make bioconda-recipe-ready`; for v0.1.0,
  [bioconda/bioconda-recipes#65367](https://github.com/bioconda/bioconda-recipes/pull/65367)
  has passed CI and is waiting for review/merge.
- Confirm the GHCR image labels and tag after the source tag is immutable.
- Run `make distribution-channels` after PyPI, Bioconda, GHCR, and Zenodo are public.
- Update `docs/distribution-release.json` with verified public and evidence links after public channels are live.
- Update `docs/scientific-claims.md` only when new evidence is committed and a corresponding gate passes.
