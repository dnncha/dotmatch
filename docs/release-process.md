# Release Process

DotMatch releases should be specific, reproducible, and tied to checked evidence.

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
make asan
make python-test
make python-package-test
make docs-ready
make repository-ready
make release-ready
make scientific-readiness-ready
make assay-evidence-ready
make alphabet-policy-ready
make citation-metadata-ready
make native-comparator-scope-ready
make distribution-record-ready
make bioconda-recipe-ready
make coverage
make native-exact-gate
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
git tag -a v<version> -m "DotMatch v<version>"
git push origin v<version>
```

Pushing `v*` tags runs `.github/workflows/release.yml`. The workflow starts
with a preflight job that runs `make test`, `make cli-test`, `make asan`,
`make python-test`, installs the public docs toolchain, `make repository-ready`,
`make release-ready`, and `make python-package-test`; artifact publication jobs
depend on that preflight.
The workflow builds:

- raw Linux wheel release artifact;
- macOS wheel;
- source distribution;
- repaired manylinux/musllinux Linux wheels for `x86_64` and `aarch64` for PyPI;
- GHCR container image index for `linux/amd64` and `linux/arm64`;
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
- broad barcode, feature quantification, amplicon consensus/variant-calling, or BCL comparisons without matching evidence.

## Distribution Follow-Up

- Confirm the Zenodo archive for the tagged release and add the release DOI to
  `CITATION.cff` when available.
- Publish the PyPI source distribution, native macOS wheel, and repaired manylinux/musllinux wheels for `x86_64` and `aarch64` through trusted publishing; do not upload raw Linux wheels. The PyPI project must have a trusted publisher matching repository `dnncha/dotmatch`, workflow `.github/workflows/release.yml`, and environment `pypi`.
- For Bioconda updates, submit or update the `bioconda-recipes` recipe after
  `make bioconda-recipe-ready`. Keep the `osx-arm64` additional-platforms opt-in
  in that recipe copy so Bioconda CI validates the Apple Silicon build. Replace
  the source SHA256 only in the upstream recipe copy. After merge and channel
  propagation, verify with
  `make distribution-channels` before announcing conda install instructions or
  BioContainers availability.
- Confirm the GHCR image labels, tag, and `linux/amd64` plus `linux/arm64` manifest descriptors after the source tag is immutable.
- Run `make distribution-channels` after PyPI, Bioconda, GHCR, and Zenodo are public.
- Update `docs/distribution-release.json` with verified public and evidence links, exact PyPI Linux wheel architectures, and exact GHCR platforms after public channels are live.
- Update `docs/scientific-claims.md` only when new evidence is committed and a corresponding gate passes.

## Current 0.7.0 handoff

The user authorized all configured publication channels. DotMatch 0.7.0 is published on GitHub, PyPI and GHCR; the previous 0.6.4 distribution record is preserved in
`releases/v0.6.4-distribution.json`. The current record verifies PyPI and GHCR against the completed release
workflow; the version-specific Zenodo archive is also verified. Bioconda and
BioContainers remain independently unverified.

The approved-release workflow tagged checked commit
`5736f9dd14fde44f3d89fa99966802af344bb0a0`; publication and public artifact
verification passed after retrying an ARM64 build-image pull. GitHub API requests are blocked from this environment; git pushes and
public workflow pages remain accessible. Hosted acceptance checks must pass
before the workflow tags and publishes the exact main commit. Do not manually
tag the local work branch.

After the public source archive exists, download the archive referenced by the
recipe and run `scripts/prepare_bioconda_handoff.py --release-tarball
release/v0.7.0.tar.gz --out bioconda-handoff`. Submit the generated DotMatch
and AssayCode recipes to bioconda-recipes. PR #69711 merged the older DotMatch 0.6.4 recipe on September 30, 2026.
Check for a new BiocondaBot update before opening a duplicate 0.7.0 PR; neither
the older merge nor its generated image tags verifies the new release.
Verify Anaconda propagation and the generated BioContainers image separately.

The experimental `crisprworks-fit` package has scientific CI and build artifacts
but no configured trusted PyPI publication workflow. Publishing DotMatch does
not publish that separate package; its alpha status and outstanding scientific
validation must remain explicit when adding a publisher.

The public GitHub source archive for 0.7.0 has SHA256
`c1b2bae12b56738c76d0b24d095c1c3c1deb42242c4553913187122de0575bd1`.
Both recipes are prepared, preserving the upstream platform declarations and
AssayCode run exports. Existing PR #68663 still targets 0.4.0; update it instead
of opening a duplicate once fork write access and GitHub API access are restored.
The working environment received HTTP 401 when pushing the recipe fork, so no
0.7.0 upstream submission or channel propagation is claimed.
