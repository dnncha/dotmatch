# Packaging Notes

DotMatch should ship with three practical install paths:

- source build with `make && make shared`;
- Docker image for reproducible command-line use;
- Python package using the ctypes wrapper and bundled or discoverable native library.

## PyPI

Initial local/GitHub packaging builds the native C core into the wheel as `dotmatch/libdotmatch.{so,dylib}` for Linux and macOS. Wheels are platform-specific but Python-ABI-neutral (`py3-none-<platform>`) because the native library is loaded through `ctypes` rather than the Python C API. The ctypes loader still accepts:

- the bundled platform library in the wheel;
- `DOTMATCH_LIB=/path/to/libdotmatch.{so,dylib}` for source-tree and custom installs.

Use `make python-package-test` to build the wheel, inspect that it contains the native library, install it into a clean virtual environment, and verify `import dotmatch` without `DOTMATCH_LIB` or `PYTHONPATH`.
The same verifier also builds the sdist, confirms it contains `src/qdalign.c` and `include/qdalign.h`, and installs that sdist into a clean virtual environment.

For PyPI, upload the sdist plus the native macOS wheel built on GitHub Actions. Linux binary wheels should go to PyPI only after they are built or repaired as manylinux/musllinux wheels. The release workflow builds repaired Linux wheel artifacts with cibuildwheel for `manylinux_x86_64` and `musllinux_x86_64`, smoke-tests `import dotmatch`, the installed console script, and `dotmatch dist ACGT AGGT`, and uploads them as GitHub release artifacts. Do not upload a raw `linux_x86_64` wheel to PyPI.

DotMatch 0.2.0 is the current release target; the `v0.2.0` release workflow publishes the source distribution, the native macOS wheel, and repaired manylinux/musllinux Linux wheels. The release workflow
uses PyPI trusted publishing from repository `dnncha/dotmatch`, workflow
`.github/workflows/release.yml`, and environment `pypi`; if that publisher is
missing or mismatched, the build artifacts are created but the publish job fails
with `invalid-publisher`.
Raw `linux_x86_64` wheels remain GitHub release artifacts only and are not uploaded to PyPI.
`make citation-metadata-ready` also checks PyPI-facing `pyproject.toml`
description, keywords, classifiers, and project URLs so the package page stays
discoverable for bioinformatics, CRISPR, FASTQ, barcode, and known-target
assignment searches.

## Bioconda

Bioconda packages DotMatch from a recipe in `bioconda-recipes`; DotMatch does
not upload a Conda package directly.
[bioconda/bioconda-recipes#65367](https://github.com/bioconda/bioconda-recipes/pull/65367)
published DotMatch 0.1.2 as the first Bioconda package.
[bioconda/bioconda-recipes#66291](https://github.com/bioconda/bioconda-recipes/pull/66291)
merged the DotMatch 0.1.8 update on 2026-06-17. The public Anaconda page currently shows 0.1.9. AssayCode is targeted for the
prepared 0.2.0 update, which still needs an immutable source tag, checksum,
recipe review, repodata visibility, and clean-install evidence. Treat the Bioconda
versions as available only after
`https://anaconda.org/bioconda/dotmatch`, repodata, and the install smoke tests
in `make distribution-channels` all verify the release version.

A release recipe template is kept under `packaging/bioconda/`. Before copying it
to `bioconda-recipes`, replace `REPLACE_WITH_RELEASE_TARBALL_SHA256` with the
SHA256 for the tagged GitHub release tarball. The checked-in
`docs/distribution-release.json` records the current channel state for the
active Bioconda handoff. Run `make bioconda-recipe-ready` before that copy so the
checked-in template stays aligned with the release version, native install
steps, CLI smoke tests, and scope notes.

The takeover-oriented distribution model uses two Bioconda coordinates without
renaming the existing package:

- `dotmatch` remains the engine, native library, compatibility CLI, and package
  upgraded by existing environments;
- `assaycode` is a noarch metapackage pinned to the matching `dotmatch` release,
  making `conda install -c bioconda assaycode` the future flagship install path
  without duplicating files or breaking existing dependencies.

The AssayCode metapackage template is
`packaging/bioconda/assaycode-meta.yaml`. Once `v0.2.0` exists, render both
upstream recipe directories from the downloaded immutable tag archive:

```bash
python scripts/prepare_bioconda_handoff.py \
  --release-tarball /path/to/v0.2.0.tar.gz \
  --out ./dotmatch-bioconda-handoff
```

Copy the resulting `recipes/dotmatch/` and `recipes/assaycode/` directories to
the `dnncha/bioconda-recipes` branch and submit them together. The renderer
computes and inserts the real DotMatch source checksum; the AssayCode recipe is
a dependency-only metapackage and therefore has no duplicate source payload.

After Bioconda merges the recipe, verify the channel with:

```bash
conda search -c bioconda dotmatch
conda search -c bioconda assaycode
make distribution-channels
```

The template also includes `extra.additional-platforms: [osx-arm64]` so the
Bioconda update opts into Apple Silicon CI/build coverage. Keep that selector in
future upstream recipe updates unless Bioconda CI demonstrates a
platform-specific blocker and the release notes clearly document that
`osx-arm64` is unavailable.

The prepared 0.2.0 Bioconda recipe installs both the additive `assaycode`
platform command and the compatibility-stable `dotmatch` command, with the
native executable bundled inside the Python package as `dotmatch-native`. It also installs the public C header, static
library, shared library, and license. Workbench and browser assets remain
outside the Bioconda recipe.

The recipe needs:

- `make`;
- `{{ compiler('c') }}` and `{{ stdlib('c') }}`;
- host `python`, `pip`, `setuptools`, `wheel`, and `zlib`;
- run `python`, plus `tomli` for Python versions before 3.11. Do not duplicate
  `zlib` in `run`: host `zlib` exports the linked `libzlib` runtime package;
- `run_exports` because the package installs a header and shared library;
- runtime tests for `dotmatch --version`, `dotmatch dist ACGT AGGT`,
  `dotmatch leq 1 ACGT AGGT`, Python import/native discovery, installed C
  artifacts, namespace help for `dotmatch assay`, `dotmatch barcode`, and
  `dotmatch panel`, tiny installed-package workflow smoke tests, and a
  GuideCounter-compatible `dotmatch guide-counter count` smoke test that writes
  counts, extended counts, and stats outputs.

The native CLI exposes `dotmatch --version`, so the Bioconda recipe and
post-release Bioconda install verifier should check version output as well as
functional CLI smoke tests.

### Bioconda 0.2.0 PR changelog draft

- Update the existing `dotmatch` package from the published 0.1.9 build to
  0.2.0; do not create or rename to a second Conda package.
- Use the immutable v0.2.0 tag and replace the SHA256 placeholder only after the
  release tarball exists.
- Install and smoke-test both `dotmatch` and `assaycode`; verify that the
  AssayCode Python namespace exposes the same DotMatch engine and version.
- Include AssayScript v2 compilation, experimental calibrated/joint decoding,
  bounded-memory sequential QC, the rewritten paper, and explicit experimental
  claim boundaries.
- Preserve the native commands, C header/static/shared libraries, workflow
  namespaces, GuideCounter compatibility, Hamming k=2/k=3 audit tests,
  `osx-arm64` opt-in, and host-zlib linkage.

## Docker

The root `Dockerfile` builds the native CLI and shared library on Debian. Example:

```bash
docker build -t dotmatch:dev .
docker run --rm dotmatch:dev --help
```

The image carries OCI labels for title, description, source, documentation,
version, license, and authorship. The release workflow smoke-tests both CLI
behavior and the `org.opencontainers.image.version` label before pushing tagged
images to `ghcr.io/dnncha/dotmatch`.

## BioContainers

BioContainers images for DotMatch are generated from the accepted Bioconda
recipe; there is no separate DotMatch Dockerfile to submit to BioContainers for
the normal release path. The 0.2.0 image is expected only after the Bioconda
recipe is accepted and propagated. The remaining local check is Docker-backed
manifest/runtime verification:

```bash
python3 scripts/check_distribution_channels.py --version 0.2.0
docker pull quay.io/biocontainers/dotmatch:0.2.0--<build>
docker run --rm quay.io/biocontainers/dotmatch:0.2.0--<build> dotmatch dist ACGT AGGT
docker run --rm quay.io/biocontainers/dotmatch:0.2.0--<build> dotmatch leq 1 ACGT AGGT
```

Do not publish a manual BioContainers image for DotMatch unless the Bioconda
automation fails after the accepted recipe is visible in Anaconda metadata and
the failure is documented in the release record.

## Post-Release Channel Verification

The prepared channel state is recorded in `docs/distribution-release.json`.
Check the package-channel record and recipe before tagging with:

```bash
make distribution-record-ready
make bioconda-recipe-ready
```

Before any public channel is verified, this record must stay in `not_released`
status with blockers and next actions. After some channels pass
`make distribution-channels`, use `partially_verified` and keep blockers only
for the remaining channels. Use `released` only when every required channel has
public evidence and the full post-release gate passes. For Bioconda, document
the exact platforms visible in repodata, including whether `osx-arm64`
propagated from the Apple Silicon recipe opt-in. Do not imply `linux-aarch64` or
any other platform availability unless those Bioconda subdirs contain DotMatch
for the release.

After publishing a tag, run:

```bash
make distribution-channels
```

This checks that the release version is visible on PyPI as a source distribution plus a macOS wheel and repaired manylinux/musllinux wheels, rejects raw
`linux_x86_64` PyPI wheels, installs with `pip install dotmatch==<version>` in a clean virtual environment, imports the Python package, runs the installed
`dotmatch` CLI, is available in Bioconda metadata, installs with
`conda create -p <env> -c conda-forge -c bioconda dotmatch=<version>` or
`micromamba`, runs the Bioconda `--version`/CLI and GuideCounter-compatible smoke tests, has a matching BioContainers
tag such as `quay.io/biocontainers/dotmatch:<version>--<build>` that runs CLI
distance and threshold smoke tests, is published as
`ghcr.io/dnncha/dotmatch:vX.Y.Z`, runs with
`docker run --rm ghcr.io/dnncha/dotmatch:v<version> --version` and a CLI distance
smoke test, and is backed by a DOI in `CITATION.cff` that resolves through
`doi.org`, and reports the same release version from Zenodo record metadata. It
is not part of `make release-ready` because it should fail until public
publication has actually happened.

## Zenodo

The repository includes `.zenodo.json` metadata for tagged software archives.
General software citation uses DOI `10.5281/zenodo.20541628`, which resolves
through Zenodo metadata for DotMatch. The 0.2.0 archive metadata remains a
post-tag verification step. Version DOI
`10.5281/zenodo.20541629` belongs to v0.1.7 and is retained only as explicit
version-specific provenance, not as the v0.2.0 DOI.
