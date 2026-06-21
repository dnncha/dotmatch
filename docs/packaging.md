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

DotMatch 0.1.8 is published on PyPI; the `v0.1.8` release workflow publishes the source distribution, the native macOS wheel, and repaired manylinux/musllinux Linux wheels. The release workflow
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
merged the DotMatch 0.1.8 update on 2026-06-17. Anaconda package metadata and a
clean install smoke test now verify DotMatch 0.1.8 on `linux-64`, `osx-64`, and
`osx-arm64`, using immutable `v0.1.8` release sources. Treat future Bioconda
versions as available only after
`https://anaconda.org/bioconda/dotmatch`, repodata, and the install smoke tests
in `make distribution-channels` all verify the release version.

A release recipe template is kept under `packaging/bioconda/`. Before copying it
to `bioconda-recipes`, replace `REPLACE_WITH_RELEASE_TARBALL_SHA256` with the
SHA256 for the tagged GitHub release tarball. For 0.1.8 that SHA256 is
`ec3819bc773431454910287559d0809aca6ec1d81959f29d3d522650edb74904`. The checked-in
`docs/distribution-release.json` records the current channel state for the
active Bioconda handoff. Run `make bioconda-recipe-ready` before that copy so the
checked-in template stays aligned with the release version, native install
steps, CLI smoke tests, and scope notes.

After Bioconda merges the recipe, verify the channel with:

```bash
conda search -c bioconda dotmatch
make distribution-channels
```

The template also includes `extra.additional-platforms: [osx-arm64]` so the
Bioconda update opts into Apple Silicon CI/build coverage. Keep that selector in
future upstream recipe updates unless Bioconda CI demonstrates a
platform-specific blocker and the release notes clearly document that
`osx-arm64` is unavailable.

The current Bioconda recipe installs the Python `dotmatch` console script as
the user-facing command, with the native executable bundled inside the Python
package as `dotmatch-native`. It also installs the public C header, static
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

### Bioconda 0.1.8 PR changelog draft

- Update DotMatch from the latest accepted Bioconda version to 0.1.8.
- Use the immutable v0.1.8 tag and replace the SHA256 after the release tarball
  is available.
- Keep the Python console-script package scope introduced in 0.1.4: `dotmatch`
  exposes the native commands plus `assay`, `barcode`, `panel`, and
  GuideCounter-compatible CRISPR counting namespaces.
- Add Hamming `k=2`/`k=3` guide-counting support and exact audit safety fields;
  keep larger-radius claims bounded to same-length Hamming fixed-window
  assignment.
- Add native exact-table shortcuts, indexed status paths, and bounded k=2
  single-unknown status stops for the CRISPR/counting hot paths.
- Add installed-package smoke tests for `dotmatch count --help`,
  `dotmatch crispr-count --help`, `dotmatch audit --help`, Hamming `k=2`
  CRISPR counting, exact Hamming `k=3` audit summaries, GuideCounter-compatible
  counts/extended-counts/stats files, barcode offset inference, and panel
  design.
- Opt into `osx-arm64` builds with `extra.additional-platforms`.
- Keep host `zlib` for FASTQ.gz/native linkage and let Conda export `libzlib`
  at runtime.

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
the normal release path. DotMatch 0.1.8 package metadata and clean Bioconda
install smoke tests pass, and `make distribution-channels` can discover a
matching `quay.io/biocontainers/dotmatch:0.1.8--<build>` tag. The remaining
local check is Docker-backed manifest/runtime verification:

```bash
python3 scripts/check_distribution_channels.py --version 0.1.8
docker pull quay.io/biocontainers/dotmatch:0.1.8--<build>
docker run --rm quay.io/biocontainers/dotmatch:0.1.8--<build> dotmatch dist ACGT AGGT
docker run --rm quay.io/biocontainers/dotmatch:0.1.8--<build> dotmatch leq 1 ACGT AGGT
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
through Zenodo metadata for DotMatch 0.1.8. Version DOI
`10.5281/zenodo.20541629` belongs to v0.1.7 and is retained only as explicit
version-specific provenance, not as the v0.1.8 DOI.
