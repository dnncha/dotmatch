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

For PyPI, upload the sdist first. Linux binary wheels should go to PyPI only
after they are built or repaired as manylinux/musllinux wheels. The release
workflow builds repaired Linux wheel artifacts with cibuildwheel for
`manylinux_x86_64` and `musllinux_x86_64`, smoke-tests `import dotmatch`, the
installed console script, and `dotmatch dist ACGT AGGT`, and uploads them as
GitHub release artifacts. Do not upload a raw `linux_x86_64` wheel to PyPI.

The release workflow is prepared for PyPI trusted publishing and publishes the source distribution and repaired manylinux/musllinux wheels on tagged releases.
Raw `linux_x86_64` wheels remain GitHub release artifacts only and are not
uploaded to PyPI.
`make citation-metadata-ready` also checks PyPI-facing `pyproject.toml`
description, keywords, classifiers, and project URLs so the package page stays
discoverable for bioinformatics, CRISPR, FASTQ, barcode, and known-target
assignment searches.

## Bioconda

Bioconda should package the native CLI first. A release recipe template is kept
under `packaging/bioconda/`. Before copying it to `bioconda-recipes`,
replace `REPLACE_WITH_RELEASE_TARBALL_SHA256` with the SHA256 for the tagged
GitHub release tarball. Run `make bioconda-recipe-ready` before that copy so the
checked-in template stays aligned with the release version, native install
steps, CLI smoke tests, and no-genome-aligner scope wording.

The recipe needs:

- `make`
- C compiler and C++ compiler for optional benchmark tools;
- zlib;
- runtime test: `dotmatch dist ACGT AGGT`.

## Docker

The root `Dockerfile` builds the native CLI and shared library on Debian. Example:

```bash
docker build -t dotmatch:dev .
docker run --rm -v "$PWD:/work" dotmatch:dev count --help
```

The image carries OCI labels for title, description, source, documentation,
version, license, and authorship. The release workflow smoke-tests both CLI
behavior and the `org.opencontainers.image.version` label before pushing tagged
images to `ghcr.io/dnncha/dotmatch`.

## Post-Release Channel Verification

The prepared channel state is recorded in `docs/distribution-release.json`.
Check the package-channel record and recipe before tagging with:

```bash
make distribution-record-ready
make bioconda-recipe-ready
```

While the first public release is still pending, this record must stay in
`not_released` status with blockers and next actions for every public channel.
After publication, replace the expected links with verified public and evidence
URLs, set channels to `verified`, and run the post-release gate.

After publishing a tag, run:

```bash
make distribution-channels
```

This checks that the release version is visible on PyPI as a source distribution plus repaired manylinux/musllinux wheels, rejects raw
`linux_x86_64` PyPI wheels, installs with `pip install dotmatch==<version>` in a clean virtual environment, imports the Python package, runs the installed
`dotmatch` CLI, is available in Bioconda metadata, installs with
`conda create -p <env> -c conda-forge -c bioconda dotmatch=<version>` or
`micromamba`, runs the Bioconda CLI smoke tests, has a matching BioContainers
tag such as `quay.io/biocontainers/dotmatch:<version>--<build>` that runs
`dotmatch --version` and a CLI distance smoke test, is published as
`ghcr.io/dnncha/dotmatch:vX.Y.Z`, runs with
`docker run --rm ghcr.io/dnncha/dotmatch:v<version> --version` and a CLI distance
smoke test, and is backed by a DOI in `CITATION.cff` that resolves through
`doi.org`. It is not part of `make release-ready` because it should fail until
public publication has actually happened.

## Zenodo

The repository includes `.zenodo.json` metadata for tagged software archives.
Do not add a DOI to `CITATION.cff` until Zenodo has minted one for an immutable
release archive.
